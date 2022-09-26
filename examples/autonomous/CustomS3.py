#!/usr/bin/python3
import os
import time, datetime
from pathlib import Path

import numpy as np

import io
# minio (S3 API) Python client reference:
# https://docs.min.io/docs/python-client-api-reference.html
from minio import Minio
from minio.error import S3Error

from Base import Base

class CustomS3(Base):
    
    def __init__(self, 
                 username,                              # S3 access_key
                 send='send',                           # Queue to send/publish to
                 receive='receive',                     # Queue to watch for signals
                 endpoint="localhost:8000",             # S3/lustre server
                 secret_key=None,                       # Authentication key
                 bucket_name='experiments',             # S3 bucket to use
                 experiment=None,                       # Name in S3 storage for the datasets
                 name='S3',                             # Name of this object (for print/log purposes)
                 save_dir='./',                         # Location for tmp save files
                 log_verbosity=5,                       # If a 'common' object is defined
                 **kwargs                               # Additional keywork arguments
                 ):
        
        super().__init__(name=name, log_verbosity=log_verbosity, **kwargs) # class Base()
        
        if secret_key is None:
            # Search in 'typical' paths
            p = Path('./S3_secret_key.txt') # Local to script execution
            if not p.exists():
                p = Path.home() / ".secret/S3_secret_key.txt" # User's home dir
        else:
            p = Path(secret_key)
        
        if not p.exists():
            self.msg_error('Specified secret_key path is not valid.', 1, 0)
            
        with open(p) as fin:
            secret_key = fin.readline().strip()
            
            
        self.client = Minio(endpoint, access_key=username, secret_key=secret_key)
        
        self.bucket_name = bucket_name
        if experiment is None:
            experiment = 'experiment_{}'.format(self.now(str_format='%Y-%m-%d_%H'))
        self.experiment = experiment
        self.send = send
        self.receive = receive
        
        self.save_dir = save_dir
            
        

    def send_path(self):
        return '{}/{}'.format(self.experiment, self.send)
    def receive_path(self):
        return '{}/{}'.format(self.experiment, self.receive)
    
    
    def getS3_floats(self):
        
        import json
        
        events = self.client.listen_bucket_notification(
            bucket_name=self.bucket_name,
            prefix='{}/'.format(self.receive_path()),
            events=["s3:ObjectCreated:Put"],
        )
        
        
        for event in events:
            #self.print_d(event)
            assert len(event["Records"]) == 1

            record = event['Records'][0]
            #self.print_d(record)
            content_type = record['s3']['object']['contentType']
            assert content_type == "AE/custom"
            
            object_name = record['s3']['object']['key']
            timestamp = record["s3"]["object"]["userMetadata"]["X-Amz-Meta-Timestamp"]
            shape = record['s3']['object']['userMetadata']['X-Amz-Meta-Shape']
            shape = tuple(json.loads(f'[{shape}]'))
            self.msg('Received S3 data: {}'.format(object_name), 4, 2)
            
            sucessful = False
            try:
                data_stream = self.client.get_object(self.bucket_name, object_name)
                data_bytes = data_stream.data
                data_np = np.frombuffer(data_bytes).reshape(shape)
                sucessful = True
                
            except Exception as ex:
                self.msg_error('Python Exception in getS3', 1, 2)
                self.print(ex)
                
            finally:
                data_stream.close()
                data_stream.release_conn()

            if not sucessful:
                self.msg_warning('S3 data retrieval failed.', 2, 2)
            
            break # Process just that one event
        
        return data_np
            
        
    
    def publishS3_floats(self, data):
        
        now_str = self.now(str_format='%Y%m%d_%H%M%S_%f')
        object_name = '{}/obj{}'.format(self.send_path(), now_str)
        
        data_bytes = data.tobytes()
        data_stream = io.BytesIO(data_bytes)
        
        
        result = self.client.put_object(
            bucket_name=self.bucket_name,
            object_name=object_name,
            data=data_stream,
            length=len(data_bytes),
            content_type="AE/custom",
            metadata={
                'timestamp': self.now(),
                'shape': data.shape,
                }
            )

        self.msg('Sent S3 data: {}'.format(object_name), 4, 2)


    def getS3_file(self):
        
        events = self.client.listen_bucket_notification(
            bucket_name=self.bucket_name,
            prefix='{}/'.format(self.receive_path()),
            events=["s3:ObjectCreated:Put"],
        )
        
        
        file_path = '{}/{}-received.npy'.format(self.save_dir, self.name)
        
        
        for event in events:
            #self.print_d(event)
            assert len(event["Records"]) == 1

            record = event['Records'][0]
            #self.print_d(record)
            content_type = record['s3']['object']['contentType']
            assert content_type == "AE/custom"
            
            object_name = record['s3']['object']['key']
            timestamp = record["s3"]["object"]["userMetadata"]["X-Amz-Meta-Timestamp"]
            self.msg('Received S3 data: {}'.format(object_name), 4, 2)
            
            sucessful = False
            try:
                self.client.fget_object(
                    self.bucket_name, 
                    object_name,
                    file_path
                    )
                data = np.load(file_path, allow_pickle=True)
                sucessful = True
                
            except Exception as ex:
                self.msg_error('Python Exception in getS3', 1, 2)
                self.print(ex)
                

            if not sucessful:
                self.msg_warning('S3 data retrieval failed.', 2, 2)
            else:
                np.save(file_path, data, allow_pickle=True)
            
            break # Process just that one event
        
        
        return data
            
        
    
    def publishS3_file(self, data):
        
        file_path = '{}/{}-sent.npy'.format(self.save_dir, self.name)
        np.save(file_path, data, allow_pickle=True)
        
        now_str = self.now(str_format='%Y%m%d_%H%M%S_%f')
        object_name = '{}/obj{}'.format(self.send_path(), now_str)
        
        
        result = self.client.fput_object(
            bucket_name=self.bucket_name,
            object_name=object_name,
            file_path=file_path,
            content_type="AE/custom",
            metadata={
                'timestamp': self.now(),
                }
            )

        self.msg('Sent S3 data: {}'.format(object_name), 4, 2)


    def publish_status_file(self, file_path, name=None):
        self.msg('Uploading status file ({})...'.format(self.now()), 4, 1)
        
        p = Path(file_path)
        name = p.name if name is None else '{}{}'.format(name, p.suffix)
        object_name = '{}/status/{}/{}'.format(self.experiment, self.send, name)
        
        
        result = self.client.fput_object(
            bucket_name=self.bucket_name,
            object_name=object_name,
            file_path=file_path,
            content_type="AE/status",
            metadata={
                'timestamp': self.now(),
                'mtime': p.stat().st_mtime,
                'ctime': p.stat().st_ctime,
                'filesize': p.stat().st_size,
                }
            )

        self.msg('Sent S3 file: {}'.format(object_name), 4, 2)


    def get_status_files(self, name='status', timestamp=False):
        
        prefix = '{}/{}/'.format(self.experiment, name)
        now_str = self.now(str_format='%Y-%m-%d_%H%M%S')
        
        self.msg('Getting status files ({})'.format(self.now()), 4, 1)
        self.msg('recursive searching: {}'.format(prefix), 4, 2)
        
        objects = self.client.list_objects(
            bucket_name=self.bucket_name,
            prefix=prefix,
            recursive=True,
            )
        
        for obj in objects:
            self.msg('downloading: {}'.format(obj.object_name), 4, 3)
            
            if timestamp:
                obj_str = obj.object_name[len(self.experiment):][len('/status/'):]
                file_path = '{}/status/{}/{}'.format( self.save_dir, now_str, obj_str )
            else:
                file_path = '{}{}'.format( self.save_dir, obj.object_name[len(self.experiment):] )
            
            try:
                self.client.fget_object(
                    self.bucket_name, 
                    obj.object_name,
                    file_path
                    )

            except Exception as ex:
                self.msg_error('Python Exception in get_status_files', 1, 2)
                self.print(ex)
            
        self.msg('Done.', 4, 2)
    

    def get(self, save=True, check_interrupted=True, force_load=False):
        '''Get the current item being published.'''
        #message = self.from_socket.recv()
        
        self.msg('Getting data/command ({})'.format(self.now()), 4, 1)
        
        if force_load or (check_interrupted and self.interrupted()):
            self.msg('Loading data/command from disk...'.format(self.now()), 4, 1)
            data = self.load()
            
            if isinstance(data, (list, tuple, np.ndarray)):
                self.msg('Loaded: list length {}'.format(len(data)), 4, 2)
            else:
                self.msg('Loaded.', 4, 2)
        
        else:
            self.msg('Waiting for data/command ({})...'.format(self.now()), 4, 1)
            #data = self.getS3_floats()
            data = self.getS3_file()
            
            if isinstance(data, (list, tuple, np.ndarray)):
                self.msg('Received: list length {}'.format(len(data)), 4, 2)
            else:
                self.msg('Received.', 4, 2)
            
            #if save:
                #np.save('{}/{}-received.npy'.format(self.save_dir, self.name), data, allow_pickle=True)

        return data


    
    def publish(self, data, save=True):
        
        self.msg('Sending data/command ({})...'.format(self.now()), 4, 1)
        #self.publishS3_floats(data)
        self.publishS3_file(data)
        self.msg('Sent.', 4, 2)

        #if save:
            #np.save('{}/{}-sent.npy'.format(self.save_dir, self.name), data, allow_pickle=True)

    def republish(self, stype='sent'):
        '''Re-send the last data/command.
        You can use this to force the loop to restart, if you interrupted it.'''
        
        self.msg('Re-publishing the last data/command that was {} ({})'.format(stype, self.now()), 4, 1)
        self.interrupted()
        
        data = self.load(stype=stype)
        self.publish(data)
        
        
    def interrupted(self):
        
        received = '{}/{}-received.npy'.format(self.save_dir, self.name)
        
        if not os.path.exists(received):
            self.msg('No saved received data from prior work-cycle.', 4, 2)
            return False
        
        received = os.path.getmtime(received)
        sent = '{}/{}-sent.npy'.format(self.save_dir, self.name)
        if os.path.isfile(sent):
            sent = os.path.getmtime(sent)
        else:
            self.msg('No saved sent data from prior work-cycle.', 4, 2)
            return True
        
        if sent>received:
            # As expected, the file we sent is newer
            self.msg('Last work-cycle completed normally (based on file dates).', 4, 2)
            self.msg('Received: {}'.format(self.time_str(received)), 5, 3)
            self.msg('Sent: {} ({})'.format(self.time_str(sent), self.time_delta(received, sent)), 5, 3)
            return False
        else:
            self.msg('Last work-cycle INTERRUPTED (based on file dates).', 3, 2)
            self.msg('Received: {}'.format(self.time_str(received)), 3, 3)
            self.msg('Sent: {} ({})'.format(self.time_str(sent), self.time_delta(received, sent)), 3, 3)
            return True


    def load(self, stype='received'):
        data = np.load('{}/{}-{}.npy'.format(self.save_dir, self.name, stype), allow_pickle=True)
        return data


    def clear(self):
        '''Remove the saved filed.'''
        self.msg('Clearing queue saved files.', 2, 1)
        
        received = '{}/{}-received.npy'.format(self.save_dir, self.name)
        if os.path.exists(received):
            self.msg('Removing {}'.format(received), 3, 2)
            os.remove(received)
        else:
            self.msg('Received data does not exist ({})'.format(received), 3, 2)

        sent = '{}/{}-sent.npy'.format(self.save_dir, self.name)
        if os.path.exists(sent):
            self.msg('Removing {}'.format(sent), 3, 2)
            os.remove(sent)
        else:
            self.msg('Sent data does not exist ({})'.format(sent), 3, 2)





# Loop for Autonomous Experimentation (AE):
########################################

VERBOSITY=4
ENDPOINT="dtn01.sdcc.bnl.gov:8000"
USERNAME='username'
SECRET_KEY='/nsls2/xf11bm/data/2021_3/UShell/S3_secret_key.txt' # Specify path
SECRET_KEY=None # Use default path
EXPERIMENT='experiment_AE_current'
#EXPERIMENT=None # Default is to use day and hour


class Queue_decision(CustomS3): # gpCAM
    def __init__(self, username=USERNAME, send='decision', receive='analyze', endpoint=ENDPOINT, secret_key=SECRET_KEY, experiment=EXPERIMENT, name='decisionS3', save_dir='./', verbosity=VERBOSITY, **kwargs):
        super().__init__(username=username, send=send, receive=receive, endpoint=endpoint, secret_key=secret_key, experiment=experiment, name=name, save_dir=save_dir, verbosity=verbosity, **kwargs)

class Queue_measure(CustomS3): # beamline
    def __init__(self, username=USERNAME, send='measure', receive='decision', endpoint=ENDPOINT, secret_key=SECRET_KEY, experiment=EXPERIMENT, name='measureS3', save_dir='./', verbosity=VERBOSITY, **kwargs):
        super().__init__(username=username, send=send, receive=receive, endpoint=endpoint, secret_key=secret_key, experiment=experiment, name=name, save_dir=save_dir, verbosity=verbosity, **kwargs)

class Queue_analyze(CustomS3): # SciAnalysis
    def __init__(self, username=USERNAME, send='analyze', receive='measure', endpoint=ENDPOINT, secret_key=SECRET_KEY, experiment=EXPERIMENT, name='analyzeS3', save_dir='./', verbosity=VERBOSITY, **kwargs):
        super().__init__(username=username, send=send, receive=receive, endpoint=endpoint, secret_key=secret_key, experiment=experiment, name=name, save_dir=save_dir, verbosity=verbosity, **kwargs)



########################################
# The usage would be:
########################################

# Inside gpCAM:
########################################
#from CustomS3 import Queue_decision
#q = Queue_decision()

#if first_iteration:
    ## Send an initial signal to the beamline
    #data = np.asarray([0])
    #q.publish(data)

#while True: # The loop that waits for new instructions...
    
    #data = q.get() # Get analysis results
    
    ## Do whatever work needs to be done
    ## time.sleep(1)
    
    #q.publish(data) # Send new command to beamline
    


# Inside beamline user.py:
########################################
#try:
    #measure_queue
#except NameError:
    #from CustomS3 import Queue_measure
    #measure_queue = Queue_measure()

#while True: # The loop that waits for new instructions...
    
    #data = measure_queue.get() # Get measurement command
    
    ## Do whatever work needs to be done
    ## time.sleep(1)
    
    #measure_queue.publish(data) # Send new results for analysis
    


# Inside SciAnalysis autonomous.py:
########################################
#from CustomS3 import Queue_analyze
#q = Queue_analyze()

##q.republish() # Restart if necessary

#while True: # The loop that waits for new instructions...
    
    #data = q.get() # Get analysis command
    
    ## Do whatever work needs to be done
    ## time.sleep(1)
    
    #q.publish(data) # Send new analysis results to gpCAM


