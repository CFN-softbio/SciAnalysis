#!/usr/bin/python3

import os
import time, datetime
import zmq
import numpy as np



class CustomQueue():
    
    def __init__(self, from_port, to_port, from_ip='localhost', to_ip='*', name='', save_dir='./', verbosity=4, **kwargs):

        # Save these in case we want to check them later
        self.from_ip = from_ip
        self.from_port = from_port
        self.to_ip = to_ip
        self.to_port = to_port
        
        self.name = name
        self.save_dir = save_dir
        self.verbosity = verbosity
        self.kwargs = kwargs
        
        self.context = zmq.Context()
        
        self.msg('Starting queues...', 2)


        # To/Tell/Publish/Push/Producer
        # Act as a "server", sending data when requested
        self.msg('Starting push queue (server)', 3, 1)
        self.msg('IP: {}, port: {}'.format(to_ip, to_port), 3, 2)
        success, attempt, start_time = False, 0, time.time()
        while not success:
            attempt += 1
            try:
                self.to_socket = self.context.socket(zmq.PUSH)
                self.to_socket.bind("tcp://{}:{:d}".format(to_ip, to_port))
                success = True
                self.msg('Success (attempt {:d}, {:.1f} s)'.format(attempt, time.time()-start_time), 3, 2)
            except zmq.error.ZMQError:
                self.msg('Retrying (attempt {:d}, {:.1f} s)'.format(attempt, time.time()-start_time), 3, 2)
                time.sleep(1)
        
        
        # From/Ask/Subscribe/Pull/Consumer
        # Act as a "client", asking for data 
        self.msg('Starting pull queue (client)', 3, 1)
        self.msg('IP: {}, port: {}'.format(from_ip, from_port), 3, 2)
        success, attempt, start_time = False, 0, time.time()
        while not success:
            attempt += 1
            try:
                self.from_socket = self.context.socket(zmq.PULL)
                self.from_socket.connect("tcp://{}:{:d}".format(from_ip, from_port))
                success = True
                self.msg('Connected to pull server (attempt {:d}, {:.1f} s)'.format(attempt, time.time()-start_time), 3, 2)
            except zmq.error.ZMQError:
                self.msg('Waiting for pull server (attempt {:d}, {:.1f} s)'.format(attempt, time.time()-start_time), 3, 2)
                time.sleep(1)
                
        
        self.msg('Started.', 3)


    def msg(self, txt, threshold=3, indent=0, indent_txt='  '):
        if self.verbosity>=threshold:
            indent = np.clip(indent, 0, 10)
            indent = indent_txt*indent
            print('ZMQ: {}> {}{}'.format(self.name, indent, txt))

    def now(self, str_format='%Y-%m-%d %H:%M:%S'):
        return time.strftime(str_format)

    def time_str(self, timestamp, str_format='%Y-%m-%d %H:%M:%S'):
        #time_tuple = time.gmtime(timestamp)
        #s = time.strftime(str_format, time_tuple)
        s = datetime.datetime.fromtimestamp(timestamp).strftime(str_format)
        return s

    def time_delta(self, first, second):
        diff = abs(second-first)
        if diff>60*60:
            diff = '{:.1f} hours'.format(diff/(60*60))
        elif diff>60:
            diff = '{:.1f} minutes'.format(diff/60)
        elif diff>0.1:
            diff = '{:.1f} s'.format(diff)
        else:
            diff = '{:.4f} s'.format(diff)
        
        if second<first:
            return '{} earlier'.format(diff)
        else:
            return '{} later'.format(diff)
        

    def reply(self):
        # Wait for next request from client
        message = self.to_socket.recv()
        self.msg('Received request: {}'.format(message))
        
        # Do some 'work'
        time.sleep(1)
        
        # Send reply back
        self.to_socket.send("REPLY")
        
    def ask(self):
        self.msg('Sending request', 3, 1)
        self.from_socket.send(b"ASK")
        message = self.from_socket.recv()
        self.msg('Received reply: {}'.format(message))


    def get(self, save=True, check_interrupted=False, force_load=False):
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
            data = self.from_socket.recv_pyobj()
            
            if isinstance(data, (list, tuple, np.ndarray)):
                self.msg('Received: list length {}'.format(len(data)), 4, 2)
            else:
                self.msg('Received.', 4, 2)
            
            if save:
                np.save('{}/{}-received.npy'.format(self.save_dir, self.name), data, allow_pickle=True)

        return data


    def publish(self, data, save=True):
        #message = '{} {}'.format(self.name, data)
        #self.to_socket.send(message)
        
        self.msg('Sending data/command ({})...'.format(self.now()), 4, 1)
        self.to_socket.send_pyobj(data)
        self.msg('Sent.', 4, 2)

        if save:
            np.save('{}/{}-sent.npy'.format(self.save_dir, self.name), data, allow_pickle=True)
        
        
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
        

    def print_connections(self):
        txt = 'Push queue: IP: {}, port: {}\nPull queue: IP: {}, port: {}'.format(self.to_ip, self.to_port, self.from_ip, self.from_port)
        self.msg(txt, 1, 0)

        
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
        
        

VERBOSITY = 4

# Simple example of three systems working in a ring/loop:
########################################
class OneQueue(CustomQueue):
    def __init__(self, from_port=5553, to_port=5551, from_ip='localhost', to_ip='*', name='One', verbosity=VERBOSITY, **kwargs):
        super().__init__(from_port=from_port, to_port=to_port, from_ip=from_ip, to_ip=to_ip, name=name, verbosity=verbosity, **kwargs)
class TwoQueue(CustomQueue):
    def __init__(self, from_port=5551, to_port=5552, from_ip='localhost', to_ip='*', name='Two', verbosity=VERBOSITY, **kwargs):
        super().__init__(from_port=from_port, to_port=to_port, from_ip=from_ip, to_ip=to_ip, name=name, verbosity=verbosity, **kwargs)
class ThreeQueue(CustomQueue):
    def __init__(self, from_port=5552, to_port=5553, from_ip='localhost', to_ip='*', name='Three', verbosity=VERBOSITY, **kwargs):
        super().__init__(from_port=from_port, to_port=to_port, from_ip=from_ip, to_ip=to_ip, name=name, verbosity=verbosity, **kwargs)



# Example of arbitrary-sized ring:
########################################
if False:
    # Servers
    c = [
        {'name': 'decision', 'ip': '127.0.0.1', 'port': 5551} ,
        {'name': 'measure', 'ip': '127.0.0.1', 'port': 5552} ,
        {'name': 'process', 'ip': '127.0.0.1', 'port': 5553} ,
        {'name': 'analyze', 'ip': '127.0.0.1', 'port': 5554} ,
        ]
    queues = []
    for i, server in enumerate(c):
        # Connect each process 'from' to the prior entry in the list, and their 'to' to the current entry.
        from_i = i-1 if i>0 else len(c)-1
        to_i = i
        arguments = {'name': c[i]['name'], 'from_ip': c[from_i]['ip'], 'from_port': c[from_i]['port'], 'to_ip': c[to_i]['ip'], 'to_port': c[to_i]['port']}
        queues.append(arguments)
        
    # Then each process can be started using:
    #from CustomQueue import *
    #q = CustomQueue(**queues[0]) # Each process uses a different index in the queues list
    #while True:
        #data = q.get()
        #q.publish(data)




# Loop for Autonomous Experimentation (AE):
########################################
local = '127.0.0.1'


xf11bm_ws1 = '10.68.80.221'
xf11bm_ws2 = '10.68.80.222'
xf11bm_ws3 = '10.68.80.223'
xf11bm_ws4 = '10.68.80.224'
xf11bm_ws5 = '10.68.80.225'
xf11bm_ws6 = '10.68.80.226'

xf11bm_srv1 = '10.68.80.25'
xf11bm_gpu1 = '10.68.80.27'
xf11bm_gpu2 = '10.68.80.29'


xf12id_ws1 = '10.12.0.201'
xf12id_ws2 = '10.12.0.202'
xf12id_ws3 = '10.12.0.203'
xf12id_srv1 = '10.12.0.4'
xf12id_srv2 = '10.12.0.5'


# Server connections
#c = {
    #'decision' : {'ip': xf11bm_srv1, 'port': 5551} ,
    #'measure' : {'ip': xf11bm_ws1, 'port': 5552} ,
    #'analyze' : {'ip': xf11bm_ws2, 'port': 5553} ,
    #}
c = {
    'decision' : {'ip': xf12id_srv1, 'port': 5551} ,
    'measure' : {'ip': xf12id_ws1, 'port': 5552} ,
    'analyze' : {'ip': xf12id_srv1, 'port': 5553} ,
    }


#c = {  'decision': {'ip': local, 'port': 5551},  'measure': {'ip': local, 'port': 5552},  'analyze': {'ip': local, 'port': 5553}  } # When all processes are on the same machine


class Queue_decision(CustomQueue): # gpCAM
    def __init__(self, from_port=c['analyze']['port'], to_port=c['decision']['port'], from_ip=c['analyze']['ip'], to_ip=c['decision']['ip'], name='decision', save_dir='./', verbosity=VERBOSITY, **kwargs):
        super().__init__(from_port=from_port, to_port=to_port, from_ip=from_ip, to_ip=to_ip, name=name, save_dir=save_dir, verbosity=verbosity, **kwargs)
        
    def get(self, save=True, check_interrupted=True, force_load=False):
        return super().get(save=save, check_interrupted=check_interrupted, force_load=force_load)

class Queue_measure(CustomQueue): # beamline
    def __init__(self, from_port=c['decision']['port'], to_port=c['measure']['port'], from_ip=c['decision']['ip'], to_ip=c['measure']['ip'], name='measure', save_dir='./', verbosity=VERBOSITY, **kwargs):
        super().__init__(from_port=from_port, to_port=to_port, from_ip=from_ip, to_ip=to_ip, name=name, save_dir=save_dir, verbosity=verbosity, **kwargs)

    def get(self, save=True, check_interrupted=True, force_load=False):
        return super().get(save=save, check_interrupted=check_interrupted, force_load=force_load)

class Queue_analyze(CustomQueue): # SciAnalysis
    def __init__(self, from_port=c['measure']['port'], to_port=c['analyze']['port'], from_ip=c['measure']['ip'], to_ip=c['analyze']['ip'], name='analyze', save_dir='./', verbosity=VERBOSITY, **kwargs):
        super().__init__(from_port=from_port, to_port=to_port, from_ip=from_ip, to_ip=to_ip, name=name, save_dir=save_dir, verbosity=verbosity, **kwargs)

    def get(self, save=True, check_interrupted=True, force_load=False):
        return super().get(save=save, check_interrupted=check_interrupted, force_load=force_load)

class Queue_analyzeFix(Queue_analyze): # SciAnalysis
    '''This version includes an ad-hoc fix to detect and stop
    any "double loops" running in the queue loop.'''

    def __init__(self, from_port=c['measure']['port'], to_port=c['analyze']['port'], from_ip=c['measure']['ip'], to_ip=c['analyze']['ip'], name='analyze', save_dir='./', verbosity=VERBOSITY, **kwargs):
        super().__init__(from_port=from_port, to_port=to_port, from_ip=from_ip, to_ip=to_ip, name=name, save_dir=save_dir, verbosity=verbosity, **kwargs)

        self.list_lengths = []

    def get(self, save=True, check_interrupted=True, force_load=False):
        
        data = super().get(save=save, check_interrupted=check_interrupted, force_load=force_load)
        
        self.list_lengths.append(len(data))
        if len(self.list_lengths)>=2:
            len_cur = self.list_lengths[-1]
            len_prev = self.list_lengths[-2]
            
            if len_cur==len_prev:
                self.msg("WARNING: Possible double-loop detected.", 3, 1)
                self.msg("Ignoring this command and waiting for next command instead.", 3, 2)
                self.get(save=save, check_interrupted=False, force_load=False)
        
        return data



########################################
# The usage would be:
########################################

# Inside gpCAM:
########################################
#from CustomQueue import Queue_decision
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
    #from CustomQueue import Queue_measure
    #measure_queue = Queue_measure()

#while True: # The loop that waits for new instructions...
    
    #data = measure_queue.get() # Get measurement command
    
    ## Do whatever work needs to be done
    ## time.sleep(1)
    
    #measure_queue.publish(data) # Send new results for analysis
    


# Inside SciAnalysis autonomous.py:
########################################
#from CustomQueue import Queue_analyze
#q = Queue_analyze()

#while True: # The loop that waits for new instructions...
    
    #data = q.get() # Get analysis command
    
    ## Do whatever work needs to be done
    ## time.sleep(1)
    
    #q.publish(data) # Send new analysis results to gpCAM




