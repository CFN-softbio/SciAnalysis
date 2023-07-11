import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import glob, os, time, datetime
import imageio, pprint

from SciAnalysis import tools
from SciAnalysis.XSAnalysis.Data import *
from SciAnalysis.XSAnalysis import Protocols
from SciAnalysis.Result import *

#db = databroker.DataBroker.named('cms')
#cat = databroker.catalog['cms']

class experiment():
    def __init__(self, name, folder=None, det='saxs', beamline='cms', ext='tiff', series_measure=False, verbose=0):
        
        if beamline != 'None' or beamline != 'none':
            if verbose>0:  print('At beamline, can use databroker')
            import databroker
        else:
            if verbose>0:  print('Not at beamline, cannot use databroker')

        if folder is None:
            folder = os.getcwd() 
#             self.folder = '/nsls2/data/cms/legacy/xf11bm/data/2022_1/'+'user'+det
        else:
            folder = folder
   
        self.dict = {'expinfo': 
                        {'expname': name,
                        'det': det,
                        'beamline': beamline,
                        'folder': folder,
                        'ext': ext
                        }
                    }
                     
        self.dict['rawinfo'] = {}
        self.dict['rawinfo']['filename'] = []
        self.dict['rawinfo']['time'] = []
        self.dict['rawinfo']['clock'] = []
        self.dict['rawinfo']['scan_id'] = []
        self.dict['rawinfo']['uid'] = []
        self.dict['rawinfo']['filenumber']=0 # total number of input files (could be number of total frames for a series measurement)
        
        self.dict['rawinfo']['series_measure'] = series_measure
        self.dict['rawinfo']['num_frames'] = []
        # if series_measure is True:
        # self.dict['rawinfo']['num_frames'] = 1
        # self.dict['rawinfo']['exposure_period'] = 0.1
    
        self.dict['analysis'] = {}
        
        #self.dict['corr'] = [] # parameters to check correlations
        #self.dict['corrdata'] = {}
 
        #self.dict['mdata_list'] = [] # parameters to pull the metadata
        self.dict['metadata'] = {}        
        self.dict['advanced'] = {} # parameters to run experimental analysis

        
    def show(self, verbose=0):
        print('\n=== Overview of experiment dictionary ===')        
        for key, val in self.dict.items():
            print('exp.dict[\'{}\']'.format(key))
            self._show(key, val, level=0, verbose=verbose)

    def _show(self, key, val, level=0, verbose=0):

        if type(val) == dict:

            if level>0: 
                for ii in np.arange(level): print('  -', end ="")
                print('  key = {}'.format(key))


            keys = list(val.keys())            

            if verbose>0:
                for ii in np.arange(level+1): print('  -', end ="")
                print('  keys = {}'.format(keys))
            else:  
                if len(keys)>0:
                    if list(val.keys())[0].isnumeric() == False: # only print when not index (e.g. '0') 
                        for ii in np.arange(level+1): print('  -', end ="")
                        print('  keys = {}'.format(keys))

            for k, v in val.items():
                if k.isnumeric()==False:
                    self._show(k, v, level=level+1, verbose=verbose)
                    
        else:
            for ii in np.arange(level): print('  -', end ="")
            
            if isinstance(val, np.ndarray)==True and val.shape==():
                print('  key = {}, {}, val = {}'.format(key, type(val), val))  
            elif isinstance(val, np.ndarray)==True and len(val)<10:
                print('  key = {}, {}, val = {}'.format(key, type(val), val))  
            elif isinstance(val, np.ndarray)==True:
                print('  key = {}, {}, val.shape = {}'.format(key, type(val), val.shape))

            elif isinstance(val, list)==False and isinstance(val, np.ndarray)==False:
                print('  key = {}, {}, val = {}'.format(key, type(val), val))
            elif isinstance(val, list)==True and len(val)<10:
                if len(str(val))<10:
                    print('  key = {}, {}, val = {}'.format(key, type(val), val))
                else:
                    print('  key = {}, {}'.format(key, type(val)))
            elif isinstance(val, list)==True:
                print('  key = {}, {}, len(shape) = {}'.format(key, type(val), len(val)))

            else:
                print('  key = {}, {}'.format(key, type(val)))


    
    def defFiles(self, fn = None, scanid=None, uid=None, stitched=False, source_dir = None, verbose=1):
        #define the files in the experiment
        #search raw tiff, return the list of scanid or uid or filenames (RL_t0, RL_t1, ...) 
        #and also look up metadata with the scanid
        #keys = ['sample_x', 'sample_temperature', 'scan_id' ] 
        
        beamline = self.dict['expinfo']['beamline']
        det = self.dict['expinfo']['det']
        folder = self.dict['expinfo']['folder']
        ext = self.dict['expinfo']['ext']

        #self.dict['rawinfo']['series_measure'] = series_measure
        series_measure = self.dict['rawinfo']['series_measure'] 

        t0 = time.time()

        if fn is None:
            fn = self.name
            
        # define the source_dir
        if source_dir is None:
            if stitched == False:
                source_dir = folder + det + '/raw/'
            else:
                source_dir = folder + det + '/stitched/'
        
        if verbose>0:
            print(source_dir)        
        #

        if uid != None and self.beamline != 'None':
            import databroker
            cat = databroker.catalog[self.beamline]
            for uidt in uid:           
                h = cat[uidt]

                self.dict['rawinfo']['filename'].append(h.metadata['start']['filename'])
                self.dict['rawinfo']['time'].append(h.metadata['start']['time']) #linux time
                self.dict['rawinfo']['clock'].append(h.metadata['start']['sample_clock'])
                self.dict['rawinfo']['scan_id'].append(h.metadata['start']['scan_id'])
                self.dict['rawinfo']['uid'].append(h.metadata['start']['uid'])
                if series_measure:
                    self.dict['rawinfo']['num_frames'].append(h.metadata['start']['measure_series_num_frames'])

                
        elif beamline != 'None':
            import databroker
            cat = databroker.catalog[beamline]
            # define infiles
            infiles = []
            if scanid is None: 
                infiles = glob.glob(os.path.join(source_dir, fn + '*.' + ext))
            else:
                for sid in range(scanid[0], scanid[-1]+1):
                    infile = os.path.join(source_dir, fn + '*' + str(sid) + '*.' + ext)
                    if os.isfile(infile):
                        infiles.append(os.path.join(source_dir, fn + '*' + str(sid) + '*.' + ext))

            #sort infiles by the scanid
            infiles = sorted(infiles, key=lambda x: x.split('_'+det)[0].split('_')[-1])        
    #         infiles.sort()
#             print(infiles)
        
            #input exp. info
            for infile in infiles:
                if series_measure:
                    infile_scanid = infile.split('_'+det)[0].split('_')[-2]
                else:
                    infile_scanid = infile.split('_'+det)[0].split('_')[-1]

                h = cat[infile_scanid]
                self.dict['rawinfo']['filename'].append(h.metadata['start']['filename'])
                self.dict['rawinfo']['time'].append(h.metadata['start']['time']) #linux time
                self.dict['rawinfo']['clock'].append(h.metadata['start']['sample_clock'])
                self.dict['rawinfo']['scan_id'].append(h.metadata['start']['scan_id'])
                self.dict['rawinfo']['uid'].append(h.metadata['start']['uid'])
                if series_measure:
                    self.dict['rawinfo']['num_frames'].append(h.metadata['start']['measure_series_num_frames'])

        else:
            # define infiles
            infiles = []
            if scanid is None: 
                infiles = glob.glob(os.path.join(source_dir, fn + '*.' + ext))
            else:
                for sid in range(scanid[0], scanid[-1]+1):
                    infile = os.path.join(source_dir, fn + '*' + str(sid) + '*.' + ext)
                    if os.isfile(infile):
                        infiles.append(os.path.join(source_dir, fn + '*' + str(sid) + '*.' + ext))

            if verbose>0: print('Loading {} files'.format(len(infiles)))

            #sort infiles by the scanid
            if ext == 'tiff':   #beamline =='cms' or beamline=='CMS':
                infiles = sorted(infiles, key=lambda x: x.split('_'+det)[0].split('_')[-1])   
            elif ext == 'tif':  #beamline =='smi' or beamline=='SMI':
                infiles = sorted(infiles, key=lambda x: x.split('id')[1].split('_')[0])   
            
    
            #input exp. info
            num_frames = 1
            for infile in infiles:
                if ext == 'tiff':
                    filename = infile.split('_'+det)[0]
                    filename = filename.split('/')[-1]
                    #scan_id = infile.split('_'+det)[0].split('_')[-1]
                    if series_measure:
                        scan_id = infile.split('_'+det)[0].split('_')[-2]
                        frame = int(infile.split('_'+det)[0].split('_')[-1])
                        if frame+1 > num_frames: num_frames = frame+1
                    else:
                        scan_id = infile.split('_'+det)[0].split('_')[-1]
                elif ext == 'tif':
                    filename = infile.split('_'+det)[0]
                    filename = filename.split('/')[-1]
                    scan_id = infile.split('id')[1].split('_')[0]

                self.dict['rawinfo']['filename'].append(filename)
                #self.dict['rawinfo']['time'].append(h.metadata['start']['time']) #linux time
                self.dict['rawinfo']['scan_id'].append(scan_id)
                self.dict['rawinfo']['num_frames'] = num_frames


        self.dict['rawinfo']['filenumber'] = len(self.dict['rawinfo']['filename'])
        if verbose>0:
            print('Loaded {} files.'.format(self.dict['rawinfo']['filenumber']))

        if self.dict['rawinfo']['series_measure']:
            self.dict['rawinfo']['filenumber'] = self.dict['rawinfo']['num_frames']
        else:
            self.dict['rawinfo']['filenumber'] = len(self.dict['rawinfo']['filename'])

        if verbose>0:
            print('(defFiles time = {:.1f}s)'.format(time.time()-t0))


    def defFiles_ScanID_ONLY(self, fn = None, num_frames=1000, scanid=None, uid=None, stitched=False, source_dir = None, verbose=1):
        #define the files in the experiment
        #search raw tiff, return the list of scanid or uid or filenames (RL_t0, RL_t1, ...) 
        #and also look up metadata with the scanid
        #keys = ['sample_x', 'sample_temperature', 'scan_id' ] 
        
        beamline = self.dict['expinfo']['beamline']
        det = self.dict['expinfo']['det']
        folder = self.dict['expinfo']['folder']
        ext = self.dict['expinfo']['ext']

        #self.dict['rawinfo']['series_measure'] = series_measure
        series_measure = self.dict['rawinfo']['series_measure'] 

        t0 = time.time()

        if fn is None:
            fn = self.name
            
        # define the source_dir
        if source_dir is None:
            if stitched == False:
                source_dir = folder + det + '/raw/'
            else:
                source_dir = folder + det + '/stitched/'
        
        if verbose>0:
            print(source_dir)        
        #

                
        if beamline != 'None':
            import databroker
            cat = databroker.catalog[beamline]
            # define infiles
            # infiles = []
            infiles = glob.glob(os.path.join(source_dir, fn + '*' + str(scanid[0]) + '*.' + ext))

            #sort infiles by the scanid
            infiles = sorted(infiles, key=lambda x: x.split('_'+det)[0].split('_')[-1])        
    #         infiles.sort()
#             print(infiles)
        
            #input exp. info
            for infile in infiles:
                if series_measure:
                    infile_scanid = infile.split('_'+det)[0].split('_')[-2]
                else:
                    infile_scanid = infile.split('_'+det)[0].split('_')[-1]

                h = cat[infile_scanid]

                # query = {'scan_id': {'$gte': scan_id, '$lte': scan_id}}
                # res = cat.search(query)
                # h = cat[list(res)[1]]

                self.dict['rawinfo']['filename'].append(fn)
                self.dict['rawinfo']['time'].append(0) #linux time
                self.dict['rawinfo']['clock'].append(0)
                self.dict['rawinfo']['scan_id'].append(h.metadata['start']['scan_id'])
                self.dict['rawinfo']['uid'].append(h.metadata['start']['uid'])
                if series_measure:
                    self.dict['rawinfo']['num_frames'].append(num_frames)

        else:
            # define infiles
            infiles = []
            if scanid is None: 
                infiles = glob.glob(os.path.join(source_dir, fn + '*.' + ext))
            else:
                for sid in range(scanid[0], scanid[-1]+1):
                    infile = os.path.join(source_dir, fn + '*' + str(sid) + '*.' + ext)
                    if os.isfile(infile):
                        infiles.append(os.path.join(source_dir, fn + '*' + str(sid) + '*.' + ext))

            if verbose>0: print('Loading {} files'.format(len(infiles)))

            #sort infiles by the scanid
            if ext == 'tiff':   #beamline =='cms' or beamline=='CMS':
                infiles = sorted(infiles, key=lambda x: x.split('_'+det)[0].split('_')[-1])   
            elif ext == 'tif':  #beamline =='smi' or beamline=='SMI':
                infiles = sorted(infiles, key=lambda x: x.split('id')[1].split('_')[0])   
            
    
            #input exp. info
            num_frames = 1
            for infile in infiles:
                if ext == 'tiff':
                    filename = infile.split('_'+det)[0]
                    filename = filename.split('/')[-1]
                    #scan_id = infile.split('_'+det)[0].split('_')[-1]
                    if series_measure:
                        scan_id = infile.split('_'+det)[0].split('_')[-2]
                        frame = int(infile.split('_'+det)[0].split('_')[-1])
                        if frame+1 > num_frames: num_frames = frame+1
                    else:
                        scan_id = infile.split('_'+det)[0].split('_')[-1]
                elif ext == 'tif':
                    filename = infile.split('_'+det)[0]
                    filename = filename.split('/')[-1]
                    scan_id = infile.split('id')[1].split('_')[0]

                self.dict['rawinfo']['filename'].append(filename)
                #self.dict['rawinfo']['time'].append(h.metadata['start']['time']) #linux time
                self.dict['rawinfo']['scan_id'].append(scan_id)
                self.dict['rawinfo']['num_frames'] = num_frames


        self.dict['rawinfo']['filenumber'] = len(self.dict['rawinfo']['filename'])
        if verbose>0:
            print('Loaded {} files.'.format(self.dict['rawinfo']['filenumber']))

        if self.dict['rawinfo']['series_measure']:
            self.dict['rawinfo']['filenumber'] = self.dict['rawinfo']['num_frames']
        else:
            self.dict['rawinfo']['filenumber'] = len(self.dict['rawinfo']['filename'])

        if verbose>0:
            print('(defFiles time = {:.1f}s)'.format(time.time()-t0))

        return infiles

    def defFiles_query(self, cycle=None, SAF=None, fn=None, timerange=None, folder=None, scanid=None, verbose=1):

        beamline = self.dict['expinfo']['beamline']
        folder = self.dict['expinfo']['folder']
        det = self.dict['expinfo']['det']

        if beamline != 'None':
            import databroker
            cat = databroker.catalog[beamline]
        else:
            print('Databroker catelog currently not working unless at beamline.')
            print('Use exp.defFiles(fn=sample) instead')
            return()

        if folder is None:
            folder = folder
            
        while folder[-1]=='/':
            folder = folder[:-1]
        folder = folder + '/'
        
        if verbose>0:
            print(folder)

        query = {
             'experiment_alias_directory': {'$in': [folder, folder[:-1]]},
                }
        
        if scanid != None:
            query['scan_id'] = {'$gte': scanid[0], '$lte': scanid[-1]}
        if fn != None:
            #query['sample_name'] = fn  #Sample name does not contain information on x, y, th, etc
            query['filename'] = {'$regex': fn }
        if SAF != None:
            query['experiment_SAF_number'] = SAF
        if timerange != None:
            query['time_range'] = timerange
        if cycle != None:
            query['experiment_cycle'] = cycle

        # ## To save the ROI intensity from the detector
        # if 'data' not in self.dict:
        #     self.dict['data'] = {} # create a dict for data loading
        self.dict['rawinfo']['det']= {}
        for i in range(1,5):
            self.dict['rawinfo']['det'][f'roi{i}']=[]
    
        if self.dict['expinfo']['det'] == 'saxs':
            detector = 'pilatus2M'
        if self.dict['expinfo']['det'] == 'waxs':
            detector = 'pilatus800'
        if self.dict['expinfo']['det'] == 'maxs':
            detector = 'pilatus8002'
            
        t0 = time.time()
        results = cat.search(query)
#         results = results[::-1]
        Nfile = len(results)
        if verbose>0: print('Len(results) = {}'.format(Nfile))        

        for ii, uid in enumerate(results):
            if verbose>0:
                if np.mod(ii, 200)==0: print('[{:.0f}%] '.format(ii/Nfile*100))
            
            if verbose>1: print(uid)
#             infile_scanid = infile.split('_'+det)[0].split('_')[-1]
            h = results[uid]

            self.dict['rawinfo']['filename'].append(h.metadata['start']['filename'])
            self.dict['rawinfo']['time'].append(h.metadata['start']['time']) #linux time
            self.dict['rawinfo']['clock'].append(h.metadata['start']['sample_clock'])
            self.dict['rawinfo']['scan_id'].append(h.metadata['start']['scan_id'])
            self.dict['rawinfo']['uid'].append(h.metadata['start']['uid'])
            
            # series measurements
            if h.metadata['start'].get('measure_type') == 'Series_measure':
                self.dict['rawinfo']['series_measure'] = True
                self.dict['rawinfo']['num_frames'] = h.metadata['start']['measure_series_num_frames']
                self.dict['rawinfo']['exposure_time'] = h.metadata['start']['exposure_time']
                
                if h.metadata['start'].get('exposure_period', None) != None: ## series measurements before 2023 does not have this field
                    self.dict['rawinfo']['exposure_period'] = h.metadata['start']['exposure_period']
            
            #### primary.read() is very slow
#             print(h.primary['data']['pilatus800_stats1_total'])
#             print(h.primary.read()['pilatus800_stats1_total'].values[0])
#             for i in range(1,5):
#                 self.dict['data']['det'][f'roi{i}'].append(h.primary.read()[f'{detector}_stats{i}_total'].values[0])

        if self.dict['rawinfo']['series_measure']:
            self.dict['rawinfo']['filenumber'] = self.dict['rawinfo']['num_frames']
        else:
            self.dict['rawinfo']['filenumber'] = len(self.dict['rawinfo']['uid'])

        print('Loaded {} files, took {}s.'.format(self.dict['rawinfo']['filenumber'], time.time()-t0))

        return results

#     def defMetadata(self, uid, keys=['scan_id']):
#         #load metadata
#         h = cat[uid]
#         for key in keys:
#             self.dict['metadata'][key] = h.start['uid'][key]
    

    def showFileInfo(self, idx=0, verbose=0):
        print('\n=== Raw Information for File {} ==='.format(idx))        
        for key in self.dict['rawinfo'].keys():
            if isinstance(self.dict['rawinfo'][key], list)==True and len(self.dict['rawinfo'][key])>0:
                print('exp.dict[\'rawinfo\'][\'{}\'][{}] = {}'.format(key, idx, self.dict['rawinfo'][key][idx]))


    def loadMetadata(self, keys=None, verbose=1):

        beamline = self.dict['expinfo']['beamline']
        if beamline != 'None':
            import databroker
            cat = databroker.catalog[beamline]
        else:
            print('Databroker catelog currently not working unless at beamline.')
            return()

        t0 = time.time()
        if keys is None:
            keys = self.dict['mdata_list']
        
        for key in keys:
            self.dict['metadata'][key] = []
        
        Nfile = len(self.dict['rawinfo']['uid'])
        for ii, uid in enumerate(self.dict['rawinfo']['uid']):
#         uid = self.dict['rawinfo']['uid'][0]
            if verbose>0:
                if np.mod(ii, 200)==0: print('[{:.0f}%] '.format(ii/Nfile*100))

            h = cat[uid]
            for key in keys:
                self.dict['metadata'][key].append(h.metadata['start'][key])
        print('(loadMetadata time = {:.1f}s)'.format(time.time()-t0))

    #TODO: change fn to self.dict
    def listMetadata(self, scanid=None, uid=None, verbosity=3):
        '''
        list of the keys in Metadata to input
        '''
        beamline = self.dict['expinfo']['beamline'] 
        if beamline != 'None':
            import databroker
            cat = databroker.catalog[beamline]
        else:
            print('Databroker catelog currently not working unless at beamline.')
            return();

    
        if scanid is None:
            scanid = self.dict['rawinfo']['uid'][0]

        h = cat[scanid]
        if verbosity>= 3:
            for key in h.metadata['start']:
                print(key)

        #return mdata_list
        if verbosity<=-1:
            mdata_list = []
            num = input("How many variables to be selected in metadata? (E.g. 2): ")
            for ii in range(int(num)):
                val = input("Input variable NO. {} (E.g. sample_x): ".format(ii))
                if val in h.metadata['start']:
                    
                    mdata_list.append(val)
                else:
                    print('Error input')
                    return
            self.dict['mdata_list'] = mdata_list
#             self.dict['metadata'] = mdata_list
            return mdata_list
                
                #         return scanids

    def showMetadata(self, scanid=None, uid=None, md_interest = None, verbosity=1):
        # md_interest = ['scan_id']

        beamline = self.dict['expinfo']['beamline'] 
        if beamline != 'None':
            import databroker
            cat = databroker.catalog[beamline]
        else:
            print('Databroker catelog currently not working unless at beamline.')
            return();
    
        if scanid is None:
            scanid = self.dict['rawinfo']['uid'][0]

        h = cat[scanid]

        if md_interest != None:
            print('### Scan {}:'.format(scanid))
            for md in md_interest:
                print('{}: {}\n'.format(md, h.metadata['start'][md]))
        else:
            print('### Scan {}:'.format(scanid))
            pprint.pprint(h.metadata['start'])



    #load the reducted data (SciAanalysis) files via the key 
    #key = 'circurlar_average'  or 'linecut_qr-q=0.1'
    #filename full path: 'circular_average/RL_t0.dat',  'qr_image/RL_t0.npz'
    #input: key = 'circular_average', filename = 'RL_t0.dat'

    
    
    def loadSciAnalysisData(self, keys=None, analysis_folder=None, verbose=0):
        t0 = time.time()

        folder = self.dict['expinfo']['folder']
        det = self.dict['expinfo']['det']

        if 'analysis' not in self.dict:
            self.dict['analysis'] = {} # create a dict for data loading
        

        if analysis_folder is None:
            analysis_folder = folder + '/' + det + '/analysis/'
        if verbose > 0: print('analysis_folder = {}'.format(analysis_folder))

        if keys is None:
            folders = glob.glob(analysis_folder + '/*/')
            keys = []
            for folder in folders:
                key = folder.split('/')[-2]
                keys.append(key)
        
        ### Filenames to load
        Nfile = self.dict['rawinfo']['filenumber']

        # for regular scan/snap measurements
        if not self.dict['rawinfo']['series_measure']:
            filenames = self.dict['rawinfo']['filename']
        
        # series measurements
        else: 
            infile = self.dict['rawinfo']['filename'][0]
            
            # to remove extension
            if infile[:-5] == '.tiff':
                infile = infile[:-5]
            
            # to remove scanid in the filename for the data before 2023 (incorrect scanid) and add the exposure_period
            jan2023 = time.mktime(datetime.datetime.strptime('01/01/2023',"%m/%d/%Y").timetuple())
            if self.dict['rawinfo']['time'][0] < jan2023:
                infile = '_'.join(infile.split('_')[:-1]) 
                exposure_period = float(infile.split('_')[-1].split('s')[0])
                self.dict['rawinfo']['exposure_period'] = exposure_period
                scan_id = self.dict['rawinfo']['scan_id'][0]
                infile = '_'.join([infile,str(scan_id+1)])

            filenames = ['_'.join([infile,str(kk).zfill(6)]) for kk in range(Nfile)]
        
        ### load files
        for nn, infile in enumerate(filenames):
            if verbose>0: 
                if np.mod(nn, 200)==0: print('[{:.0f}%] '.format(nn/Nfile*100))
            if verbose>1: print('Searching analysis results for {}'.format(infile))

            for key in keys:
                
                if key not in self.dict['analysis']:
                    self.dict['analysis'][key] = {}
                
                if verbose>1: print(os.path.join(analysis_folder+key, infile+'*.dat'))

                if 'average' in key:
                    files = glob.glob(os.path.join(analysis_folder+key, infile+'*.dat'))
                    if len(files) == 0:
                        return print('There is no data in the folder {}. '.format(key))
                    else:
                        file=files[0]                    
                    headers = pd.read_csv(file, delim_whitespace=True, nrows=0).columns[1:]
                    dat = pd.read_csv(file, delim_whitespace=True, header=None, skiprows=1, names=headers)
                    
                    self.dict['analysis'][key][str(nn)] = {}
                    self.dict['analysis'][key][str(nn)][headers[0]] = dat[headers[0]].values
                    self.dict['analysis'][key][str(nn)][headers[2]] = dat[headers[2]].values
                    

                if 'linecut' in key:
                    files = glob.glob(os.path.join(analysis_folder + key, infile + '*.dat'))
                    if len(files) == 0:
                        return print('There is no data in the folder {}. '.format(key))
                    else:
                        file=files[0]                       
                    headers = pd.read_csv(file, delim_whitespace=True, nrows=0).columns[1:]
                    dat = pd.read_csv(file, delim_whitespace=True, header=None, skiprows=1, names=headers)
                    
                    self.dict['analysis'][key][str(nn)] = {}
                    self.dict['analysis'][key][str(nn)][headers[0]] = dat[headers[0]].values
                    self.dict['analysis'][key][str(nn)][headers[1]] = dat[headers[1]].values
                    
                if 'image' in key:
                    files = glob.glob(os.path.join(analysis_folder + key, infile + '*.png'))
                    if len(files) == 0:
                        return print('There is no data in the folder {}. '.format(key))
                    else:
                        file=files[0]                    
                    self.dict['analysis'][key][str(nn)] = imageio.imread(file)
                
                if 'roi' in key: # load roi results from xml file
                    files = glob.glob(os.path.join(analysis_folder + 'results', infile + '*.xml'))
                    if len(files) == 0:
                        return print('There is no data in the folder {}. '.format(key))
                    else:
                        file=files[0]                    
                    names, values = Results().extract_results_from_xml(file, protocol='roi', verbosity=3)
                    for kk,name in enumerate(names):
                        if name not in self.dict['analysis'][key]:
                            self.dict['analysis'][key][name]= []
                        self.dict['analysis'][key][name].append(values[kk])  

        print('loadSciAnalysisData time = {:.1f}s'.format(time.time()-t0))             
                    

    #TODO
    def doPlot(self, protocol, key):
        pass
        

        

    #TODO
    #plotting all heat maps
    def plotHeatmap(self, key, y_axes=None, flag_log = 1, plot_xrange=None):
        '''
        plot the heatmap for the 2D arrays
        '''

        # print('0')

        xy_axes = list(self.dict['analysis'][key][str(0)].keys())  
        x_axis = self.dict['corrdata']['2Darray'][key][xy_axes[0]]
        I_array = self.dict['corrdata']['2Darray'][key]['I_array']
        if flag_log==1:
            I_array = np.log10(I_array)

        y_idx = np.arange(self.dict['rawinfo']['filenumber'])
        if plot_xrange is None:
            plot_xrange = np.arange(len(x_axis))

        # print('1')
        if y_axes is not None:
            y_axis = self.dict['corrdata']['2Darray'][key][y_axes[0]]            
        
        fig = plt.figure() 
        plt.clf()
        plt.rcParams['figure.figsize'] = [8, 4] #[height, width] of the figure
        plt.rcParams['xtick.labelsize'] = 12
        plt.rcParams['ytick.labelsize'] = 12
        plt.rcParams['axes.labelsize'] = 12
        plot_buffers = plot_buffers = [0.2, 0.1, 0.2, 0.1]
        left_buf, right_buf, bottom_buf, top_buf = plot_buffers
        fig_width = 1.0-right_buf-left_buf
        fig_height = 1.0-top_buf-bottom_buf
        fig.add_axes( [left_buf, bottom_buf, fig_width, fig_height] )

        # print('2')
        ## 2D plot
        plt.pcolormesh(x_axis[plot_xrange], y_idx, I_array[:,plot_xrange])

        plt.xlabel('{}.(A^-1 or azi angle)'.format(xy_axes[0]))
        plt.ylabel('index')
        if y_axes is not None:
            plt.yticks(y_idx, y_axis)
            plt.ylabel(y_axes[0])
        plt.colorbar()
        # print('3')

        ## plot the 1D
        y_params = list(self.dict['corrdata']['2Darray'][key].keys())

        y_params.remove(xy_axes[0])
        y_params.remove('file_index')
        y_params.remove('I_array')
        
        cmap = plt.get_cmap('jet',len(y_params))
        
        fig, axs = plt.subplots(len(y_params), 1, sharex=True, figsize = [8,4])

        # Remove horizontal space between axes
        fig.subplots_adjust(hspace=0)

        if len(y_params) == 1:
            axs.plot(y_idx, self.dict['corrdata']['2Darray'][key][y_params[0]], color = cmap(0), marker='o')
            axs.set_ylabel(y_params[0])
            axs.set_xlabel('index')
            
        elif len(y_params) > 1:
            for kk, ax in enumerate(axs.flat):
                ax.plot(y_idx, self.dict['corrdata']['2Darray'][key][y_params[kk]], color = cmap(kk), marker='o')
                ax.set_ylabel(y_params[kk])
            ax.set_xlabel('index')
#         plt.tight_layout()
            
#         fig,ax = plt.subplots()
#         for kk, y_param in enumerate(y_params):
#             if kk == 0:
#                 ax.plot(y_axes, self.dict['analysis']['2Darray'][key][y_param], color = cmap(kk))
#                 ax.set_ylabel(y_param, color = cmap(kk))
#             if kk >=1:
#                 ax_new = ax.twinx()
#                 ax_new.plot(y_axes, self.dict['analysis']['2Darray'][key][y_param], color = cmap(kk))
#                 ax_new.set_ylabel(y_param, color = cmap(kk))

#             ax.set_xlabel('idx')
        # plt.show()
        

    def plotWaterfall(self, key, y_spacing = None, gridlines=True, flag_log = [0, 0], plot_xrange=None, lw = 1):
        '''
        plot waterfall for the 2D arrays
        '''
        xy_axes = list(self.dict['analysis'][key][str(0)].keys())  
        x_axis = self.dict['corrdata']['2Darray'][key][xy_axes[0]]
        I_array = self.dict['corrdata']['2Darray'][key]['I_array']
        if flag_log[1]==1:
            I_array = np.log10(I_array)
        if flag_log[0]==1:
            x_axis = np.log10(x_axis)

        if plot_xrange is None:
            plot_xrange = np.arange(len(x_axis))
    
        y_idx = np.arange(self.dict['rawinfo']['filenumber'])
        
        
        fig = plt.figure() 
        plt.clf()
        plt.rcParams['figure.figsize'] = [8, 4] #[height, width] of the figure
        plt.rcParams['xtick.labelsize'] = 12
        plt.rcParams['ytick.labelsize'] = 12
        plt.rcParams['axes.labelsize'] = 12
        plot_buffers = plot_buffers = [0.2, 0.1, 0.2, 0.1]
        left_buf, right_buf, bottom_buf, top_buf = plot_buffers
        fig_width = 1.0-right_buf-left_buf
        fig_height = 1.0-top_buf-bottom_buf
        fig.add_axes( [left_buf, bottom_buf, fig_width, fig_height] )
        
        ## 2D plot
        if y_spacing is None:
            y_spacing = np.mean(I_array)*2
        N = len(I_array)


        cmap = mpl.colormaps['viridis']
        colors = cmap(np.linspace(0.0, 1.0, N))

        for nn, curve in enumerate(I_array):
            plt.plot(x_axis[plot_xrange], curve[plot_xrange] + y_spacing*nn, color=colors[nn], lw=lw)

        plt.xlabel('{}.(A^-1 or azi angle)'.format(xy_axes[0]))
        plt.ylabel('Intensity')
        if gridlines:
            plt.grid()
        plt.title('Experiment {}, {}'.format(self.name, key))

    
    
    def addCorr(self, fitting=False, verbosity=3):
        
        #print out the list of data            
        if verbosity>=3:
            for key in self.dict['analysis'].keys():
                print(key)
        data = input('select the dataset:')

        #print out the list of corr parameters
        if verbosity>=3:
            for key in self.dict['metadata'].keys():
                print(key)       
        mdata = input('select the metada:')
        
        #TODO:locate the fitting keys
        #print out the list of fitting results
        if fitting==True:
            if verbosity>=3:
                xml = loadxml
                
                for key in xml.keys:
                    print(key)       
            mdata = input('select the metada:')
            
        
        if fitting == True:
            self.dict['corr'].append([data, medata, fitting])
        else:
            self.dict['corr'].append([data, medata])
        
    def listCorr(self, verbosity=3):        
        for corr in self.dict['corr']:
            print(corr)
        
    def doCorr(self, corrs=None):
        '''
            input : 
                corr = [['2Darray'],
                        ['circular_average','temperature'],
                        # For peak q0, plot peak intensity (or width etc) vs. sample position x
                        ['circular_average_fitting',  'sample_x', 'peak_q',] 
                        # 2D mapping in smx and smy : For peak q0, plot peak intensity (or width etc) vs (x,y)
                            ]

                2Darray 
                average * tempreature: 2D array 
                circular_average_fitting*sample_x*peak_q]

            results:
                self.dict['corrdata']['circular_average__temperature'] #save 2D array, (optional export PNG)
                self.dict['corrdata']['peak_q__sample_x']
                self.dict['corrdata']['2Darray]

        '''    
        if corrs is None:
            corrs = self.dict['corr']
            
        for corr in corrs:
            
            #2D array
            if corr == ['2Darray']:
                self.dict['corrdata'][corr[0]] = {}
                for key in self.dict['analysis'].keys():
                    if ('average' in key) or ('linecut' in key):
                        
                        # Data in a 2D dict
                        # 2D array: 'Iq_array' (nrow x ncol)
                        # 1D array: 'q' (len = ncol)
                        # 1D array: 'scan_id', 'time', and other parameters in the 'mdata_list' (len = nrow)
                        #            for series measurement: 'time_series'

                        print('doCorr {}'.format(key))
                        self.dict['corrdata'][corr[0]][key] = {}
                        
                        xy_axes = list(self.dict['analysis'][key][str(0)].keys())
                        self.dict['corrdata'][corr[0]][key][xy_axes[0]] = self.dict['analysis'][key][str(0)][xy_axes[0]]

                        qIq = np.array([v for _,one_scan in self.dict['analysis'][key].items() 
                                          for k,v in one_scan.items()])
                        self.dict['corrdata'][corr[0]][key]['I_array'] = qIq[1::2,:]
                        
                        ## Include corr parameter 
                        # add basic param (always), such as 'file_index'
                        self.dict['corrdata'][corr[0]][key]['file_index'] = np.array(list(self.dict['analysis'][key].keys()))

                        if self.dict['rawinfo']['series_measure']: # Add 'time_series' for series_measurement
                            self.dict['corrdata'][corr[0]][key]['time_series'] = np.arange(self.dict['rawinfo']['num_frames']) * self.dict['rawinfo']['exposure_period']
                        else:
                            corr_params = ['scan_id', 'time'] # Add 'scan_id' and 'time' for scan/snap measurements
                            for basic_param in corr_params:
                                self.dict['corrdata'][corr[0]][key][basic_param] = np.array(self.dict['rawinfo'][basic_param])       

                        ## Include corresponding metadata
                        for param in self.dict['mdata_list']:
                            self.dict['corrdata'][corr[0]][key][param] = np.array(self.dict['metadata'][param])
          
                            
            # if 'reduction' in corr:

            #     print('reduction is working')
            #     for item in set(self.dict['metadata'][corr[-1]]):
            #         reduction_protocol = data + corr[-1] + str(item)
            #         self.dict['analysis'][reduction_protocol] = []
            #         for ii, data in enumerate(self.dict['metadata'][corr]):
            #             if data == item:

            #                 for key in self.dict['analysis'].keys():
            #                     if ('average' in key) or ('linecut' in key):
                                    
            #                         # Data in a 2D dict
            #                         # 2D array: 'Iq_array' (nrow x ncol)
            #                         # 1D array: 'q' (len = ncol)
            #                         # 1D array: 'scan_id', 'time', and other parameters in the 'mdata_list' (len = nrow)

            #                         self.dict['analysis'][protocol[0]] = {}
            #                         print(key)
            #                         self.dict['analysis'][protocol[0]][key] = {}
                                    
                                    
            #                         xy_axes = list(self.dict['analysis'][key][0].keys())
                                    
            #                         self.dict['analysis'][protocol[0]][key][xy_axes[0]] = self.dict['analysis'][key][0][xy_axes[0]] 
                                    
            #                         ## Include basic parameters (always)
            #                         basic_params = ['scan_id', 'time']
            #                         for basic_param in basic_params:
            #                             self.dict['analysis'][protocol[0]][key][basic_param] = np.array(self.dict['rawinfo'][basic_param])       
                                
            #                         ## Include corresponding metadata
            #                         for param in self.dict['mdata_list']:
            #                             self.dict['analysis'][protocol[0]][key][param] = np.array(self.dict['metadata'][param])
                                    
            #                         # Iq 2D array
            #                         self.dict['analysis'][protocol[0]][key]['I_array'] = np.array([
            #                             self.dict['analysis'][key][i][xy_axes[1]] for i in range(self.dict['rawinfo']['filenumber'])])        
#                         #load x as the first row
#                         self.dict['analysis'][protocol][key] = self.dict['analysis'][key]
                        
#                         I_list = []#2D array for data stacking
#                         for ii in self.dict['analysis'][key]:
#                             #TODO:
#                             I_list.append(self.dict['analysis'][key])
                            
                            
#                             self.dict['analysis'][protocol]['sequence'][ii] = ii
                            
#                         self.dict['analysis'][protocol]['analysis'] = np.asarray(I_list)                                
#                         self.dict['analysis'][protocol]['scan_id'] = self.dict['rawinfo']['scan_id'] 
            #1D waterflow plots will be the plotting issue. 
            
            
#             #2D mapping in smx and smy
#             if protocol == ['mapping']:
#                 #TODO : check smx and smy or smy2
#                 #use one ROI or peak intensity as the 2D mapping images.
#                 pass
                
# #                         self.dict['analysis']['peak_q__sample_x']
                    
#             #1D plot, x=metadata, y=fitting results
#             elif len(protocol) >= 3:
                
#                 #define the name of correlated protocols
#                 pname=[]
#                 for pp in protocol:
#                     pname = pp+'*'+pname
#                     if 'average' or 'linecut' in pp: #identify the reducted data as pdata
#                         pdata = pp
                
#                 for pp in protocol:
                    
#                     if 'average' or 'linecut' in pp: #add the raw data and fitting curve into 2D array
#                         pass
#                     elif pp in self.dict['metadata']:
#                         #x-axis
#                         self.dict['analysis'][pname][pp] = self.dict['analysis'][pp]
                        
#                     else: #pp = peak_position, peak_width, the same as .xml
#                         self.dict['analysis'][pname]['y_axis'] = self.dict['analysis'][key]
#                         self.dict['analysis'][pname]['scan_id'] = self.dict['rawinfo']['scan_id'] 

#                         for index, data in enumerate(self.dict['analysis'][pdata]):
#                             self.dict['analysis'][pname]['sequence'][index] = index
#                             #read xml for the results
#                             fittingresult = self.loadxml(pp)
#                             self.dict['analysis'][pname][pp][index] = fittingresult
    
    
    
    # TODO
    def doExpAnalysis(self, exp_protocol=None):
        '''
        exp analysis
        '''          
        pass
    
    def listExpProtocol(self, verbosity=3):        
        for protocol in self.dict['exp_protocol']:
            print(protocol)
    
    # TODO
    def addExpProtocol(self, fitting=False, verbosity=3):
        
        #print out the list of data            
        if verbosity>=3:
            for key in self.dict['analysis'].keys():
                print(key)
        data = input('select the dataset:')

        #print out the list of protocols
        if verbosity>=3:
            for key in self.dict['metadata'].keys():
                print(key)       
        mdata = input('select the metada:')
        
        #TODO:locate the fitting keys
        #print out the list of fitting results
        if fitting==True:
            if verbosity>=3:
                xml = loadxml
                
                for key in xml.keys:
                    print(key)       
            mdata = input('select the metada:')
            
        
        if fitting == True:
            self.dict['exp_protocol'].append([data, medata, fitting])
        else:
            self.dict['exp_protocol'].append([data, medata])

    def doDict(self):
        #print out of all steps to make self.dict

        print('This is the instructon to make the analysis dictionary.')
        print('Step1: defFiles   : define the list of files in the experiment')
        print('Step2: loadMetadata  ')
        print('Step3: loadSciAnalysisData  : load reduced data from SciAnalysis')
        print('Step4: addProtocol   :  add individual analysis protocol')
        print('Step5: doAnalysis    :  ')
        print('Step6: publishH5     :  ')

        print('DONE. Enjoy the plotting')


      
        
        
        
    def publishH5(self, mdata_dict):
#         dicttoH5
        pass



#exp = experiment('BestExpEver_SAXS')
#exp = experiment('BestExpEver_MAXS')



