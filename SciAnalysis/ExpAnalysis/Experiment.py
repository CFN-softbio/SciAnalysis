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
    def __init__(self, name, folder=None, det='saxs', beamline='cms'):
        
        self.name = name
        self.det = det
        self.beamline = beamline
        if beamline is not None:
            print('At beamline, can use databroker')
            import databroker
        else:
            print('Not at beamline, cannot use databroker')
        
        if folder is None:
            self.folder = os.getcwd() 
#             self.folder = '/nsls2/data/cms/legacy/xf11bm/data/2022_1/'+'user'+det
        else:
            self.folder = folder
   
        self.dict = {'expinfo': 
                        {'expname': self.name,
                        'detector': self.det,
                        'beamline': self.beamline,
                        'folder': folder
                        }
                    }
                     
        #self.dict['expinfo'] = {}
        self.dict['expinfo']['filename'] = []
        self.dict['expinfo']['time'] = []
        self.dict['expinfo']['clock'] = []
        self.dict['expinfo']['scan_id'] = []
        self.dict['expinfo']['uid'] = []
        self.dict['expinfo']['filenumber']=0 # total number of input files (could be number of total frames for a series measurement)
        
        self.dict['expinfo']['series_measure'] = False
        # if series_measure is True:
        # self.dict['expinfo']['num_frames'] = 1
        # self.dict['expinfo']['exposure_period'] = 0.1
    
        self.dict['data'] = {}
        
        self.dict['corr'] = [] # parameters to check correlations
        self.dict['corrdata'] = {}
 
        self.dict['mdata_list'] = [] # parameters to pull the metadata
        self.dict['metadata'] = {}
        
        self.dict['exp_protocol'] = [] # parameters to run experimental analysis
        self.dict['analysis'] = {}
        
    def show(self, verbose=0):
        print('=== Overview of experiment dictionary ===')        
        for key, val in self.dict.items():
            print('exp.dict[\'{}\']'.format(key))
            self._show(key, val, level=0, verbose=verbose)

    def _show(self, key, val, level=0, verbose=0):
           
        if type(val) == dict:
            keys = list(val.keys())            

            if verbose>0:
                for ii in np.arange(level): print('    -', end ="")
                print('    keys = {}'.format(keys))
            else:  
                if len(keys)>0:
                    if list(val.keys())[0].isnumeric() == False: # only print when not index (e.g. '0') 
                        for ii in np.arange(level): print('    -', end ="")
                        print('    keys = {}'.format(keys))

            for k, v in val.items():
                if k.isnumeric()==False:
                    self._show(k, v, level=level+1, verbose=verbose)
                    
        else:
            for ii in np.arange(level): print('    -', end ="")
            
            if isinstance(val, np.ndarray)==True and val.shape==():
                print('    key = {}, val = {}'.format(key, val))  
            elif isinstance(val, np.ndarray)==True and len(val)<10:
                print('    key = {}, val = {}'.format(key, val))  
            elif isinstance(val, list)==False and isinstance(val, np.ndarray)==False:
                print('    key = {}, val = {}'.format(key, val))
            else:
                print('    key = {}, val.shape = {}'.format(key, val.shape))



    
    def defFiles(self, fn = None, scanid=None, uid=None, stitched=False, burstmode=False, verbose=1):
        #define the files in the experiment
        #search raw tiff, return the list of scanid or uid or filenames (RL_t0, RL_t1, ...) 
        #and also look up metadata with the scanid
        #keys = ['sample_x', 'sample_temperature', 'scan_id' ] 
        
        t0 = time.time()

        if fn is None:
            fn = self.name
            
        # define the source_dir
        if stitched == False:
            source_dir = self.folder + self.det + '/raw/'
        else:
            source_dir = self.folder + self.det + '/stitched/'
        if verbose>0:
            print(source_dir)        
        #

        if uid != None and self.beamline is not None:
            import databroker
            cat = databroker.catalog[self.beamline]
            for uidt in uid:           
                h = cat[uidt]
                self.dict['expinfo']['filename'].append(h.metadata['start']['filename'])
                self.dict['expinfo']['time'].append(h.metadata['start']['time']) #linux time
                self.dict['expinfo']['clock'].append(h.metadata['start']['sample_clock'])
                self.dict['expinfo']['scan_id'].append(h.metadata['start']['scan_id'])
                self.dict['expinfo']['uid'].append(h.metadata['start']['uid'])

                
        elif self.beamline is not None:
            import databroker
            cat = databroker.catalog[self.beamline]
            # define infiles
            infiles = []
            if scanid is None: 
                infiles = glob.glob(os.path.join(source_dir, fn + '*.tiff'))
            else:
                for sid in range(scanid[0], scanid[-1]+1):
                    infile = os.path.join(source_dir, fn + '*' + str(sid) + '*.tiff')
                    if os.isfile(infile):
                        infiles.append(os.path.join(source_dir, fn + '*' + str(sid) + '*.tiff'))

            #sort infiles by the scanid
            infiles = sorted(infiles, key=lambda x: x.split('_'+self.det)[0].split('_')[-1])        
    #         infiles.sort()
#             print(infiles)
        
            #input exp. info
            for infile in infiles:

                infile_scanid = infile.split('_'+self.det)[0].split('_')[-1]

                h = cat[infile_scanid]
                self.dict['expinfo']['filename'].append(h.metadata['start']['filename'])
                self.dict['expinfo']['time'].append(h.metadata['start']['time']) #linux time
                self.dict['expinfo']['clock'].append(h.metadata['start']['sample_clock'])
                self.dict['expinfo']['scan_id'].append(h.metadata['start']['scan_id'])
                self.dict['expinfo']['uid'].append(h.metadata['start']['uid'])

        else:
            # define infiles
            infiles = []
            if scanid is None: 
                infiles = glob.glob(os.path.join(source_dir, fn + '*.tiff'))
            else:
                for sid in range(scanid[0], scanid[-1]+1):
                    infile = os.path.join(source_dir, fn + '*' + str(sid) + '*.tiff')
                    if os.isfile(infile):
                        infiles.append(os.path.join(source_dir, fn + '*' + str(sid) + '*.tiff'))

            #sort infiles by the scanid
            infiles = sorted(infiles, key=lambda x: x.split('_'+self.det)[0].split('_')[-1])        
            #input exp. info
            for infile in infiles:
                filename = infile.split('_'+self.det)[0]
                filename = filename.split('/')[-1]
                scan_id = infile.split('_'+self.det)[0].split('_')[-1]


                self.dict['expinfo']['filename'].append(filename)
                #self.dict['expinfo']['time'].append(h.metadata['start']['time']) #linux time
                self.dict['expinfo']['scan_id'].append(scan_id)


        self.dict['expinfo']['filenumber'] = len(self.dict['expinfo']['filename'])
        if verbose>0:
            print('Loaded {} files.'.format(self.dict['expinfo']['filenumber']))

        if self.dict['expinfo']['series_measure']:
            self.dict['expinfo']['filenumber'] = self.dict['expinfo']['num_frames']
        else:
            self.dict['expinfo']['filenumber'] = len(self.dict['expinfo']['filename'])

        if verbose>0:
            print('(defFiles time = {:.1f}s)'.format(time.time()-t0))
                
    def defFiles_query(self, cycle=None, SAF=None, fn=None, timerange=None, folder=None, scanid=None, verbose=1):

        if self.beamline is not None:
            import databroker
            cat = databroker.catalog[self.beamline]
        else:
            print('Databroker catelog currently not working unless at beamline.')
            print('Use exp.defFiles(fn=sample) instead')
            return()

        if folder is None:
            folder = self.folder
            
        while folder[-1]=='/':
            folder = folder[:-1]
        folder = folder + '/'
        
        if verbose>0:
            print(folder)

        query = {
             'experiment_alias_directory': {'$in': [folder, folder[:-1]]},
                }
        
        if scanid is not None:
            query['scan_id'] = {'$gte': scanid[0], '$lte': scanid[-1]}
        if fn is not None:
            #query['sample_name'] = fn  #Sample name does not contain information on x, y, th, etc
            query['filename'] = {'$regex': fn }
        if SAF is not None:
            query['experiment_SAF_number'] = SAF
        if timerange is not None:
            query['time_range'] = timerange
        if cycle is not None:
            query['experiment_cycle'] = cycle

        # ## To save the ROI intensity from the detector
        # if 'data' not in self.dict:
        #     self.dict['data'] = {} # create a dict for data loading
        self.dict['data']['det']= {}
        for i in range(1,5):
            self.dict['data']['det'][f'roi{i}']=[]
    
        if self.dict['type'] == 'saxs':
            detector = 'pilatus2M'
        if self.dict['type'] == 'waxs':
            detector = 'pilatus800'
        if self.dict['type'] == 'maxs':
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
#             infile_scanid = infile.split('_'+self.det)[0].split('_')[-1]
            h = results[uid]

            self.dict['expinfo']['filename'].append(h.metadata['start']['filename'])
            self.dict['expinfo']['time'].append(h.metadata['start']['time']) #linux time
            self.dict['expinfo']['clock'].append(h.metadata['start']['sample_clock'])
            self.dict['expinfo']['scan_id'].append(h.metadata['start']['scan_id'])
            self.dict['expinfo']['uid'].append(h.metadata['start']['uid'])
            
            # series measurements
            if h.metadata['start'].get('measure_type') == 'Series_measure':
                self.dict['expinfo']['series_measure'] = True
                self.dict['expinfo']['num_frames'] = h.metadata['start']['measure_series_num_frames']
                self.dict['expinfo']['exposure_time'] = h.metadata['start']['exposure_time']
                
                if h.metadata['start'].get('exposure_period', None) is not None: ## series measurements before 2023 does not have this field
                    self.dict['expinfo']['exposure_period'] = h.metadata['start']['exposure_period']
            
            #### primary.read() is very slow
#             print(h.primary['data']['pilatus800_stats1_total'])
#             print(h.primary.read()['pilatus800_stats1_total'].values[0])
#             for i in range(1,5):
#                 self.dict['data']['det'][f'roi{i}'].append(h.primary.read()[f'{detector}_stats{i}_total'].values[0])

        if self.dict['expinfo']['series_measure']:
            self.dict['expinfo']['filenumber'] = self.dict['expinfo']['num_frames']
        else:
            self.dict['expinfo']['filenumber'] = len(self.dict['expinfo']['uid'])

        print('Loaded {} files, took {}s.'.format(self.dict['expinfo']['filenumber'], time.time()-t0))

        return results

#     def defMetadata(self, uid, keys=['scan_id']):
#         #load metadata
#         h = cat[uid]
#         for key in keys:
#             self.dict['metadata'][key] = h.start['uid'][key]
    




    def loadMetadata(self, keys=None, verbose=1):

        if self.beamline is not None:
            import databroker
            cat = databroker.catalog[self.beamline]
        else:
            print('Databroker catelog currently not working unless at beamline.')
            return()

        t0 = time.time()
        if keys is None:
            keys = self.dict['mdata_list']
        
        for key in keys:
            self.dict['metadata'][key] = []
        
        Nfile = len(self.dict['expinfo']['uid'])
        for ii, uid in enumerate(self.dict['expinfo']['uid']):
#         uid = self.dict['expinfo']['uid'][0]
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

        if self.beamline is not None:
            import databroker
            cat = databroker.catalog[self.beamline]
        else:
            print('Databroker catelog currently not working unless at beamline.')
            return();

    
        if scanid is None:
            scanid = self.dict['expinfo']['uid'][0]

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

        if self.beamline is not None:
            import databroker
            cat = databroker.catalog[self.beamline]
        else:
            print('Databroker catelog currently not working unless at beamline.')
            return();
    
        if scanid is None:
            scanid = self.dict['expinfo']['uid'][0]

        h = cat[scanid]

        if md_interest is not None:
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
        if 'data' not in self.dict:
            self.dict['data'] = {} # create a dict for data loading
        

        if analysis_folder is None:
            analysis_folder = self.folder + '/' + self.det + '/analysis/'
        if verbose > 0: print('analysis_folder = {}'.format(analysis_folder))

        if keys is None:
            folders = glob.glob(analysis_folder + '/*/')
            keys = []
            for folder in folders:
                key = folder.split('/')[-2]
                keys.append(key)
        
        ### Filenames to load
        Nfile = self.dict['expinfo']['filenumber']

        # for regular scan/snap measurements
        if not self.dict['expinfo']['series_measure']:
            filenames = self.dict['expinfo']['filename']
        
        # series measurements
        else: 
            infile = self.dict['expinfo']['filename'][0]
            
            # to remove extension
            if infile[:-5] == '.tiff':
                infile = infile[:-5]
            
            # to remove scanid in the filename for the data before 2023 (incorrect scanid) and add the exposure_period
            jan2023 = time.mktime(datetime.datetime.strptime('01/01/2023',"%m/%d/%Y").timetuple())
            if self.dict['expinfo']['time'][0] < jan2023:
                infile = '_'.join(infile.split('_')[:-1]) 
                exposure_period = float(infile.split('_')[-1].split('s')[0])
                self.dict['expinfo']['exposure_period'] = exposure_period
                scan_id = self.dict['expinfo']['scan_id'][0]
                infile = '_'.join([infile,str(scan_id+1)])

            filenames = ['_'.join([infile,str(kk).zfill(6)]) for kk in range(Nfile)]
        
        ### load files
        for nn, infile in enumerate(filenames):
            if verbose>0: 
                if np.mod(nn, 200)==0: print('[{:.0f}%] '.format(nn/Nfile*100))
            if verbose>1: print('Searching analysis results for {}'.format(infile))

            for key in keys:
                
                if key not in self.dict['data']:
                    self.dict['data'][key] = {}
                
                if verbose>1: print(os.path.join(analysis_folder+key, infile+'*.dat'))

                if 'average' in key:
                    files = glob.glob(os.path.join(analysis_folder+key, infile+'*.dat'))
                    if len(files) == 0:
                        return print('There is no data in the folder {}. '.format(key))
                    else:
                        file=files[0]                    
                    headers = pd.read_csv(file, delim_whitespace=True, nrows=0).columns[1:]
                    dat = pd.read_csv(file, delim_whitespace=True, header=None, skiprows=1, names=headers)
                    
                    self.dict['data'][key][str(nn)] = {}
                    self.dict['data'][key][str(nn)][headers[0]] = dat[headers[0]].values
                    self.dict['data'][key][str(nn)][headers[2]] = dat[headers[2]].values
                    

                if 'linecut' in key:
                    files = glob.glob(os.path.join(analysis_folder + key, infile + '*.dat'))
                    if len(files) == 0:
                        return print('There is no data in the folder {}. '.format(key))
                    else:
                        file=files[0]                       
                    headers = pd.read_csv(file, delim_whitespace=True, nrows=0).columns[1:]
                    dat = pd.read_csv(file, delim_whitespace=True, header=None, skiprows=1, names=headers)
                    
                    self.dict['data'][key][str(nn)] = {}
                    self.dict['data'][key][str(nn)][headers[0]] = dat[headers[0]].values
                    self.dict['data'][key][str(nn)][headers[1]] = dat[headers[1]].values
                    
                if 'image' in key:
                    files = glob.glob(os.path.join(analysis_folder + key, infile + '*.png'))
                    if len(files) == 0:
                        return print('There is no data in the folder {}. '.format(key))
                    else:
                        file=files[0]                    
                    self.dict['data'][key][str(nn)] = imageio.imread(file)
                
                if 'roi' in key: # load roi results from xml file
                    files = glob.glob(os.path.join(analysis_folder + 'results', infile + '*.xml'))
                    if len(files) == 0:
                        return print('There is no data in the folder {}. '.format(key))
                    else:
                        file=files[0]                    
                    names, values = Results().extract_results_from_xml(file, protocol='roi', verbosity=3)
                    for kk,name in enumerate(names):
                        if name not in self.dict['data'][key]:
                            self.dict['data'][key][name]= []
                        self.dict['data'][key][name].append(values[kk])  

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

        xy_axes = list(self.dict['data'][key][str(0)].keys())  
        x_axis = self.dict['corrdata']['2Darray'][key][xy_axes[0]]
        I_array = self.dict['corrdata']['2Darray'][key]['I_array']
        if flag_log==1:
            I_array = np.log10(I_array)

        y_idx = np.arange(self.dict['expinfo']['filenumber'])
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
        xy_axes = list(self.dict['data'][key][str(0)].keys())  
        x_axis = self.dict['corrdata']['2Darray'][key][xy_axes[0]]
        I_array = self.dict['corrdata']['2Darray'][key]['I_array']
        if flag_log[1]==1:
            I_array = np.log10(I_array)
        if flag_log[0]==1:
            x_axis = np.log10(x_axis)

        if plot_xrange is None:
            plot_xrange = np.arange(len(x_axis))
    
        y_idx = np.arange(self.dict['expinfo']['filenumber'])
        
        
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
            for key in self.dict['data'].keys():
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
                for key in self.dict['data'].keys():
                    if ('average' in key) or ('linecut' in key):
                        
                        # Data in a 2D dict
                        # 2D array: 'Iq_array' (nrow x ncol)
                        # 1D array: 'q' (len = ncol)
                        # 1D array: 'scan_id', 'time', and other parameters in the 'mdata_list' (len = nrow)
                        #            for series measurement: 'time_series'

                        print('doCorr {}'.format(key))
                        self.dict['corrdata'][corr[0]][key] = {}
                        
                        xy_axes = list(self.dict['data'][key][str(0)].keys())
                        self.dict['corrdata'][corr[0]][key][xy_axes[0]] = self.dict['data'][key][str(0)][xy_axes[0]]

                        qIq = np.array([v for _,one_scan in self.dict['data'][key].items() 
                                          for k,v in one_scan.items()])
                        self.dict['corrdata'][corr[0]][key]['I_array'] = qIq[1::2,:]
                        
                        ## Include corr parameter 
                        # add basic param (always), such as 'file_index'
                        self.dict['corrdata'][corr[0]][key]['file_index'] = np.array(list(self.dict['data'][key].keys()))

                        if self.dict['expinfo']['series_measure']: # Add 'time_series' for series_measurement
                            self.dict['corrdata'][corr[0]][key]['time_series'] = np.arange(self.dict['expinfo']['num_frames']) * self.dict['expinfo']['exposure_period']
                        else:
                            corr_params = ['scan_id', 'time'] # Add 'scan_id' and 'time' for scan/snap measurements
                            for basic_param in corr_params:
                                self.dict['corrdata'][corr[0]][key][basic_param] = np.array(self.dict['expinfo'][basic_param])       

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

            #                 for key in self.dict['data'].keys():
            #                     if ('average' in key) or ('linecut' in key):
                                    
            #                         # Data in a 2D dict
            #                         # 2D array: 'Iq_array' (nrow x ncol)
            #                         # 1D array: 'q' (len = ncol)
            #                         # 1D array: 'scan_id', 'time', and other parameters in the 'mdata_list' (len = nrow)

            #                         self.dict['analysis'][protocol[0]] = {}
            #                         print(key)
            #                         self.dict['analysis'][protocol[0]][key] = {}
                                    
                                    
            #                         xy_axes = list(self.dict['data'][key][0].keys())
                                    
            #                         self.dict['analysis'][protocol[0]][key][xy_axes[0]] = self.dict['data'][key][0][xy_axes[0]] 
                                    
            #                         ## Include basic parameters (always)
            #                         basic_params = ['scan_id', 'time']
            #                         for basic_param in basic_params:
            #                             self.dict['analysis'][protocol[0]][key][basic_param] = np.array(self.dict['expinfo'][basic_param])       
                                
            #                         ## Include corresponding metadata
            #                         for param in self.dict['mdata_list']:
            #                             self.dict['analysis'][protocol[0]][key][param] = np.array(self.dict['metadata'][param])
                                    
            #                         # Iq 2D array
            #                         self.dict['analysis'][protocol[0]][key]['I_array'] = np.array([
            #                             self.dict['data'][key][i][xy_axes[1]] for i in range(self.dict['expinfo']['filenumber'])])        
#                         #load x as the first row
#                         self.dict['analysis'][protocol][key] = self.dict['data'][key]
                        
#                         I_list = []#2D array for data stacking
#                         for ii in self.dict['data'][key]:
#                             #TODO:
#                             I_list.append(self.dict['data'][key])
                            
                            
#                             self.dict['analysis'][protocol]['sequence'][ii] = ii
                            
#                         self.dict['analysis'][protocol]['data'] = np.asarray(I_list)                                
#                         self.dict['analysis'][protocol]['scan_id'] = self.dict['expinfo']['scan_id'] 
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
#                         self.dict['analysis'][pname][pp] = self.dict['data'][pp]
                        
#                     else: #pp = peak_position, peak_width, the same as .xml
#                         self.dict['analysis'][pname]['y_axis'] = self.dict['data'][key]
#                         self.dict['analysis'][pname]['scan_id'] = self.dict['expinfo']['scan_id'] 

#                         for index, data in enumerate(self.dict['data'][pdata]):
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
            for key in self.dict['data'].keys():
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



