import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import pandas as pd
import glob, os, time
import databroker
import imageio

from SciAnalysis import tools
from SciAnalysis.XSAnalysis.Data import *
from SciAnalysis.XSAnalysis import Protocols
from SciAnalysis.Result import *

#db = databroker.DataBroker.named('cms')
cat = databroker.catalog['cms']

class experiment():
    def __init__(self, name, folder=None, det='saxs'):
        
        #self.name = name
        self.det = det
        
        if folder==None:
            self.folder = os.getcwd() 
#             self.folder = '/nsls2/data/cms/legacy/xf11bm/data/2022_1/'+'user'+det
        else:
            self.folder = folder
        
        self.dict = {'exp': name,
                     'type': det,
                    }                     
        self.dict['expinfo'] = {}
        self.dict['expinfo']['filename'] = []
        self.dict['expinfo']['time'] = []
        self.dict['expinfo']['clock'] = []
        self.dict['expinfo']['scan_id'] = []
        self.dict['expinfo']['uid'] = []
        self.dict['expinfo']['filenumber']=0 # total number of files

        self.dict['metadata'] = {}
        self.dict['mdata_list'] = []
        self.dict['analysis'] = {}
        self.dict['protocol'] = []
    
    def defFiles(self, fn = None, scanid=None, uid=None, stitched=False, burstmode=False):
        #define the files in the experiment
        #search raw tiff, return the list of scanid or uid or filenames (RL_t0, RL_t1, ...) 
        # and also look up metadata with the scanid
#         keys = ['sample_x', 'sample_temperature', 'scan_id' ] 
        
        if fn==None:
            fn = self.name
            
        # define the source_dir
        if stitched == False:
            source_dir = self.folder + self.det + '/raw/'
        else:
            source_dir = self.folder + self.det + '/stitched/'
        
        #
        if uid != None:
            for uidt in uid:
                h = cat[uidt]
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
        self.dict['expinfo']['filenumber'] = len(self.dict['expinfo']['uid'])
                
    def defFiles_query(self, cycle=None, SAF=None, fn=None, timerange=None, folder=None, scanid=None):

        if folder is None:
            folder = self.folder
            
        while folder[-1]=='/': folder = folder[:-1]

        query = {
             'experiment_alias_directory': folder,
#              'scan_id':{'$gte': scanid[0], '$lte': scanid[-1]},
#              'sample_name':fn,
#              'time_range'
#              'experiment_SAF_number': SAF,
                }
        
        if scanid is not None:
            query['scan_id'] = {'$gte': scanid[0], '$lte': scanid[-1]}
        if fn is not None:
            query['sample_name'] = fn
        if SAF is not None:
            query['experiment_SAF_number'] = SAF
        if timerange is not None:
            query['time_range'] = time_range
        if cycle is not None:
            query['experiment_cycle'] = cycle

        ## To save the ROI intensity from the detector
        if 'data' not in self.dict:
            self.dict['data'] = {} # create a dict for data loading
        self.dict['data']['det']= {}
        for i in range(1,5):
            self.dict['data']['det'][f'roi{i}']=[]
    
        if self.dict['type'] == 'saxs':
            detector = 'pilatus2M'
        if self.dict['type'] == 'waxs':
            detector = 'pilatus800'
        if self.dict['type'] == 'maxs':
            detector = 'pilatus8002'
            
#         print(query)
        t0 = time.time()
        results = cat.search(query)
#         results = results[::-1]
        
        for uid in results:
            
            print(uid)
#             infile_scanid = infile.split('_'+self.det)[0].split('_')[-1]
            h = results[uid]
            self.dict['expinfo']['filename'].append(h.metadata['start']['filename'])
            self.dict['expinfo']['time'].append(h.metadata['start']['time']) #linux time
            self.dict['expinfo']['clock'].append(h.metadata['start']['sample_clock'])
            self.dict['expinfo']['scan_id'].append(h.metadata['start']['scan_id'])
            self.dict['expinfo']['uid'].append(h.metadata['start']['uid'])
            
            
            #### primary.read() is very slow
#             print(h.primary['data']['pilatus800_stats1_total'])
#             print(h.primary.read()['pilatus800_stats1_total'].values[0])
#             for i in range(1,5):
#                 self.dict['data']['det'][f'roi{i}'].append(h.primary.read()[f'{detector}_stats{i}_total'].values[0])
        self.dict['expinfo']['filenumber'] = len(self.dict['expinfo']['uid'])
        print('Loaded {} files, took {}s.'.format(self.dict['expinfo']['filenumber'], time.time()-t0))

        return results

#     def defMetadata(self, uid, keys=['scan_id']):
#         #load metadata
#         h = cat[uid]
#         for key in keys:
#             self.dict['metadata'][key] = h.start['uid'][key]
    
    def loadMetadata(self, keys=None):
        if keys is None:
            keys = self.dict['mdata_list']
        
        for key in keys:
            self.dict['metadata'][key] = []
        #
        for uid in self.dict['expinfo']['uid']:
#         uid = self.dict['expinfo']['uid'][0]
            h = cat[uid]
            for key in keys:
                self.dict['metadata'][key].append(h.metadata['start'][key])

    #TODO: change fn to self.dict
    def listMetadata(self, scanid=None, uid=None, verbosity=3):
        '''
        list of the keys in Metadata to input
        '''
    
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

    #load the reducted data (SciAanalysis) files via the key 
    #key = 'circurlar_average'  or 'linecut_qr-q=0.1'
    #filename full path: 'circular_average/RL_t0.dat',  'qr_image/RL_t0.npz'
    #input: key = 'circular_average', filename = 'RL_t0.dat'

    
    
    def loadSciAnalysisData(self, keys=None):
        
        if 'data' not in self.dict:
            self.dict['data'] = {} # create a dict for data loading
        
        analysis_folder = self.folder + '/' + self.det + '/analysis/'
        if keys is None:
            folders = glob.glob(analysis_folder + '/*/')
            keys = []
            for folder in folders:
                key = folder.split('/')[-2]
                keys.append(key)
        
#         print(folders)
#         print(keys)

        for nn, infile in enumerate(self.dict['expinfo']['filename']):

            for key in keys:
                
                if key not in self.dict['data']:
                    self.dict['data'][key] = {}

                if 'average' in key:
                    file = glob.glob(os.path.join(analysis_folder + key + '/' + infile + '*.dat'))[0]
                    
                    headers = pd.read_csv(file, delim_whitespace=True, nrows=0).columns[1:]
                    dat = pd.read_csv(file, delim_whitespace=True, header=None, skiprows=1, names=headers)
                    
                    self.dict['data'][key][nn] = {}
                    self.dict['data'][key][nn][headers[0]] = dat[headers[0]].values
                    self.dict['data'][key][nn][headers[2]] = dat[headers[2]].values
                    

                if 'linecut' in key:
                    file = glob.glob(os.path.join(analysis_folder + key + '/' + infile + '*.dat'))[0]
                    headers = pd.read_csv(file, delim_whitespace=True, nrows=0).columns[1:]
                    dat = pd.read_csv(file, delim_whitespace=True, header=None, skiprows=1, names=headers)
                    
                    self.dict['data'][key][nn] = {}
                    self.dict['data'][key][nn][headers[0]] = dat[headers[0]].values
                    self.dict['data'][key][nn][headers[1]] = dat[headers[1]].values
                    
                if 'image' in key:
                    file = glob.glob(os.path.join(analysis_folder + key + '/' + infile + '*.png'))[0]
                    self.dict['data'][key][nn] = imageio.imread(file)
                
                if 'roi' in key: # load roi results from xml file
                    file = glob.glob(os.path.join(analysis_folder + 'results' + '/' + infile + '*.xml'))[0]
                    names, values = Results().extract_results_from_xml(file, protocol='roi', verbosity=3)
                    for kk,name in enumerate(names):
                        if name not in self.dict['data'][key]:
                            self.dict['data'][key][name]= []
                        self.dict['data'][key][name].append(values[kk])                       
                        
        
                    


    #TODO
    def doPlot(self, protocol, key):
        pass
        
        

    #TODO
    #plotting all heat maps
    def plotHeatmap(self, key, y_axes=None, flag_log = 1, plot_range=np.arange(0,100)):
        '''
        plot the heatmap for the 2D arrays
        '''
        xy_axes = list(self.dict['data'][key][0].keys())  
        x_axis = self.dict['analysis']['2Darray'][key][xy_axes[0]]
        I_array = self.dict['analysis']['2Darray'][key]['I_array']
        if flag_log==1:
            I_array = np.log10(I_array)

        y_idx = np.arange(self.dict['expinfo']['filenumber'])
        
        if y_axes is not None:
            y_axis = self.dict['analysis']['2Darray'][key][y_axes[0]]    
        
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
        plt.pcolormesh(x_axis[plot_range], y_idx, I_array[:,plot_range])

        plt.xlabel('{}.(A^-1 or azi angle)'.format(xy_axes[0]))
        plt.ylabel('idx')
        if y_axes is not None:
            plt.yticks(y_idx, y_axis)
            plt.ylabel(y_axes[0])
        plt.colorbar()

        ## plot the 1D
        y_params = list(self.dict['analysis']['2Darray'][key].keys())
        y_params.remove(xy_axes[0])
        y_params.remove('I_array')
        
        cmap = plt.get_cmap('jet',len(y_params))
        
        fig, axs = plt.subplots(len(y_params), 1, sharex=True)
        # Remove horizontal space between axes
        fig.subplots_adjust(hspace=0)

        for kk, ax in enumerate(axs.flat):
            ax.plot(y_idx, self.dict['analysis']['2Darray'][key][y_params[kk]], color = cmap(kk), marker='o')
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


    def plotWaterfall(self, key, y_axes=None, gridlines=True, flag_log = [0, 0], plot_range=np.arange(0,100)):
        '''
        plot waterfall for the 2D arrays
        '''
        xy_axes = list(self.dict['data'][key][0].keys())  
        x_axis = self.dict['analysis']['2Darray'][key][xy_axes[0]]
        I_array = self.dict['analysis']['2Darray'][key]['I_array']
        if flag_log[1]==1:
            I_array = np.log10(I_array)
        if flag_log[0]==1:
            x_axis = np.log10(x_axis)
        y_idx = np.arange(self.dict['expinfo']['filenumber'])
        
        if y_axes is not None:
            y_axis = self.dict['analysis']['2Darray'][key][y_axes[0]]    
        
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
        spacing = np.mean(I_array)*2
        N = len(I_array)
        for nn, curve in enumerate(I_array):
            plt.plot(x_axis[plot_range], curve[plot_range] + spacing*nn, color=[0.1+nn/N, 1-nn/N, 0.8-nn/N/2])

        plt.xlabel('{}.(A^-1 or azi angle)'.format(xy_axes[0]))
        if gridlines:
            plt.grid()

    
    
    def addProtocol(self, fitting=False, verbosity=3):
        
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
            self.dict[protocol].append([data, medata, fitting])
        else:
            self.dict[protocol].append([data, medata])
        
    def listProtocol(self, verbosity=3):        
        for protocol in self.dict[protocol]:
            print(protocol)
        
    def doAnalysis(self, protocols=None):
        '''
            input : 
                protocols = [['2Darray'],
                             ['circular_average','temperature'],
                             # For peak q0, plot peak intensity (or width etc) vs. sample position x
                                ['circular_average_fitting',  'sample_x', 'peak_q',] 
                             # 2D mapping in smx and smy : For peak q0, plot peak intensity (or width etc) vs (x,y)
                            ]

                2Darray 
                average * tempreature: 2D array 
                circular_average_fitting*sample_x*peak_q]

            results:
                self.dict['analysis']['circular_average__temperature'] #save 2D array, (optional export PNG)
                self.dict['analysis']['peak_q__sample_x']

        '''    
        if protocols == None:
            protocols = self.dict['protocol']
            
        for protocol in protocols:
            
            #2D array
            if protocol == ['2Darray']:
                for key in self.dict['data'].keys():
                    if ('average' in key) or ('linecut' in key):
                        
                        # Data in a 2D dict
                        # 2D array: 'Iq_array' (nrow x ncol)
                        # 1D array: 'q' (len = ncol)
                        # 1D array: 'scan_id', 'time', and other parameters in the 'mdata_list' (len = nrow)

                        self.dict['analysis'][protocol[0]] = {}
                        print(key)
                        self.dict['analysis'][protocol[0]][key] = {}
                        
                        
                        xy_axes = list(self.dict['data'][key][0].keys())
                        
                        self.dict['analysis'][protocol[0]][key][xy_axes[0]] = self.dict['data'][key][0][xy_axes[0]] 
                        
                        ## Include basic parameters (always)
                        basic_params = ['scan_id', 'time']
                        for basic_param in basic_params:
                            self.dict['analysis'][protocol[0]][key][basic_param] = np.array(self.dict['expinfo'][basic_param])       
                    
                        ## Include corresponding metadata
                        for param in self.dict['mdata_list']:
                            self.dict['analysis'][protocol[0]][key][param] = np.array(self.dict['metadata'][param])
                        
                        # Iq 2D array
                        self.dict['analysis'][protocol[0]][key]['I_array'] = np.array([
                            self.dict['data'][key][i][xy_axes[1]] for i in range(self.dict['expinfo']['filenumber'])])                    
                            

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



