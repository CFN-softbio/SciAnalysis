##Octo 20, 2019 @CFN 
## This code is essentially to save dictionary as h5 file and load h5 file as dictionary
## Two important functions are    dicttoh5,  h5todict
## Examples:
# save dictionary as h5 file: dicttoh5( dictionary,  output_filename, h5path= key,  mode='a'  )        
# load h5 as dictionary:  h5todict( input_filename )
## This code is developed based on silx


# Note: Imports are kept to a minimum to prevent errors on systems that may not have
# all packages installed.
#import numpy as np
#from PIL import Image
#import shutil,glob, os
#from os import listdir
#from os.path import isfile, join
#import h5py, sys, enum
#import pandas as pds
#import collections
#import logging
#logger = logging.getLogger(__name__)


# Note: Currently SciAnalysis only uses the dicttoh5 function, which only requires
# the imports below. If additional functions of this python file are used, more
# of the import commands (above) will need to be turned on.
import numpy as np
import h5py, sys, enum
import collections

## a global string type     
string_types = (basestring,) if sys.version_info[0] == 2 else (str,)
def dicttoh5(treedict, h5file, h5path='/',mode="w", overwrite_data=False, create_dataset_args=None):
    '''Save dictionary as a h5 file
    Input:
        treedict:  dictionary, could be nested format
        h5file: the filename the output h5 file
        h5path,  group in the output h5 file, default '/' (the root key)
        mode:   Can be 
                    ``"r+"`` (read/write, file must exist)
                    ``"w"`` (write, existing file is lost)
                    ``"w-"`` (write, fail ifexists)
                    ``"a"`` (read/write if exists, create otherwise).
        overwrite_data: bool, 
                if True, existing groups and datasets can be overwritten
                if False, they are skipped. This parameter is only relevant if ``h5file_mode`` is ``"r+"`` or ``"a"``
        create_dataset_args: dictionary, specify filters andcompression parameters
    Output:
        None
    Example1:
        d = { 'data': np.arange(10), 'label':'data', 'results': {  'collected time': '2019/10/10', 'fit_res': { 'p0': 10 } } }
        dicttoh5( d, 'test.h5', h5path='/', mode='a')
        
        --> if do h5todict( 'test.h5' )
            get: -->
                    {'data': array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
                      'label': 'data',                      
                      'results': {'collected time': '2019/10/10', 'fit_res': {'p0': array(10)}}}
                
    Example2:
        d = { 'data': np.arange(10), 'label':'data', 'results': {  'collected time': '2019/10/10', 'fit_res': { 'p0': 10 } } }
        dicttoh5( d, 'test.h5', h5path='/cms/', mode='a')
        
        --> if do h5todict( 'test.h5' )
            get: -->
                    {'cms': 
                        {'data': array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
                          'label': 'data',                      
                          'results': {'collected time': '2019/10/10', 'fit_res': {'p0': array(10)}}}    
                    }
    
    '''   
    if create_dataset_args is None:
        create_dataset_args = {'compression': "gzip", 'shuffle': True,  'fletcher32': True}    
    compssion,shuffle,fletcher32 = ( create_dataset_args['compression'], 
                                     create_dataset_args['shuffle'], 
                                     create_dataset_args['fletcher32'] 
                                   )
    if not h5path.endswith("/"):
        h5path += "/"
    with _SafeH5FileWrite(h5file, mode=mode) as h5f:
        for key in treedict:
            #print(key)
            if isinstance(treedict[key], dict) and len(treedict[key]):
                # non-empty group: recurse
                dicttoh5(treedict[key], h5f, h5path + key,overwrite_data=overwrite_data, create_dataset_args=create_dataset_args) 
            elif treedict[key] is None or (isinstance(treedict[key], dict) and not len(treedict[key])):
                if (h5path + key) in h5f:
                    if overwrite_data is True:
                        del h5f[h5path + key]
                    else:
                        logger.warning('key (%s) already exists. '  'Not overwriting.' % (h5path + key))
                        continue                
                h5f.create_group(h5path + key)
            else:
                ds = _prepare_hdf5_dataset(treedict[key])               
                if ds.shape == () or create_dataset_args is None:
                    if h5path + key in h5f:
                        if overwrite_data is True:
                            del h5f[h5path + key]
                        else:
                            logger.warning('key (%s) already exists. ' 'Not overwriting.' % (h5path + key))
                            continue
                    h5f.create_dataset(h5path + key, data=ds)
                else:
                    if h5path + key in h5f:
                        if overwrite_data is True:
                            del h5f[h5path + key]
                        else:
                            logger.warning('key (%s) already exists. ' 'Not overwriting.' % (h5path + key))
                            continue
                    h5f.create_dataset(h5path + key, data=ds, **create_dataset_args)                     

def h5todict(h5file, path="/", exclude_names=None, asarray=True):
    '''Load h5 to a dictionary 
    Input:         
        h5file: the filename the input h5 file
        path,  group in the input h5 file, default '/' (the root key)
        exclude_names:  will not load that exclude_names group
        asarray: bool, 
                if True, (default) to read scalar as arrays
                if False, to read them as scalar
                if False, they are skipped. This parameter is only relevant if ``h5file_mode`` is ``"r+"`` or ``"a"``         
    Output:
        dictionary
    Example:
        d = { 'data': np.arange(10), 'label':'data', 'results': {  'collected time': '2019/10/10', 'fit_res': { 'p0': 10 } } }
        dicttoh5( d, 'test.h5', h5path='/cms/', mode='a')
        
        --> if do h5todict( 'test.h5' )
            get: -->
                    {'cms': 
                        {'data': array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
                          'label': 'data',                      
                          'results': {'collected time': '2019/10/10', 'fit_res': {'p0': array(10)}}}    
                    }
    
    '''     
    with _SafeH5FileRead(h5file) as h5f:
        ddict = {}        
        for key in h5f[path]:
            if _name_contains_string_in_list(key, exclude_names):
                continue
            if is_group(h5f[path + "/" + key]):
                ddict[key] = h5todict(h5f,path + "/" + key, exclude_names=exclude_names, asarray=asarray)
            else:                
                data = h5f[path + "/" + key][()]
                if asarray:  # Convert HDF5 dataset to numpy array
                    if isinstance( data , bytes):
                        data = bstring_to_string( data )
                    elif ( isinstance( data, np.ndarray )  or isinstance( data, list )  ):
                        if isinstance( data[0] , bytes):
                            data = bstring_to_string( data )  
                    elif ( isinstance( data, float )  or isinstance( data, int )   or isinstance( data, str ) ):
                        pass                                                    
                    else:
                        data = np.array(data, copy=False)                      
                ddict[key] = data
    return ddict

def bstring_to_string( bstring ):
    '''Y.G., Dev@CFN Nov 20, 2019 convert a btring to string
     
    Parameters
    ----------
        bstring: bstring or list of bstring 
    Returns
    -------  
        string:    
    '''
    s =  np.array( bstring )
    if not len(s.shape):
        s=s.reshape( 1, )
        return s[0].decode('utf-8')
    else:
        return np.char.decode( s )    

def _prepare_hdf5_dataset(array_like):
    '''A python object into a numpy array in a HDF5 friendly format.'''
    if isinstance(array_like, string_types):
        array_like = np.string_(array_like)    
    if not isinstance(array_like, (np.ndarray, np.string_)):
        array = np.array(array_like)
    else:
        array = array_like    
    if not isinstance(array, np.string_):
        data_kind = array.dtype.kind
        if data_kind.lower() in ["s", "u"]:
            array = np.asarray(array, dtype=np.string_)
    return array

class _SafeH5FileWrite(object):
    ''' A class for safe write h5 file, close the file on exiting.'''
    def __init__(self, h5file, mode="w"):
        self.raw_h5file = h5file
        self.mode = mode
    def __enter__(self):
        if not isinstance(self.raw_h5file, h5py.File):
            self.h5file = h5py.File(self.raw_h5file, self.mode)
            self.close_when_finished = True
        else:
            self.h5file = self.raw_h5file
            self.close_when_finished = False
        return self.h5file
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.close_when_finished:
            self.h5file.close()

class _SafeH5FileRead(object):
    ''' A class for safe read h5 file, close the file on exiting.'''
    def __init__(self, h5file):
        self.raw_h5file = h5file 
    def __enter__(self):
        if not is_file(self.raw_h5file):
            self.h5file = h5py.File(self.raw_h5file, "r") 
            self.close_when_finished = True
        else:
            self.h5file = self.raw_h5file
            self.close_when_finished = False 
        return self.h5file
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.close_when_finished:
            self.h5file.close()            
               
def _name_contains_string_in_list(name, strlist):
    if strlist is None:
        return False
    for filter_str in strlist:
        if filter_str in name:
            return True
    return False

class H5Type(enum.Enum):    
    '''Identify a set of HDF5 format'''
    DATASET = 1
    GROUP = 2
    FILE = 3
    SOFT_LINK = 4
    EXTERNAL_LINK = 5
    HARD_LINK = 6

def _get_classes_type():   
    '''Returns a mapping between Python classes and HDF5 concepts.    
    '''        
    _CLASSES_TYPE = collections.OrderedDict()
    _CLASSES_TYPE[h5py.Dataset] = H5Type.DATASET
    _CLASSES_TYPE[h5py.File] = H5Type.FILE
    _CLASSES_TYPE[h5py.Group] = H5Type.GROUP
    _CLASSES_TYPE[h5py.SoftLink] = H5Type.SOFT_LINK
    _CLASSES_TYPE[h5py.HardLink] = H5Type.HARD_LINK
    _CLASSES_TYPE[h5py.ExternalLink] = H5Type.EXTERNAL_LINK
    return _CLASSES_TYPE

def get_h5_class(obj=None, class_=None):
    '''Returns the HDF5 type relative to the object or to the class.'''
    if class_ is None:
        class_ = obj.__class__
    classes = _get_classes_type()
    t = classes.get(class_, None)
    if t is not None:
        return t
    if obj is not None:
        if hasattr(obj, "h5_class"):
            return obj.h5_class
    for referencedClass_, type_ in classes.items():
        if issubclass(class_, referencedClass_):
            classes[class_] = type_
            return type_
    classes[class_] = None
    return None

def h5type_to_h5py_class(type_):
    if type_ == H5Type.FILE:
        return h5py.File
    if type_ == H5Type.GROUP:
        return h5py.Group
    if type_ == H5Type.DATASET:
        return h5py.Dataset
    if type_ == H5Type.SOFT_LINK:
        return h5py.SoftLink
    if type_ == H5Type.HARD_LINK:
        return h5py.HardLink
    if type_ == H5Type.EXTERNAL_LINK:
        return h5py.ExternalLink
    return None

def get_h5py_class(obj):
    if hasattr(obj, "h5py_class"):
        return obj.h5py_class
    type_ = get_h5_class(obj)
    return h5type_to_h5py_class(type_)

def is_group(obj):
    t = get_h5_class(obj)
    return t in [H5Type.GROUP, H5Type.FILE]

def is_file(obj):
    t = get_h5_class(obj)
    return t == H5Type.FILE 
