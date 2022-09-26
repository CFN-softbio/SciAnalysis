# Global settings which influence how SciAnalysis operates.


SUPPRESS_EXCEPTIONS = False
#SUPPRESS_EXCEPTIONS = True # Suppress Python exceptions (errors). This allows the script to keep running even if there is an error processing one particular file.


#MATPLOTLIB_BACKEND = None # Leave as default
MATPLOTLIB_BACKEND = 'Agg' # For 'headless' plotting (e.g. over an SSH connection, or to avoid bugs with joblib parallelization)


DEFAULT_SAVE_RESULTS = ['xml', 'sql', 'plots', 'txt']
#DEFAULT_SAVE_RESULTS = ['xml', 'sql', 'plots', 'txt', 'npy', 'npz', 'pkl', 'hdf5']
