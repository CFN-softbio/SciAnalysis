#!/usr/bin/python
# -*- coding: utf-8 -*-
# vi: ts=4 sw=4


# Load BLEX
#import sys, os
#BLEX_PATH='/GPFS/xf12id1/analysis/CFN/BLEX/SMI_2018C3/'
#BLEX_PATH in sys.path or sys.path.append(BLEX_PATH)
#from sample import *


# %run -i /GPFS/xf12id1/analysis/CFN/BLEX/SMI_2018C3/81-beam.py
# %run -i /GPFS/xf12id1/analysis/CFN/BLEX/SMI_2018C3/91-fit_scan.py
# %run -i /GPFS/xf12id1/analysis/CFN/BLEX/SMI_2018C3/94-sample.py



# Sample types
########################################

class SampleGISAXS(SampleGISAXS_Generic):
    
    def __init__(self, name, base=None, **md):
        super().__init__(name=name, base=base, **md)
        self.naming_scheme = ['name', 'extra', 'th', 'exposure_time']
        self.incident_angles_default = [0.08, 0.10, 0.12, 0.15, 0.20]
        


        
class Sample(SampleGISAXS):
    pass





RE.md['experiment_proposal_number'] = 302843
RE.md['experiment_SAF_number'] = 303365
RE.md['experiment_group'] = 'Chris Murray (U.Penn.)'
RE.md['experiment_user'] = 'TBD'
RE.md['experiment_project'] = 'nanoparticles'
RE.md['experiment_alias_directory'] = '/GPFS/xf12id1/data/CFN/2018_3/KYager/'
RE.md['experiment_type'] = 'GISMAXS'
RE.md['beam_intensity_expected'] = 19300

smi.bsx_pos = 1.2 # Centered
smi.bsx_pos = 1.45 # Reduced scattering from bs
smi.detectors_measure = [pil1M, rayonix]
#smi.detselect(pil1M)
smi.detselect([pil1M, rayonix])
smi.SAXS.setCalibration([457.0, 569.0], 5.300, [-5.4957, -60.0])

print('\n\n\nReminders:')
print('    Define your detectors using, e.g.: detselect(pilatus2M)')
print('    Reload your user-specific script, e.g.: %run -i /GPFS/xf11bm/data/2017_2/user_group/user.py')
print('\n')



# User samples
########################################        
if True:
    # For testing (and as examples...)
    # %run -i /opt/ipython_profiles/profile_collection/startup/98-user.py
    
    hol = GIBar(base=stg)
    hol.addSampleSlotPosition( Sample('test_sample_01'), 1, 56.1 )
    hol.addSampleSlotPosition( Sample('test_sample_02'), 3, 56.1+10 )
    
    sam = hol.getSample(1)    
    
    
    

# Note for user
########################################
# %run -i /GPFS/xf12id1/data/CFN/2018_3/KYager/user.py







# Development notes:
########################################
# SMI bsui ipython startup folder:
# /home/xf12id/.ipython/profile_collection/startup
#
# %run -i /GPFS/xf12id1/data/CFN/2018_3/KYager/user.py