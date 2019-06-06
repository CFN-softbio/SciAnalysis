#!/usr/bin/python
# -*- coding: utf-8 -*-
# vi: ts=4 sw=4




################################################################################
#  Short-term settings (specific to a particular user/experiment) can
# be placed in this file. You may instead wish to make a copy of this file in
# the user's data directory, and use that as a working copy.
################################################################################



from ophyd import EpicsSignal
from bluesky.suspenders import SuspendFloor, SuspendCeil


if True:
    # Define suspenders to hold data collection if x-ray
    # beam is not available.
    
    ring_current = EpicsSignal('SR:OPS-BI{DCCT:1}I:Real-I')
    sus = SuspendFloor(ring_current, 100, resume_thresh=101)
    RE.install_suspender(sus)

    #absorber_pos = EpicsSignal( 'XF:11BMB-ES{SM:1-Ax:ArmR}Mtr.RBV')
    #sus_abs_low = SuspendFloor(absorber_pos, -56, resume_thresh=-55)
    #sus_abs_hi = SuspendCeil(absorber_pos, -54, resume_thresh=-55)
    #RE.install_suspender(sus_abs_low)
    #RE.install_suspender(sus_abs_hi)

    #RE.clear_suspenders()


if False:
    # The following shortcuts can be used for unit conversions. For instance,
    # for a motor operating in 'mm' units, one could instead do:
    #     sam.xr( 10*um )
    # To move it by 10 micrometers. HOWEVER, one must be careful if using
    # these conversion parameters, since they make implicit assumptions.
    # For instance, they assume linear axes are all using 'mm' units. Conversely,
    # you will not receive an error if you try to use 'um' for a rotation axis!
    m = 1e3
    cm = 10.0
    mm = 1.0
    um = 1e-3
    nm = 1e-6
    
    inch = 25.4
    pixel = 0.172 # Pilatus
    
    deg = 1.0
    rad = np.degrees(1.0)
    mrad = np.degrees(1e-3)
    urad = np.degrees(1e-6)
    

INTENSITY_EXPECTED_050 = 18800.0
INTENSITY_EXPECTED_025 = INTENSITY_EXPECTED_050*0.5


def get_default_stage():
    return stg


class SampleTSAXS(SampleTSAXS_Generic):
    
    def __init__(self, name, base=None, **md):
        super().__init__(name=name, base=base, **md)
        self.naming_scheme = ['name', 'extra', 'temperature', 'exposure_time']
        
        self.md['exposure_time'] = 30.0


class SampleGISAXS(SampleGISAXS_Generic):
    
    def __init__(self, name, base=None, **md):
        super().__init__(name=name, base=base, **md)
        self.naming_scheme = ['name', 'extra', 'th', 'exposure_time']

 
class Sample(SampleGISAXS):

    def __init__(self, name, base=None, **md):
       super().__init__(name=name, base=base, **md)
       
       #self.naming_scheme = ['name', 'extra', 'clock', 'temperature', 'th', 'exposure_time']
       #self.naming_scheme = ['name', 'extra', 'temperature', 'th', 'exposure_time']
       #self.naming_scheme = ['name', 'extra', 'th', 'exposure_time']
       #self.naming_scheme = ['name', 'extra', 'y', 'th', 'clock', 'exposure_time']
       self.naming_scheme = ['name', 'extra', 'x', 'th', 'exposure_time']
       
       self.md['exposure_time'] = 30.0
       self.incident_angles_default = [0.08, 0.10, 0.12, 0.15, 0.20]


class GIBarCustom(GIBar_long_thermal):
    
    def __init__(self, name='GIBar', base=None, **kwargs):
        
        super().__init__(name=name, base=base, **kwargs)
        
        self._positional_axis = 'x'
        
        # smy = 8.36
        self._axes['y'].origin = 2.06+6.3
        
        self.mark('right edge', x=+152.4)
        self.mark('left edge', x=0)
        self.mark('center', x=76.2, y=0)

        
    def do(self, num_spots=8, x_offset=0.2):
        
        self.alignSamples()
        
        for sample in self.getSamples():
            sample.gotoOrigin(['x','y','th'])
            for ix in range(num_spots):
                sample.xr(x_offset)
                sample.measureIncidentAngles(exposure_time=30)
               


def cooling_off():
    caput('XF:11BMB-ES{IO}AO:1-SP', 0) # Cooling off 
def cooling_on():
    caput('XF:11BMB-ES{IO}AO:1-SP', 5) # Cooling on


cms.SAXS.setCalibration([764, 1680-579], 5.03, [-60, -71]) # 13.5 keV



RE.md['experiment_group'] = 'User group (University)'
RE.md['experiment_alias_directory'] = '/nsls2/xf11bm/data/2019_X/UShell/'
RE.md['experiment_proposal_number'] = 'TBD'
RE.md['experiment_SAF_number'] = 'TBD'
RE.md['experiment_user'] = 'TBD'
RE.md['experiment_type'] = 'GISAXS' # TSAXS, GISAXS, GIWAXS, etc.
RE.md['experiment_project'] = 'TBD'




if True:
    
    md = {
    'owner' : 'Name' ,
    'series' : 'TBD' ,
    'substrate' : 'Si' ,
    'substrate_thickness' : 0.5 ,
    }

    
    hol = GIBarCustom(base=stg) # P, I, D = 60, 5, 5 (Long heater bar)

    eo = 0
    hol.addSampleSlotPosition( Sample('run14_DG-HPL2Days_H3wt%_1_C-Si_Substrate', **md), 1, 9.80+eo)
    hol.addSampleSlotPosition( Sample('run14_DG-HPL2Days_H3wt%_2_C-Si_Substrate', **md), 2, 22.68+eo)
  

  
    
if False:
    
    md = {
    'owner' : 'Name' ,
    'series' : 'TBD' ,
    }
    RE.md['experiment_type'] = 'TSAXS'

    
    hol = CapillaryHolderCustom(base=stg) # P, I, D = 60, 2, 3 (Capillary holder heater, with cooling loop)

    hol.addSampleSlot( Sample('run16_BlackOpsRamp_1', **md), 8)
    hol.addSampleSlot( Sample('run16_BlackOpsRamp_2', **md), 9)
   
   
        

    
# To reload this script in your bsui session:
# 1. Edit this file
# 2. Save this file
# 3. Run this command in your bsui session:
#   %run -i /nsls2/xf11bm/data/2019_X/UShell/user.py
