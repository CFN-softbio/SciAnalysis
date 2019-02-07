#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import sys # To get commandline arguments

#SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
#SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)
om SciAnalysis import tools
#from SciAnalysis.XSAnalysis.Data import *
from SciAnalysis.ImAnalysis.Data import *



if len(sys.argv)>1:
    # Use commandline argument
    source_dir = ''
    output_dir = './'
    infile = sys.argv[1]
    
else:
    # Use hardcoded file
    root_dir = '/home/kyager/ppl/AtikurRahman/layered_BCPs_01/'
    data_dir = 'L-36/'
    #data_dir = 'normal/'
    
    protocol = 'local_avg_realspace'
    #protocol = 'thumbnails'
    source_dir = os.path.join( root_dir, 'raw', data_dir )
    output_dir = os.path.join( root_dir, 'ImAnalysis', data_dir )
    
    tools.make_dir(output_dir)
    
    #infile = 'L-36-3cyc_L-74-3-cyc_O2-at-end_1_q01.tif'
    
    import glob
    #infiles = glob.glob(source_dir + '*.tif')
    infiles = glob.glob(source_dir + 'L-36-3cyc_L-74-3-cyc_O2-at-end_1_q01.tif')
    
    #infiles = glob.glob(source_dir + '6.tif')
    #infiles = glob.glob(source_dir + '2012-09-25 Lam Blend01-2p5KRPM_q02.tif')
    
    infiles.sort()
    


if True:
    
    from SciAnalysis.ImAnalysis import Protocols
    process = Protocols.ProcessorIm()
    #protocols = [ Protocols.thumbnails(crop=1.3) ]
    #protocols = [ Protocols.circular_average(ylog=True) ]
    #protocols = [ Protocols.thumbnails() ]

    load_args = { 'crop_bottom' : 124, # pixels
                    'scale' : 500.0/353, # nm/pixel
                    }
    
    run_args = { 'verbosity' : 3,
                    'local_partition_image_size' : 30, # pixels
                    'local_partition_step' : 5.0, # relative to image_size
                    #'blur' : 1.0, # pixels
                    'q0' : 0.18, # nm^-1
                    'dq' : 0.18*0.6, # nm^-1
                    'symmetry' : 4,
                    #'rcParams': {'axes.labelsize': 45,},
                    }
                
    
    protocols = [ Protocols.fft(**run_args) ]
    #protocols = [ Protocols.local_avg_realspace(**run_args) ]
    
    
    
    process.run(infiles, protocols, output_dir=output_dir, force=True, load_args=load_args)
    
        
    
    
    
exit()


if False:
    
    root_dir = '/home/kyager/ppl/EmilyCranston/Kevin_DeFrance/2015_10Oct_18-SAXS_all/Hydrogels/Static tests & cals/'
    source_dir = root_dir
    
    protocol = 'circular_average'
    output_dir = root_dir + '/SciAnalysis/' + protocol + '/'
    tools.make_dir(output_dir)
    
    
    calibration = Calibration(wavelength_A=1.5418) # 8.04 keV (Cu K-alpha e.g. Nanostar)
    calibration.set_image_size(2048)
    calibration.set_pixel_size(width_mm=140.0)
    #calibration.set_distance(1.11508)
    calibration.set_distance(1.15508)
    calibration.set_beam_position(1019.0, 1029.5)
    
    mask_dir = '/home/kyager/current/code/SciAnalysis/main/SciAnalysis/XSAnalysis/masks/'
    mask = Mask(mask_dir+'Bruker_Nanostar_SAXS-2015Oct-mask.png')
    
    data = Data2DScattering(calibration=calibration, mask=mask)
    #data.load_BrukerASCII(source_dir+'AgBH_kineticCal.dat')
    data.load_tiff(source_dir+'AgBH_kineticCal.tiff')
    #data.blur(2.0)
    
    data._linecut_test(0.1076, 0.001)
    data.set_z_display( [None, None, 'gamma', 0.3] )
    data.plot(show=True)
    #data.plot_image('out.png')
    
    if False:
        line = data.circular_average_q_bin(error=True)
        line.smooth(2.0, bins=10)
        line.plot(show=True, plot_range=[0, 0.27, 0, None], ylog=False, error=False, error_band=False)
        #line.save_data('out_q.data')
    
    
    
    

if False:
    
    # Output a plot of a mask
    #m = Mask('CHX_Eiger4M_blemish-mask.hd5')
    m = Mask('ffield_threshed_EIGER4M.hd5')
    
    idx = abs(m.data - 1.0) > 0.05
    m.data[idx] = 0
    m.data[~idx] = 1
    
    data = Data2DScattering()
    data.data = m.data
    data.plot(show=True)
    data.plot_image(save='./a.png')


if False:
    
    mask_dir = '/home/kyager/current/code/SciAnalysis/main/SciAnalysis/XSAnalysis/masks/'
    mask = Mask()
    mask.load(mask_dir+'Eiger4M_all_gaps-mask.png')
    #mask.load(mask_dir+'CHX_Eiger4M-bad_flatfield_10percent-mask.png')
    mask.load(mask_dir+'CHX_Eiger4M-bad_flatfield_05percent-mask.png')
    mask.load(mask_dir+'CHX_pipe-2015Oct-mask.png')
    mask.load(mask_dir+'CHX_bs_streaks-2015Oct-mask.png')
    
    
    
    
    data = Data2DScattering(source_dir+infile, mask=mask)
    #data.threshold_pixels(4294967295-1) # Eiger inter-module gaps
    #data.blur(0.7)
    
    data.set_z_display( [None, None, 'gamma', 0.3] )
    data.set_z_display( [1, 64, 'gamma', 0.3] )

    #data.origin = [839, 1826]
    #data.set_z_display( [0, 4e6, 'r', 2.0] )

    kwargs = { 'size':10 }
    outfile = infile[:-3]+'.png'
    #data.plot(save=output_dir+outfile, show=True)
    data.plot(show=True, **kwargs)
    
    #data.plot_image(save=output_dir+outfile)
    #data.plot_image(save='./'+outfile)


if False:
    import glob
    infiles = glob.glob(source_dir+'*_master.h5')

    tools.make_dir(output_dir)
    
    
    # TODO: Make this a processing function
    infiles.sort()
    for infile in infiles:
        
        # TODO: Check if file exists (and 'force' kwarg)
       
       
        print('Running %s on %s...'%(protocol, infile))
        
        try:

            data = Data2DScattering(infile)

            data.threshold_pixels(4294967295-1) # Eiger inter-module gaps

            outfile = output_dir+infile[len(source_dir):-3]+'.png'
            #print('\n\n')
            #print( outfile )
            #print('\n\n')
        
            data.plot(save=outfile, show=False, **kwargs)
            
        except KeyboardInterrupt:
            raise
        
        except:
            print('    Error when processing %s on %s'%(protocol, infile))
            #print(sys.exc_info()[0])
            
            
