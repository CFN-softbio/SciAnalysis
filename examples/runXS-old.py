#!/usr/bin/python3
# -*- coding: utf-8 -*-

import sys # To get commandline arguments

#SciAnalysis_PATH='/home/kyager/current/code/SciAnalysis/main/'
#SciAnalysis_PATH in sys.path or sys.path.append(SciAnalysis_PATH)

from SciAnalysis import tools
from SciAnalysis.XSAnalysis.Data import *

if len(sys.argv)>1:
    # Use commandline argument
    source_dir = ''
    output_dir = './'
    infile = sys.argv[1]
    
else:
    # Use hardcoded file
    root_dir = '/media/extend/CHX/'
    data_dir = 'LDRD_Meso_correlation/2015_10Oct-ref_patterns/21/'
    
    protocol = 'view'
    source_dir = root_dir+'data/'+data_dir
    output_dir = root_dir+'analysis/'+data_dir+'/'+protocol+'/'
    
    tools.make_dir(output_dir)
    
    infile = 'AgBH_22_master.h5'
    infile = 'Coral_porous_glass_0.3s_30_master.h5'
    infile = 'chip411_found_star_x-0.356_y13.468_922_master.h5'




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
        
        
        
if True:
    root_dir = '/media/extend/CHX/'
    data_dir = 'LDRD_Meso_correlation/2015_10Oct-ref_patterns/24/'
    source_dir = root_dir+'data/'+data_dir
    output_dir = root_dir+'analysis/'+data_dir

    import glob
    infiles = glob.glob(source_dir+'FL*_master.h5')
    
    infiles.sort()
    
    
    # Experimental setup
    calibration = Calibration()
    calibration.set_energy(8.8984) # CHX
    calibration.set_image_size(2070, 2167) # Eiger 4M
    calibration.set_pixel_size(pixel_size_um=75.0)
    calibration.set_distance(4.755)
    calibration.set_beam_position(838.5, 1833.4) # Approximate
        
    mask_dir = '/home/kyager/current/code/SciAnalysis/main/SciAnalysis/XSAnalysis/masks/'
    mask = Mask()
    #mask.load(mask_dir+'Eiger4M_all_gaps-mask.png')
    #mask.load(mask_dir+'CHX_Eiger4M-bad_flatfield_10percent-mask.png')
    mask.load(mask_dir+'CHX_Eiger4M-bad_flatfield_05percent-mask.png')
    mask.load(mask_dir+'CHX_pipe-2015Oct-mask.png')
    #mask.load(mask_dir+'CHX_bs_streaks-2015Oct-mask.png')
            
    
    from SciAnalysis.XSAnalysis import Protocols
    
    
    process = Protocols.ProcessorXS(calibration=calibration, mask=mask)
    protocols = [ Protocols.thumbnails(crop=1.3) ]
    #protocols = [ Protocols.circular_average(ylog=True) ]
    #protocols = [ Protocols.fit_calibration(ylog=True) ]
    
    process.run(infiles, protocols, output_dir=output_dir, force=True)
    
        
    
            
if False:
    

    
    root_dir = '/media/extend/CHX/'
    data_dir = 'LDRD_Meso_correlation/2015_10Oct-ref_patterns/21/'
    source_dir = root_dir+'data/'+data_dir
    output_dir = root_dir+'analysis/'+data_dir
    
    
    import glob
    #infiles = glob.glob(source_dir+'AgBH_series_24_master.h5')
    #infiles = glob.glob(source_dir+'*_master.h5')
    infiles = glob.glob(source_dir+'chip411_found_star3_x-0.264_y13.472_1052_master.h5')
    infiles = glob.glob(source_dir+'chip411_found_star_x-0.356_y13.464_937_master.h5')
    
    infiles.sort()
    
    
    
    # Experimental setup
    calibration = Calibration()
    calibration.set_energy(8.8984) # CHX
    calibration.set_image_size(2070, 2167) # Eiger 4M
    calibration.set_pixel_size(pixel_size_um=75.0)
    calibration.set_distance(4.755)
    calibration.set_beam_position(838.5, 1833.4) # Approximate
    
    
    mask_dir = '/home/kyager/current/code/SciAnalysis/main/SciAnalysis/XSAnalysis/masks/'
    mask = Mask()
    #mask.load(mask_dir+'Eiger4M_all_gaps-mask.png')
    #mask.load(mask_dir+'CHX_Eiger4M-bad_flatfield_10percent-mask.png')
    mask.load(mask_dir+'CHX_Eiger4M-bad_flatfield_05percent-mask.png')
    mask.load(mask_dir+'CHX_pipe-2015Oct-mask.png')
    #mask.load(mask_dir+'CHX_bs_streaks-2015Oct-mask.png')
    
    
    #calibration.set_beam_position(836.0, 1620.0) # Approximate
    #data = Data2DScattering(infiles[0], calibration=calibration, mask=mask)
    #data.plot(show=True)
    #line = data.circular_average(error=False)
    #line.smooth(2.0, bins=10)
    #line.plot(show=True, plot_range=[0, 0.27, 0, None], ylog=False, error=False, error_band=False)
    #line.save_data('out_q.data')
    
    #calibration.set_beam_position(837.2, 1621.5+12) # Starting guess for AgBH position
    
    data = Data2DScattering(infiles[0], calibration=calibration, mask=mask)
    data.threshold_pixels(4294967295-1) # Eiger inter-module gaps
    data.overlay_ring(0.1076, 0.0002) # AgBH
    data.overlay_ring(0.02, 0.0002)
    data.overlay_ring(0.01, 0.0002)
    data.overlay_ring(0.005, 0.0002)
    data.plot(show=True)
    exit()
    
    from SciAnalysis.XSAnalysis import Protocols
    
    
    process = Protocols.ProcessorXS(calibration=calibration, mask=mask)
    #protocols = [ Protocols.thumbnails(crop=1.3) ]
    #protocols = [ Protocols.circular_average(ylog=True) ]
    protocols = [ Protocols.fit_calibration(ylog=True) ]
    
    process.run(infiles, protocols, output_dir=output_dir, force=True)
    
    