# SciAnalysis

## Author : Kevin G. Yager
## Ported to git by : Julien Lhermitte
## Contributors: Ruipeng Li, Esther Tsai, Yugang Zhang

SciAnalysis is a set of Python scripts for batch processing of image data,
including x-ray scattering detector images. The code was written primarily by
Kevin Yager.

http://gisaxs.com/index.php/SciAnalysis

INSTALL: 

 * python setup.py develop
 
 OR, simply download the package and see examples/beamlines/NSLSII_11BM_CMS/UShell/SciAnalysis_jn.ipynb for tutorial.


---

Example protocals for X-ray scattering data:
 * Protocols.qr_image(blur=None, colorbar=True, save_results=['npz','hdf5'], transparent=False, label_filename=True) #plot_buffers = [0.1, 0.1, 0.1, 0.1], dpi=200
 * Protocols.linecut_angle(q0=2.30, dq=0.01, extra='_q2p30', show_region=True)
 * Protocols.linecut_qr(qz=0.025, dq=0.02, ylog=True, show_region=True, gridlines=True); #dq is half-width
 * Protocols.linecut_qz(name='linecut_qz_new', ylog=True, qr=0, dq=0.02, show_region=True, plot_range=[0.2, 0.4, None, None])
 * Protocols.linecut_qz_fit(qr=0.0, dq=0.01, show_region=True, label_filename=True, trim_range=[0.01, 0.4], fit_range=[0.093, 0.115], plot_range=[0.01, 0.4, 0, None], q0=[0.11]) 
 * Protocols.circular_average_q2I_fit(plot_range=[0.8, 1.3, 0, None], qn_power=0.0, trim_range=[0.1, 3.5], fit_range=[0.95, 1.4], num_curves=2, q0=[1.00, 1.2], sigma=0.02, show_curves=1, label_filename=True), 
 * Protocols.sector_average(angle=70, dangle=10, plot_range=[1.2, 3.7, 0, 1200], show_region=True) #pie-shaped
 * Protocols.roi(show_region=True, qx=1, dqx=0.02, qz=1, dqz=0.02, prepend='stats_')

run_args = { 'verbosity' : 3,
            #'save_results' : ['xml', 'plots', 'txt', 'hdf5'],
            'rcParams': {'axes.labelsize': 25,
                            'xtick.labelsize': 20,
                            'ytick.labelsize': 20,
                            'xtick.major.pad': 10,
                            'ytick.major.pad': 10,
                            },
            }


