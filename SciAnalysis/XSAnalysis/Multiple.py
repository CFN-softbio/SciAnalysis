#!/usr/bin/python
# -*- coding: utf-8 -*-
# vi: ts=4 sw=4
'''
:mod:`SciAnalysis.XSAnalysis.Multiple` - Data analysis protocols
================================================
.. module:: SciAnalysis.XSAnalysis.Multiple
   :synopsis: Convenient protocols for data analysis of multiple files.
.. moduleauthor:: Dr. Kevin G. Yager <kyager@bnl.gov>
                    Brookhaven National Laboratory
'''

################################################################################
#  Data analysis protocols that operate on file sequences (multiple files).
################################################################################
# Known Bugs:
#  N/A
################################################################################
# TODO:
#  Search for "TODO" below.
################################################################################

from .Data import *
from ..tools import *
from .Protocols import *

class _deprecated_ProtocolMultiple(Protocol):
 
    def preliminary(self, infiles, output_dir, **run_args):

        if run_args['processor'] is None:
            print('ERROR: {} needs a reference to the Processor (so that it can load files).'.format(self.name))
            return None
        
        if len(infiles)<1:
            print('ERROR: {} needs 1 or more files ({} files supplied)).'.format(self.name, len(infiles)))
            return None

        # Make sure all the requested infiles actually exist
        for infile in infiles:
            if infile is not os.path.isfile(infile):
                print('ERROR: {}; not all infiles exist ({} not found)'.format(self.name, infile))
                return None
        
        
        # Identify name_base (for saving output)
        if 'name_base' in run_args:
            name_base = run_args['name_base']
        elif 'name_re' in run_args:
            f = Filename(infiles[0])
            filepath, filename, filebase, ext = f.split()
            
            import re
            match_re = re.compile(run_args['name_re'])
            m = match_re.match(filename)
            if m:
                name_base = m.groups()[0]
                print(name_base)
        else:
            import os
            prefix = os.path.commonprefix(infiles)
            filepath, name_base = os.path.split(prefix)

        if 'append_protocol_name' in run_args:
            output_dir = output_dir + self.name
            make_dir(output_dir)
        outfile = self.get_outfile(name_base, output_dir, ext=run_args['file_extension'])
        
        if (not run_args['force']) and os.path.isfile(outfile):
            print(' Skipping (internal check) {} for {}'.format(self.name, data.name))
            return results

        if run_args['verbosity']>=3:
            print(' Running for {} files; saving to: {}'.format(len(infiles), outfile))
            
        return outfile
        
 
 
class sum_images(ProtocolMultiple):
    
    def __init__(self, name='sum_images', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.npy'
        self.run_args = {
                        'crop' : None,
                        'blur' : None,
                        'resize' : None,
                        'file_extension' : '-sum.npy',
                        'append_protocol_name' : True,
                        'force' : False,
                        'verbosity' : 3,
                        }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, datas, output_dir, basename, **run_args):
        
        results = {}
        
        outfile = self.get_outfile(basename, output_dir, ext=run_args['file_extension'])


        for data in datas:
            self.transform(data, **run_args)
        
        for data in datas[1:]:
            datas[0].data += data.data
        
        
        
        results['files_saved'] = [
            { 'filename': '{}'.format(outfile) ,
             'description' : 'sum of multiple images (npy format)' ,
             'type' : 'data' # 'data', 'plot'
            } ,
            ]
            
        np.save(outfile, data.data)
        
        
        return results
        
        
    def transform(self, data, **run_args):
        
        if run_args['crop'] is not None:
            data.crop(run_args['crop'])
        if run_args['blur'] is not None:
            data.blur(run_args['blur'])
        if run_args['resize'] is not None:
            data.resize(run_args['resize'])        
            
        return data
    
    
class stitch_images_position(ProtocolMultiple):
    '''Stitches images together into a single effective detector
    image, where the beam position may be different for the different
    images.'''
    
    def __init__(self, name='stitched', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.npy'
        self.run_args = {
                        'file_extension' : '.tiff',
                        'verbosity' : 3,
                        }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, datas, output_dir, basename, **run_args):
        
        results = {}
        
        outfile = self.get_outfile(basename, output_dir, ext=run_args['file_extension'])
        
        
        # Determine total image size
        h, w = datas[0].data.shape
        min_x, min_y = datas[0].calibration.x0, datas[0].calibration.y0
        max_x, max_y = min_x, min_y
        for data, position in zip(datas, run_args['positions']):
            data.calibration.use_beam_position(position)
            x0, y0 = data.calibration.x0, data.calibration.y0
            if run_args['verbosity']>=5:
                print("    Position '{}' (x0, y0) = ({:,.1f}, {:,.1f}) for {}".format(position, x0, y0, data.name))
            if x0<min_x: min_x = x0
            if x0>max_x: max_x = x0
            if y0<min_y: min_y = y0
            if y0>max_y: max_y = y0

        span_x, span_y = int(np.ceil(max_x-min_x)), int(np.ceil(max_y-min_y))
        if run_args['verbosity']>=5:
            print("    Each detector image is: ({:d}, {:d})".format(w, h))
            print("    Maximum spread in beam positions: ({:d}, {:d})".format(span_x, span_y))
            print("    Full image is: ({:d}, {:d})".format(w+span_x, h+span_y))
            
            
        # Create an image array big enough to hold all detector images
        Intensity_map = np.zeros( (h+span_y, w+span_x) )
        count_map = np.zeros( (h+span_y, w+span_x) )
        
        # Copy each image into the full array
        for data, position in zip(datas, run_args['positions']):
            data.calibration.use_beam_position(position)
            x0, y0 = data.calibration.x0, data.calibration.y0
            
            xi = max_x - x0
            yi = max_y - y0
            xf = xi+w
            yf = yi+h
            Intensity_map[yi:yf, xi:xf] += data.data[:,:]
            
            if data.mask is None:
                count_map[yi:yf, xi:xf] += np.ones(data.data.shape)
            else:
                count_map[yi:yf, xi:xf] += data.mask.data
            

        Intensity_map = np.nan_to_num( Intensity_map/count_map )
        Intensity_map = Intensity_map.astype(np.uint32)


        results['files_saved'] = [
            { 'filename': '{}'.format(outfile) ,
             'description' : 'merging of multiple images into common detector image (tiff format)' ,
             'type' : 'data' # 'data', 'plot'
            } ,
            ]
        img = PIL.Image.fromarray(Intensity_map)
        img.save(outfile)
        

        
        return results
            
    
class merge_images_position(ProtocolMultiple):
    '''Merges images into reciprocal-space, where the different images
    have different detector positions.'''
    
    def __init__(self, name='merge_images', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.npy'
        self.run_args = {
                        'file_extension' : '-merged.npy',
                        'verbosity' : 3,
                        'save_axes' : True,
                        'save_mask' : True,
                        'save_maps' : True,
                        }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, datas, output_dir, basename, **run_args):
        
        results = {}
        
        outfile = self.get_outfile(basename, output_dir, ext=run_args['file_extension'])

        # Prepare a large-area q matrix
        q_range = run_args['q_range']
        dq = run_args['dq']
        qxs = np.arange(q_range[0], q_range[1], dq)
        qzs = np.arange(q_range[2], q_range[3], dq)
        QXs, QZs = np.meshgrid(qxs, qzs)
        
        Intensity_map = np.zeros( (len(qzs), len(qxs)) )
        count_map = np.zeros( (len(qzs), len(qxs)) )
        
        if run_args['verbosity']>=5:
            print('      Expecting array size num_qx = {}, num_qz = {}'.format(len(qxs), len(qzs)))
            print('      Coordinate matrices sized {}'.format(QXs.shape))
            #print('      Coordinate matrices sized {}'.format(QZs.shape))
            print('      Data matrices sized {}'.format(Intensity_map.shape))
            #print('      Data matrices sized {}'.format(count_map.shape))


        
        import re
        name_re = re.compile('.+_(pos\d)_.+')
        for data in datas:
            m = name_re.match(data.name)
            if m:
                data.calibration.use_beam_position( m.groups()[0] )
            else:
                data.calibration.use_beam_position('default')
            data.calibration.clear_maps()
            
            remesh_data, num_per_pixel = data.remesh_q_bin_explicit(qx_min=q_range[0], qx_max=q_range[1], num_qx=len(qxs), qz_min=q_range[2], qz_max=q_range[3], num_qz=len(qzs), **run_args)

            if run_args['verbosity']>=5:
                print('      remesh_data matrix sized {}'.format(remesh_data.shape))
            
            Intensity_map += remesh_data
            count_map += num_per_pixel
            
            

        Intensity_map = np.nan_to_num( Intensity_map/count_map )
        
        results['files_saved'] = [
            { 'filename': '{}'.format(outfile) ,
             'description' : 'merging of multiple images into common q-space (npy format)' ,
             'type' : 'data' # 'data', 'plot'
            } ,
            ]
        np.save(outfile, Intensity_map)
        

        if run_args['save_axes']:
            outfile = self.get_outfile(basename, output_dir, ext='-axes.npz')
            
            results['files_saved'].append( 
                { 'filename': '{}'.format(outfile) ,
                'description' : 'qx and qz axes of output data (npz format)' ,
                'type' : 'metadata' , # 'data', 'plot', 'metadata'
                'comment' : "reload using: npzfile = np.load('{}{}')".format(basename, '-axes.npz')
                } , )
            np.savez_compressed(outfile, qxs=qxs, qzs=qzs)


        if run_args['save_maps']:
            QYs = datas[0].calibration.compute_qy(QXs, QZs)
            outfile = self.get_outfile(basename, output_dir, ext='-maps.npz')
            np.savez_compressed(outfile, QX=QXs, QY=QYs, QZ=QZs)
            

        if run_args['save_mask']:
            
            mask = np.where(count_map==0, np.zeros_like(count_map), np.ones_like(count_map))
            
            outfile = self.get_outfile(basename, output_dir, ext='-mask.png')
            #np.save(outfile, mask)
            
            try:
                import scipy.misc
                scipy.misc.imsave(outfile, mask)
            except:
                import imageio
                imageio.imwrite(outfile, mask)

            results['files_saved'].append( 
                { 'filename': '{}'.format(outfile) ,
                'description' : 'mask file for the data' ,
                'type' : 'metadata' , # 'data', 'plot', 'metadata'
                } , )
            
        
        
        return results
        

class merge_images_swaxs(merge_images_position):
    '''Merges images into reciprocal-space, where the different images
    have different detector positions.
    This protocol is designed for merging data from SAXS and WAXS.'''

    def __init__(self, name='merge_images', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.npy'
        self.run_args = {
                        'file_extension' : '-merged.npy',
                        'verbosity' : 3,
                        'save_axes' : True,
                        'save_mask' : True,
                        'save_maps' : True,
                        'method': 'linear', # 'linear', 'nearest',  or 'cubic'
                        }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, datas, output_dir, basename, **run_args):
        
        results = {}
        
        outfile = self.get_outfile(basename, output_dir, ext=run_args['file_extension'])

        # Prepare a large-area q matrix
        q_range = run_args['q_range']
        dq = run_args['dq']
        qxs = np.arange(q_range[0], q_range[1]+dq, dq)
        qzs = np.arange(q_range[2], q_range[3]+dq, dq)
        QXs, QZs = np.meshgrid(qxs, qzs)
        
        Intensity_map = np.zeros( (len(qzs), len(qxs)) )
        count_map = np.zeros( (len(qzs), len(qxs)) )
        
        if run_args['verbosity']>=5:
            print('      Expecting array size num_qx = {}, num_qz = {}'.format(len(qxs), len(qzs)))
            print('      Coordinate matrices sized {}'.format(QXs.shape))
            #print('      Coordinate matrices sized {}'.format(QZs.shape))
            print('      Data matrices sized {}'.format(Intensity_map.shape))
            #print('      Data matrices sized {}'.format(count_map.shape))
        
        for ii, data in enumerate(datas):
            if run_args['verbosity']>=5:
                print("[{}] data.name = {}".format(ii, data.name))
                #print("    q_per_pixel = {}".format(data.calibration.get_q_per_pixel()))
            
            # Using run_args['dq']
            remesh_data, num_per_pixel = data.remesh_q_interpolate_explicit(qx_min=q_range[0], qx_max=q_range[1], qz_min=q_range[2], qz_max=q_range[3],  **run_args)

            if run_args['verbosity']>=5:
                print('      remesh_data matrix sized {}'.format(remesh_data.shape))
                print("      remesh_data MAX {}".format(np.max(remesh_data)))
                print("      remesh_data SUM {}".format(np.sum(remesh_data)))
                print("      num_per_pixel SUM {}".format(np.sum(num_per_pixel)))
            
            Intensity_map += remesh_data
            count_map += num_per_pixel
                    
        #Intensity_map = np.nan_to_num( Intensity_map/count_map )
        
        results['files_saved'] = [
            { 'filename': '{}'.format(outfile) ,
             'description' : 'merging of multiple images into common q-space (npy format)' ,
             'type' : 'data' # 'data', 'plot'
            } ,
            ]
        np.save(outfile, Intensity_map)
        
        
        # Save to PNG (Temporary)
        print("-------")
        from PIL import Image
        outfile2 = self.get_outfile(basename, output_dir, ext='.png')
        img_out = np.nan_to_num((Intensity_map))
        print(np.min(img_out))
        print(np.max(img_out))
        rescaled = (255.0 / img_out.max() * (img_out - img_out.min())).astype(np.uint8)
        print(np.sum(rescaled))
        Image.fromarray(rescaled).save(outfile2)
        

        if run_args['save_axes']:
            outfile = self.get_outfile(basename, output_dir, ext='-axes.npz')
            
            results['files_saved'].append( 
                { 'filename': '{}'.format(outfile) ,
                'description' : 'qx and qz axes of output data (npz format)' ,
                'type' : 'metadata' , # 'data', 'plot', 'metadata'
                'comment' : "reload using: npzfile = np.load('{}{}')".format(basename, '-axes.npz')
                } , )
            np.savez_compressed(outfile, qxs=qxs, qzs=qzs)


        if run_args['save_maps']:
            QYs = datas[0].calibration.compute_qy(QXs, QZs)
            outfile = self.get_outfile(basename, output_dir, ext='-maps.npz')
            np.savez_compressed(outfile, QX=QXs, QY=QYs, QZ=QZs)
            

        if run_args['save_mask']:
            
            mask = np.where(count_map==0, np.zeros_like(count_map), np.ones_like(count_map))
            
            outfile = self.get_outfile(basename, output_dir, ext='-mask.png')
            #np.save(outfile, mask)
            
            try:
                import scipy.misc
                scipy.misc.imsave(outfile, mask)
            except:
                import imageio
                imageio.imwrite(outfile, mask)

            results['files_saved'].append( 
                { 'filename': '{}'.format(outfile) ,
                'description' : 'mask file for the data' ,
                'type' : 'metadata' , # 'data', 'plot', 'metadata'
                } , )
            
        return results        

            
    
class merge_images_gonio_phi(ProtocolMultiple):
    '''Merges images into reciprocal-space, where the different images
    have different detector angular positions (for detectors on a
    goniometer circle centered around the sample)..'''
    
    def __init__(self, name='merge_images', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.npy'
        self.run_args = {
                        'file_extension' : '-merged.npy',
                        'verbosity' : 3,
                        'phi_scaling' : -1.0,
                        'det_theta_g' : 0.0,
                        'save_axes' : True,
                        'save_mask' : True,
                        'save_maps' : True,
                        }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, datas, output_dir, basename, **run_args):
        
        results = {}
        
        outfile = self.get_outfile(basename, output_dir, ext=run_args['file_extension'])

        # Prepare a large-area q matrix
        q_range = run_args['q_range']
        dq = run_args['dq']
        qxs = np.arange(q_range[0], q_range[1], dq)
        qzs = np.arange(q_range[2], q_range[3], dq)
        QXs, QZs = np.meshgrid(qxs, qzs)
        
        Intensity_map = np.zeros( (len(qzs), len(qxs)) )
        count_map = np.zeros( (len(qzs), len(qxs)) )
        
        if run_args['verbosity']>=5:
            print('      Expecting array size num_qx = {}, num_qz = {}'.format(len(qxs), len(qzs)))
            print('      Coordinate matrices sized {}'.format(QXs.shape))
            #print('      Coordinate matrices sized {}'.format(QZs.shape))
            print('      Data matrices sized {}'.format(Intensity_map.shape))
            #print('      Data matrices sized {}'.format(count_map.shape))



        
        for data in datas:
            phi = self.get_phi(data, **run_args)*run_args['phi_scaling']
            data.calibration.set_angles(det_phi_g=phi, det_theta_g=run_args['det_theta_g'])
            
            remesh_data, num_per_pixel = data.remesh_q_bin_explicit(qx_min=q_range[0], qx_max=q_range[1], num_qx=len(qxs), qz_min=q_range[2], qz_max=q_range[3], num_qz=len(qzs), **run_args)

            if run_args['verbosity']>=5:
                print('      remesh_data matrix sized {}'.format(remesh_data.shape))
            
            Intensity_map += remesh_data
            count_map += num_per_pixel
            
            

        Intensity_map = np.nan_to_num( Intensity_map/count_map )
        
        results['files_saved'] = [
            { 'filename': '{}'.format(outfile) ,
             'description' : 'merging of multiple images into common q-space (npy format)' ,
             'type' : 'data' # 'data', 'plot'
            } ,
            ]
        np.save(outfile, Intensity_map)
        
        if False:
            # Old method: save using pkl format
            import pickle
            outfile = self.get_outfile(basename, output_dir, ext='-axes.pkl')
            
            results['files_saved'].append( 
                { 'filename': '{}'.format(outfile) ,
                'description' : 'qx and qz axes of output data (Python pickle format)' ,
                'type' : 'metadata' , # 'data', 'plot', 'metadata'
                'comment' : "reload using: qxs, qzs = pickle.load(open('{}{}','rb'))".format(basename, '-axes.pkl')
                } , )
            with open(outfile, 'wb') as fout:
                pickle.dump([qxs, qzs], fout)
                
            # Reload using:
            #qxs, qzs = pickle.load(open('infile.pkl','rb'))        

        if run_args['save_axes']:
            outfile = self.get_outfile(basename, output_dir, ext='-axes.npz')
            
            results['files_saved'].append( 
                { 'filename': '{}'.format(outfile) ,
                'description' : 'qx and qz axes of output data (npz format)' ,
                'type' : 'metadata' , # 'data', 'plot', 'metadata'
                'comment' : "reload using: npzfile = np.load('{}{}')".format(basename, '-axes.npz')
                } , )
            np.savez_compressed(outfile, qxs=qxs, qzs=qzs)

        if run_args['save_maps']:
            QYs = datas[0].calibration.compute_qy(QXs, QZs)
            outfile = self.get_outfile(basename, output_dir, ext='-maps.npz')
            np.savez_compressed(outfile, QX=QXs, QY=QYs, QZ=QZs)
            

        if run_args['save_mask']:
            
            mask = np.where(count_map==0, np.zeros_like(count_map), np.ones_like(count_map))
            
            outfile = self.get_outfile(basename, output_dir, ext='-mask.png')
            #np.save(outfile, mask)
            import scipy.misc
            scipy.misc.imsave(outfile, mask)            

            results['files_saved'].append( 
                { 'filename': '{}'.format(outfile) ,
                'description' : 'mask file for the data (PNG file)' ,
                'type' : 'metadata' , # 'data', 'plot', 'metadata'
                } , )
            
        
        
        return results
        
        
    def get_phi(self, data, **run_args):
        
        phi_re = re.compile(run_args['phi_re'])
        m = phi_re.match(data.name)
        if m:
            phi = float(m.groups()[0])
            if run_args['verbosity']>=5:
                print('    Extracted phi = {} from {}'.format(phi, data.name))
            return phi
        else:
            if run_args['verbosity']>=1:
                print('ERROR: Could not extract phi for {}'.format(data.name))
            return None
        
        
        
    
