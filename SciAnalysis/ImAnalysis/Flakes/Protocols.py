#!/usr/bin/python
# -*- coding: utf-8 -*-
# vi: ts=4 sw=4

from ..Protocols import *

import sys
import pickle

from PIL import Image
import matplotlib.patches as patches

from scipy import ndimage
from scipy.stats import entropy
from scipy.spatial.distance import cdist
import skimage
import skimage.exposure
import skimage.measure as measure
import cv2







class thumbnails_contrast(thumbnails):
    
    def __init__(self, name='thumbnails', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.jpg'
        self.run_args = {
                        'crop' : 0.5,
                        'blur' : 1.0,
                        'resize' : 0.5,
                        }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        outfile = self.get_outfile(data.name, output_dir)
        data.plot_image(save=outfile, size=10*run_args['resize'], **run_args)
        
        return results
        

    
    
    
class find_flakes(thumbnails):
    
    def __init__(self, name='find_flakes', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
                        'background' : None ,
                        'image_threshold' : 10 , # Threshold to define the foreground in the residual map
                        'size_threshold' : 0 , # Size (in pixels) to consider a contiguous region to be a foregrond object
                        'verbosity' : 3 ,
                        'show' : False ,
                        'dpi' : 'auto' ,
                        'resize' : 1.0 ,
                        'bbox_expanded' : 0.1 ,
                        'overlays' : 4 , # Larger number adds more annotations
                        }
        self.run_args.update(kwargs)
        
        
        # WARNING: Potential namespace collision with Python 'utils' package
        try:
            # Boyu Wang (Stony Brook University) code:
            #  https://github.com/Boyu-Wang/material_segmentation
            from MaterialSegmentation import utils as MatSegUtils
            #from MaterialSegmentation import MatSegUtils
        except:
            code_PATH='/home/qpress/current/code/MaterialSegmentation/main/'
            code_PATH in sys.path or sys.path.insert(0, code_PATH)
            import utils as MatSegUtils
            #import MatSegUtils
            
        self.MatSegUtils = MatSegUtils
        

    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        scale = np.average((data.x_scale, data.y_scale)) # um/pixel
        scale2 = scale*scale # um^2/pixel^2
        
        # Load image
        ########################################
        if run_args['verbosity']>=5:
            print('  Loading image')
        
        image = Image.open(data.infile).convert('RGB')
        im_rgb = np.array(image).astype('float')
        im_gray = np.array(image.convert('L', (0.2989, 0.5870, 0.1140, 0))).astype('float') # vals in range [0..255]
        # im_hsv = np.array(image.convert('HSV')).astype('float') # vals in range [0..255]
        # Alternate HSV conversion (consistent with Matlab)
        im_hsv = skimage.color.rgb2hsv(im_rgb)
        im_hsv[:,:,2] = im_hsv[:,:,2]/255.0 # vals in range [0..1]

        h, w, c = im_rgb.shape
        if run_args['verbosity']>=5:
            print('    {}×{} = {:.2f} Mpix; {:.1f}μm × {:.1f}μm = {:.3f} mm^2'.format(w, h, w*h/1e6, w*scale, h*scale, w*h*scale2/1e6))

        
        if run_args['background']:
            background_image = Image.open(run_args['background']).convert('RGB')
            bk_rgb = np.array(background_image).astype('float')
            bk_gray = np.array( background_image.convert('L', (0.2989, 0.5870, 0.1140, 0)) ).astype('float')
            # bk_hsv = np.array(background_image.convert('HSV')).astype('float')
            # Alternate HSV conversion (consistent with Matlab)
            bk_hsv = skimage.color.rgb2hsv(bk_rgb)
            bk_hsv[:,:,2] = bk_hsv[:,:,2]/255.0
            
            im_maps = im_gray, im_hsv, im_rgb, bk_gray, bk_hsv
            
        else:
            im_maps = im_gray, im_hsv, im_rgb, None, None
        
        
        
        # Flake segmentation
        ########################################
            
        if run_args['verbosity']>=4:
            print('  Segmenting image to find flakes')
            
        
        #res_map, image_labelmap, flake_centroids, flake_sizes, num_flakes = self.MatSegUtils.perform_robustfit_multichannel(im_hsv, im_gray, run_args['image_threshold'], run_args['size_threshold'])

        res_map, image_labelmap, flake_centroids, flake_sizes, num_flakes = self.perform_background_threshold(im_hsv, im_gray, background_image, **run_args)
        
        
        
        if run_args['verbosity']>=4:
            print('    {} flakes'.format(num_flakes))
        
        if run_args['background']:
            if run_args['verbosity']>=4:
                print('  Segmenting background')
            _, _, bk_flake_centroids, _, n = self.MatSegUtils.perform_robustfit(bk_gray, 10, 0)
            
            if n>0:
                # Remove regions in background
                dis = cdist(flake_centroids, bk_flake_centroids)
                dis_min = dis.min(1)
                # to_remove = np.nonzero(dis_min < 5)[0]
                to_keep = np.nonzero(dis_min >= 5)[0]

                flake_sizes = flake_sizes[to_keep]
                flake_centroids = flake_centroids[to_keep]
                num_flakes = to_keep.shape[0]
                if run_args['verbosity']>=4:
                    print('    {} background objects; {} flakes remain'.format(n, num_flakes))
                
                # Reprocess label map
                new_image_labels = np.zeros(image_labelmap.shape)
                cnt = 0
                for idx in to_keep:
                    cnt += 1
                    new_image_labels[image_labelmap==idx+1] = cnt
                image_labelmap = new_image_labels            
            
        results['num_flakes'] = num_flakes


        if run_args['verbosity']>=10:
            # Diagnostics
            import scipy.misc
            
            outfile = self.get_outfile(data.name, output_dir, ext='-residuals.png')
            #scipy.misc.toimage(res_map, cmin=0, cmax=255).save(outfile)
            Image.fromarray(res_map.astype(np.uint8)).save(outfile)

            arr = image_labelmap>0
            outfile = self.get_outfile(data.name, output_dir, ext='-binary.png')
            #scipy.misc.toimage(arr, cmin=0, cmax=1).save(outfile)
            Image.fromarray( (np.clip(arr, 0, 1)*255).astype(np.uint8) ).save(outfile)
            
            arr = np.dstack([arr]*3)
            im = np.where(arr, im_rgb, np.zeros_like(im_rgb))
            outfile = self.get_outfile(data.name, output_dir, ext='-flakes.png')
            #scipy.misc.toimage(im, cmin=0, cmax=255).save(outfile)
            Image.fromarray(im.astype(np.uint8)).save(outfile)
            
            outfile = self.get_outfile(data.name, output_dir, ext='-colored.png')
            self.colored_objects(image_labelmap, outfile, **run_args)

        
        # Prepre output image
        ########################################
        aspect = w/h
        if not isinstance(run_args['dpi'], (int, float)):
            # Try to retain the 'real' size of the image
            run_args['dpi'] = w/10
        self.fig = plt.figure( figsize=(10,10/aspect), facecolor='white' )
        self.ax = self.fig.add_axes( [0, 0, 1, 1] )
        
        if run_args['image_contrast'] is not None:
            in_range = ( run_args['image_contrast'][0]*255, run_args['image_contrast'][1]*255 )
            data_rgb = skimage.exposure.rescale_intensity(im_rgb, in_range=in_range, out_range='dtype')
        self.ax.imshow(data_rgb)
        


        # Analyze each flake
        ########################################
        if run_args['verbosity']>=4:
            print('  Analyzing {} flakes'.format(num_flakes))
        
        flakes = []
        kernel = np.ones((5,5), np.uint8)
        for i in range(num_flakes):
            if run_args['verbosity']>=5:
                print('    Flake {}/{} ({:.1f}%)'.format(i+1, num_flakes, 100.*i/num_flakes))
                
            flake_i = self.flake_analysis(image_labelmap, im_maps, index=i+1, kernel=kernel, **run_args)

            flake_i['name'] = data.name
            flake_i['infile'] = data.infile
            flake_i['size'] = flake_sizes[i] # Surface area in pixels
            flake_i['scale'] = scale
            flake_i['scale2'] = scale2

            
            flake_i = self.flake_plot(flake_i, scale, **run_args)

            flakes.append(flake_i)

        

        # Finalize
        ########################################
        
        # Save image
        #xi, xf, yi, yf = self.ax.axis()
        self.ax.axis( [0, w, h, 0] )
        
        if run_args['show']:
            plt.show()
        outfile = self.get_outfile(data.name, output_dir)
        plt.savefig(outfile, dpi=run_args['dpi']*run_args['resize'])
        plt.close(self.fig.number)
        
        
        # Save flake information
        to_save = {}
        to_save['res_map'] = res_map
        to_save['image_labelmap'] = image_labelmap
        to_save['flakes'] = flakes
        outfile = self.get_outfile(data.name, output_dir, ext='.pkl')
        with open(outfile, 'wb') as fout:
            pickle.dump(to_save, fout)
        
        return results



    def perform_background_threshold(self, im_hsv, im_gray, bk_rgb, **run_args):
        
        #res_map, image_labelmap, flake_centroids, flake_sizes, num_flakes = self.MatSegUtils.perform_robustfit_multichannel(im_hsv, im_gray, run_args['image_threshold'], run_args['size_threshold'])
        
        # Use gray, hue, saturation
        im_ghs = np.concatenate([np.expand_dims(im_gray/255,2), im_hsv[:,:,:2]], axis=2) # vals in range [0..1]
        
        
        bk_gray = np.array( bk_rgb.convert('L', (0.2989, 0.5870, 0.1140, 0)) ).astype('float') # vals in range [0..255]
        bk_rgb = np.array(bk_rgb).astype('float')
        bk_hsv = skimage.color.rgb2hsv(bk_rgb)
        bk_hsv[:,:,2] = bk_hsv[:,:,2]/255.0 # vals in range [0..1]
        bk_ghs = np.concatenate([np.expand_dims(bk_gray/255,2), bk_hsv[:,:,:2]], axis=2) # vals in range [0..1]


        h, w, c = im_ghs.shape
        res_map = np.zeros([h,w])
        for ic in range(c):
            im_c = im_ghs[:,:,ic]
            #val_stats(im_c, name=ic)
            bk_c = bk_ghs[:,:,ic]

            if ic==1:
                # Hue channel is special due to the cyclic nature of its definition (where 0 and 1 are same number)
                r = np.min( [np.abs(im_c-bk_c) , 1-np.abs(im_c-bk_c)], axis=0 )
                res_map +=  r
            else:
                res_map +=  np.abs(im_c-bk_c)
            
        thresh = c*run_args['image_threshold']/255
        outlier_map = res_map>thresh
        outlier_map = outlier_map.astype(np.uint8)
        
        # connected component detection
        nCC, image_labelmap, _, flake_centroids = cv2.connectedComponentsWithStats(outlier_map)
        # flake_centroid: [m, 2] array, indicates row, column of the centroid
        flake_centroids = np.flip(flake_centroids, 1)
        flake_centroids = flake_centroids[1:].astype('int')

        _, flake_sizes = np.unique(image_labelmap, return_counts=True)
        # remove the background size
        flake_sizes = flake_sizes[1:]
        if run_args['size_threshold'] > 0:
            # remove small connect component
            large_flakes = flake_sizes > run_args['size_threshold']
            large_flake_idxs = np.nonzero(large_flakes)[0]
            new_image_labels = np.zeros([h, w])
            cnt = 0
            for idx in large_flake_idxs:
                cnt += 1
                new_image_labels[image_labelmap==idx+1] = cnt
            image_labelmap = new_image_labels
            num_flakes = large_flakes.sum()
            flake_centroids = flake_centroids[large_flake_idxs]
            flake_sizes = flake_sizes[large_flake_idxs]
        else:
            num_flakes = nCC - 1        
        
        return res_map, image_labelmap, flake_centroids, flake_sizes, num_flakes
        

        
    def flake_plot(self, flake_i, scale, **run_args):
        '''Adds information about the flake to the plot.
        Also computes some additional things about the flake.'''
        
        scale2 = scale*scale
        
        flake_i['size_um'] = flake_i['size_pixels']*scale2
        flake_i['radius_pixels'] = np.sqrt(flake_i['size']/np.pi)
        flake_i['radius_um'] = flake_i['radius_pixels']*scale
        flake_i['perimeter_um'] = flake_i['perimeter_pixels']*scale 
        flake_i['P_over_A'] = flake_i['perimeter_um']/flake_i['size_um']
        flake_i['Pr_over_A'] = flake_i['perimeter_um']*flake_i['radius_um']/flake_i['size_um']
        flake_i['P_over_C'] = flake_i['perimeter_um']/(2*np.pi*flake_i['radius_um'])
        
        flake_i['contour_perimeter_um'] = flake_i['contour_perimeter_pixels']*scale
        flake_i['contour_size_um'] = flake_i['contour_size_pixels']*scale2
        flake_i['convex_perimeter_um'] = flake_i['convex_perimeter_pixels']*scale
        flake_i['P_over_chP'] = flake_i['perimeter_um']/flake_i['convex_perimeter_um']
        flake_i['convex_size_um'] = flake_i['convex_size_pixels']*scale2
        flake_i['convex_fill'] = flake_i['size_pixels']/flake_i['convex_size_pixels']
        
        
        if run_args['overlays']>=1:
            c = flake_i['contour']
            self.ax.plot(c[:,0], c[:,1], '-', linewidth=0.5, color='r', dashes=[4,2], alpha=0.3)
        if run_args['overlays']>=7:
            c = flake_i['convex_hull']
            self.ax.plot(c[:,1], c[:,0], '-', linewidth=0.5, color='g')
        
        if run_args['overlays']>=3:
            # Circle overlay denoting size
            y, x = flake_i['center_of_mass']
            size = flake_i['radius_pixels']
            circ = patches.Circle(xy=(x,y), radius=size, linewidth=0.8, edgecolor='r', facecolor='none', alpha=0.3)
            self.ax.add_patch(circ)
            s = r'${:.1f} \, \mathrm{{\mu m}}$'.format(flake_i['radius_um'])
            self.ax.text(x, y+size, s, size=6, color='r', horizontalalignment='center', verticalalignment='top', alpha=0.7)
            
        if run_args['overlays']>=4:
            # Cross hair
            rect = patches.Rectangle( (x-size/2, y), size, 0, linewidth=0.8, edgecolor='r', facecolor='none', alpha=0.3) # Horizontal bar
            self.ax.add_patch(rect)
            rect = patches.Rectangle( (x, y-size/2), 0, size, linewidth=0.8, edgecolor='r', facecolor='none', alpha=0.3) # Vertical bar
            self.ax.add_patch(rect)
        
        
        # Bounding box
        y1, y2, x1, x2 = flake_i['bbox']
        bbox_size = abs(y2-y1)*abs(x2-x1) # pixels^2
        flake_i['bbox_size'] = bbox_size
        flake_i['bbox_fill'] = flake_i['size']/bbox_size
        if run_args['overlays']>=2:
            rect = patches.Rectangle( (x1,y1), x2-x1, y2-y1, linewidth=1.2, edgecolor='orange', facecolor='none', alpha=0.5)
            self.ax.add_patch(rect)

        if run_args['overlays']>=3:
            s = '{:.3f}'.format(flake_i['flake_contrast'])
            self.ax.text(x2, y2, s, size=6, color='orange', horizontalalignment='left', verticalalignment='top', alpha=0.9)

        if run_args['overlays']>=5:
            # Location annotations
            s = r'$({:d}, {:d})$'.format(int(x), int(y))
            self.ax.text(x2, y1, s, size=6, color='orange', horizontalalignment='right', verticalalignment='bottom', alpha=0.9)
            s = '{:d}'.format(i+1)
            self.ax.text(x2, y1, s, size=6, color='orange', horizontalalignment='right', verticalalignment='top', alpha=0.9)

        if run_args['overlays']>=7:
            # Extra information about flake
            s1 = r'$f_{{\mathrm{{bbox}}}} = {:.1f}\%$'.format(flake_i['bbox_fill']*100.)
            s2 = r'$\frac{{P}}{{C}} = {:.2f}$'.format(flake_i['Pr_over_A'], flake_i['P_over_C'])
            s3 = r'$f_{{\mathrm{{convex}}}} = {:.1f}\%$'.format(flake_i['convex_fill']*100.)
            s = '\n'.join((s1,s2,s3))
            self.ax.text(x1, y2, s, size=6, color='orange', horizontalalignment='left', verticalalignment='top', alpha=0.9)

        if run_args['overlays']>=6:
            # Expanded bounding box
            y1, y2, x1, x2 = flake_i['bbox_expanded']
            rect = patches.Rectangle( (x1,y1), x2-x1, y2-y1, linewidth=1.0, linestyle='dashed', edgecolor='orange', facecolor='none', alpha=0.25)
            self.ax.add_patch(rect)

        if run_args['verbosity']>=5:
            print('      (x, y) = ({:.1f}, {:.1f})'.format(x, y))
            print('      {} pixels ({:.1f}% of {} pix bbox, {:.1f}% of {:d} pix convex hull)'.format(flake_i['size'], 100.*flake_i['bbox_fill'], bbox_size, 100.*flake_i['convex_fill'], int(flake_i['convex_size_pixels'])))
            print('      A = {:.1f} μm^2; r = {:.1f} μm'.format(flake_i['size_um'], flake_i['radius_um']))
            print('      P = {:.1f} μm; P/A = {:.1f} μm^-1; Pr/A = {:.2f}; P/C = {:.2f}'.format(flake_i['perimeter_um'], flake_i['P_over_A'], flake_i['Pr_over_A'], flake_i['P_over_C']))
            print('      P/P_contour = {:.2f}; P/P_ch = {:.2f}; P_contour/P_ch = {:.2f}'.format(flake_i['perimeter_um']/flake_i['contour_perimeter_um'], flake_i['P_over_chP'], flake_i['contour_perimeter_um']/flake_i['convex_perimeter_um']))
                    
        
        return flake_i
    
    

    def flake_analysis(self, image_labelmap, im_maps, index, kernel, **run_args):
        
        h, w = image_labelmap.shape
        im_gray, im_hsv, im_rgb, bk_gray, bk_hsv = im_maps

        
        flake_i = {}
        
        flake_region = (image_labelmap==index).astype(np.uint8)
        flake_i['index'] = index
        flake_i['size_pixels'] = np.sum(flake_region)
        flake_i['center_of_mass'] = ndimage.measurements.center_of_mass(flake_region)
        flake_i['perimeter_pixels'] = measure.perimeter(flake_region)
        
        
        # Bounding box [y1, y2, x1, x2]
        f_mask_r, f_mask_c = np.nonzero(flake_region)
        f_mask_r_min = min(f_mask_r)
        f_mask_r_max = max(f_mask_r)
        f_mask_height = f_mask_r_max - f_mask_r_min
        f_mask_c_min = min(f_mask_c)
        f_mask_c_max = max(f_mask_c)
        f_mask_width = f_mask_c_max - f_mask_c_min
        flake_i['bbox'] = [f_mask_r_min, f_mask_r_max+1, f_mask_c_min, f_mask_c_max+1]
        
        # Expanded bounding box
        expand = run_args['bbox_expanded']
        flake_large_bbox = [max(0, f_mask_r_min - int(expand*f_mask_height)), min(h, f_mask_r_max + int(expand*f_mask_height)),
                                    max(0, f_mask_c_min - int(expand*f_mask_width)), min(w, f_mask_c_max + int(expand*f_mask_width))]
        f_mask_large = max(f_mask_height, f_mask_width)
        flake_large_bbox_border = [max(0, f_mask_r_min - int(expand*f_mask_large)), min(h, f_mask_r_max + int(expand*f_mask_large)),
                                    max(0, f_mask_c_min - int(expand*f_mask_large)), min(w, f_mask_c_max + int(expand*f_mask_large))]
        yc, xc = (f_mask_r_min+f_mask_r_max)*0.5, (f_mask_c_min+f_mask_c_max)*0.5
        flake_large_bbox_square = [max(0, yc - int((0.5+expand)*f_mask_large)), min(h, yc + int((0.5+expand)*f_mask_large)),
                                    max(0, xc - int((0.5+expand)*f_mask_large)), min(w, xc + int((0.5+expand)*f_mask_large))]
        
        # Different definitions of the box:
        #flake_i['bbox_expanded'] = flake_large_bbox # bbox expanded by a certain % in both directions
        #flake_i['bbox_expanded'] = flake_large_bbox_border # Uniform sized border
        flake_i['bbox_expanded'] = flake_large_bbox_square # Square box
        
        
        # Flake contour (call signature depends on cv2 version)
        #_, flake_contour, _ = cv2.findContours(flake_region, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        flake_contour, _ = cv2.findContours(flake_region, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        flake_i['contour'] = np.squeeze(flake_contour[0], 1)
        flake_i['contour_perimeter_pixels'] = self.total_length(flake_i['contour'])
        flake_i['contour_size_pixels'] = self.polygon_area(flake_i['contour'])

        # Compute convex hull of the contour
        flake_contour = np.flip(flake_i['contour'], 1)
        flake_i['convex_hull'] = np.squeeze(cv2.convexHull(flake_contour), 1)
        flake_i['convex_perimeter_pixels'] = self.total_length(flake_i['convex_hull'])
        flake_i['convex_size_pixels'] = self.polygon_area(flake_i['convex_hull'])
        
        # Calculate additional shape features
        flake_shape_len_area_ratio = flake_contour.shape[0] / (flake_i['size_pixels']*1.)
        contours_center_dis = cdist(np.expand_dims(flake_i['center_of_mass'],0), flake_contour)
        flake_shape_contour_hist = np.histogram(contours_center_dis, bins=15)[0]
        flake_shape_contour_hist = flake_shape_contour_hist / flake_shape_contour_hist.sum()
        flake_shape_fracdim = self.MatSegUtils.fractal_dimension(flake_region)
        flake_i['flake_shape_fea'] = np.array([flake_shape_len_area_ratio] + list(flake_shape_contour_hist) + [flake_shape_fracdim])
        
        flake_i['flake_shape_fea_names'] = ['P/A'] + ['hist {}'.format(i) for i in range(15)] + ['fractal dimension']


        # Calculate color features of the flake
        if run_args['background']:
            flake_i['flake_contrast'] = (im_gray[flake_region>0].mean() - bk_gray[flake_region>0].mean())/255.0
        else:
            background_region = (image_labelmap==0)
            flake_i['flake_contrast'] = (im_gray[flake_region>0].mean() - im_gray[background_region].mean())/255.0
        
        inner_region = cv2.erode(flake_region, kernel, iterations=1)

        if run_args['background']:
            flake_color_fea = [im_gray[flake_region>0].mean() - bk_gray[flake_region>0].mean(), 
                            im_hsv[flake_region>0, 2].mean() - bk_hsv[flake_region>0, 2].mean()] + \
                            [im_gray[flake_region>0].mean(), im_gray[flake_region>0].std()] + \
                            list(im_hsv[flake_region>0].mean(0)) + list(im_hsv[flake_region>0].std(0)) + \
                            list(im_rgb[flake_region>0].mean(0)) + list(im_rgb[flake_region>0].std(0))
        else:
            flake_color_fea = [im_gray[flake_region>0].mean() - im_gray[background_region].mean(), 
                            im_hsv[flake_region>0, 2].mean() - im_hsv[background_region, 2].mean()] + \
                            [im_gray[flake_region>0].mean(), im_gray[flake_region>0].std()] + \
                            list(im_hsv[flake_region>0].mean(0)) + list(im_hsv[flake_region>0].std(0)) + \
                            list(im_rgb[flake_region>0].mean(0)) + list(im_rgb[flake_region>0].std(0))
            
        # flake_color_entropy = entropy(im_gray[flake_region>0].astype('uint8'), disk(5))
        flake_color_entropy = cv2.calcHist([im_gray[flake_region>0].astype('uint8')],[0],None,[256],[0,256])
        flake_color_entropy = entropy(flake_color_entropy, base=2)
        flake_inner_color_fea = [0] * 16
        flake_inner_color_entropy = 0
        if inner_region.sum() > 0:
            if run_args['background']:
                flake_inner_color_fea = [im_gray[inner_region>0].mean() - bk_gray[inner_region>0].mean(), 
                            im_hsv[inner_region>0, 2].mean() - bk_hsv[inner_region>0, 2].mean()] + \
                            [im_gray[inner_region>0].mean(), im_gray[inner_region>0].std()] + \
                            list(im_hsv[inner_region>0].mean(0)) + list(im_hsv[inner_region>0].std(0)) + \
                            list(im_rgb[inner_region>0].mean(0)) + list(im_rgb[inner_region>0].std(0))
            else:
                flake_inner_color_fea = [im_gray[inner_region>0].mean() - im_gray[background_region].mean(), 
                            im_hsv[inner_region>0, 2].mean() - im_hsv[background_region, 2].mean()] + \
                            [im_gray[inner_region>0].mean(), im_gray[inner_region>0].std()] + \
                            list(im_hsv[inner_region>0].mean(0)) + list(im_hsv[inner_region>0].std(0)) + \
                            list(im_rgb[inner_region>0].mean(0)) + list(im_rgb[inner_region>0].std(0))

            # flake_inner_color_entropy = entropy(im_gray[inner_region>0].astype('uint8'), disk(5))
            flake_inner_color_entropy = cv2.calcHist([im_gray[inner_region>0].astype('uint8')],[0],None,[256],[0,256])
            flake_inner_color_entropy = entropy(flake_inner_color_entropy, base=2)
        flake_i['flake_color_fea'] = np.array(flake_color_fea + [flake_color_entropy] + flake_inner_color_fea + [flake_inner_color_entropy])
        
        
        feature_names_color = [
            'g contrast',
            'v contrast',
            'gray',
            'gray std',
            'H',
            'S',
            'V',
            'H std',
            'S std',
            'V std',
            'R',
            'G',
            'B',
            'R std',
            'G std',
            'B std',
            'entropy'
            ]
        feature_names_color = feature_names_color + ['{}_inner'.format(f) for f in feature_names_color]
        flake_i['flake_color_fea_names'] = feature_names_color

        
        return flake_i
    
    
    
    def total_length(self, line):
        '''Computes the length of a sequence of line segments, given as a sequence of (x,y) pairs.
        The curve is 'closed', meaning that the length between the first and last points is included
        in the sum.'''
        
        line2 = np.roll(line, 1, axis=0)
        segment_lengths = np.sqrt(np.square(line[:,0]-line2[:,0])+np.square(line[:,1]-line2[:,1]))
        
        return np.sum(segment_lengths)
                
    def polygon_area(self, line):
        # Simple implementation of Shoelace formula:
        #  https://en.wikipedia.org/wiki/Shoelace_formula
        # Code from:
        #  https://stackoverflow.com/questions/24467972/calculate-area-of-polygon-given-x-y-coordinates
        x = line[:,0]
        y = line[:,1]
        return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))
            
    def colored_objects(self, labeled_array, outfile, **run_args):
        # Colored image
        im = PIL.Image.fromarray( np.uint8(labeled_array) )
        im = im.convert('RGB')
        pix = im.load()
        
        color_list = [ (1,0,0), (0,1,0), (0,0,1), (1,1,0), (1,0,1),(0,1,1),(1,1,1),]
        color_list2 = [ (0.7*c[0], 0.7*c[1], 0.7*c[2]) for c in color_list ]
        color_list3 = [ (0.5*c[0], 1.0*c[1], 1.0*c[2]) for c in color_list ]
        color_list4 = [ (1.0*c[0], 0.5*c[1], 1.0*c[2]) for c in color_list ]
        color_list5 = [ (1.0*c[0], 1.0*c[1], 0.5*c[2]) for c in color_list ]
        color_list6 = [ (1.0*c[0], 0.7*c[1], 0.5*c[2]) for c in color_list ]
        color_list7 = [ (1.0*c[0], 0.5*c[1], 0.7*c[2]) for c in color_list ]
        color_list8 = [ (0.7*c[0], 1.0*c[1], 0.5*c[2]) for c in color_list ]
        color_list9 = [ (0.5*c[0], 1.0*c[1], 0.7*c[2]) for c in color_list ]
        color_list10 = [ (0.7*c[0], 0.5*c[1], 1.0*c[2]) for c in color_list ]
        color_list11 = [ (0.5*c[0], 0.7*c[1], 1.0*c[2]) for c in color_list ]
        color_list = color_list + color_list2 + color_list3 + color_list4 + color_list5 + color_list6 + color_list7 + color_list8 + color_list9 + color_list10 + color_list11
        
        h, w = labeled_array.shape
        for ix in range(w):
            for iy in range(h):
                object_index = labeled_array[iy,ix]
                if object_index==0:
                    c = (0, 0, 0)
                else:
                    c = color_list[ int(object_index%len(color_list)) ]
                    c = ( int(c[0]*255), int(c[1]*255), int(c[2]*255) )
                pix[ix,iy] = c            
        
        im.save(outfile)        



def get_in_range(data, im_contrast, image_contrast_trim=None, **run_args):
    
    if image_contrast_trim is not None:
        image_contrast_trim = np.clip(image_contrast_trim, 0, 0.95)
        
        avg = np.average(data.data_rgb)
        avg /= 255
        amt = image_contrast_trim
        im_contrast = ( avg*amt , 1.0-(1.0-avg)*amt )

    in_range = ( im_contrast[0]*255, im_contrast[1]*255 )    
    
    return in_range


class flake_images(Protocol):
    
    def __init__(self, name='flake_images', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
                        'verbosity' : 3 ,
                        'bbox_pad' : 0.5,
                        'image_contrast' : (0, 1),
                        'image_contrast2' : None,
                        'image_contrast_trim' : None,
                        'saved_dir' : './find_flakes/',
                        'overlays' : 10,
                        }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, data, output_dir, **run_args):
        
        results = {}
        
        if 'scale' not in run_args:
            run_args['scale'] = np.average((data.x_scale, data.y_scale)) # um/pixel
        run_args['scale2'] = run_args['scale']*run_args['scale'] # um^2/pixel^2
        
            
        with open(os.path.join(run_args['saved_dir'], data.name+'.pkl'), 'rb') as fin:
            saved_information = pickle.load(fin)
            
        
            
        for flake_i in saved_information['flakes']:
            self.plot_flake(data, flake_i, output_dir=output_dir, **run_args)
        
        results['flakes_plotted'] = len(saved_information['flakes'])
        
        return results
    
    
    
    def plot_flake(self, data, flake_i, output_dir, **run_args):
        
        plt.rcParams['xtick.labelsize'] = 15
        h, w, c = data.data_rgb.shape
        
        y, x = flake_i['center_of_mass']
        if run_args['verbosity']>=4:
            print('  Flake id{}, (x, y) = ({:d}, {:d})'.format(flake_i['index'], int(x), int(y)))

        ta = 0.5
        ba = 0.5 # Modify this to expand/contract the bottom half of the layout
        wa, ha = 1, ta+ba
        aspect = wa/ha
        self.fig = plt.figure( figsize=(10,10/aspect), facecolor='white' )
        #                               x    y    xw    yh
        self.ax1 = self.fig.add_axes( [0, ba/ha, 0.5, ta/ha] )
        self.ax2 = self.fig.add_axes( [0.5, ba/ha, 0.5, ta/ha] )
        
        spacer = 0.1
        self.ax3 = self.fig.add_axes( [0, 0, 0.5, (ba/ha)*(1-spacer)], zorder=-10 )
        #self.ax4 = self.fig.add_axes( [0.5, 0, 0.5, (ba/ha)*(1-spacer)], zorder=-10 )
        
        self.ax5 = self.fig.add_axes( [0.58, 0.06, 0.41, (ba/ha)*(1-spacer)*0.7], zorder=-10 )
        
        
        y1, y2, x1, x2 = flake_i['bbox']

        # Make the crop border a bit bigger than the flake bounding box
        box_size = (1+run_args['bbox_pad'])*max( abs(x2-x1), abs(y2-y1) )
        x1p = int(np.clip((x1+x2)*0.5 - box_size/2, 0, w))
        x2p = int(np.clip((x1+x2)*0.5 + box_size/2, 0, w))
        y1p = int(np.clip((y1+y2)*0.5 - box_size/2, 0, h))
        y2p = int(np.clip((y1+y2)*0.5 + box_size/2, 0, h))
        box = y1p, y2p, x1p, x2p
        
        self.plot_ax1(data, box, **run_args)
        
        self.plot_ax2(data, box, flake_i, **run_args)
        
        self.plot_ax3(flake_i, **run_args)
        
        self.plot_ax5(flake_i, **run_args)
        

        if 'show' in run_args and run_args['show']:
            plt.show()
        outfile = self.get_outfile(data.name, output_dir, ext='-flake{:03d}.png'.format(flake_i['index']))
        self.fig.savefig(outfile)
        plt.close(self.fig.number)

    
    def plot_ax1(self, data, box, **run_args):
        
        y1p, y2p, x1p, x2p = box
        
        #in_range = ( run_args['image_contrast'][0]*255, run_args['image_contrast'][1]*255 )
        in_range = get_in_range(data, run_args['image_contrast'], **run_args)
        
        flake = data.data_rgb[y1p:y2p , x1p:x2p, :]
        flake = skimage.exposure.rescale_intensity(flake, in_range=in_range, out_range='dtype')
        self.ax1.imshow(flake)
        #xi, xf, yi, yf = self.ax1.axis()
        #self.ax1.axis([0, box_size, box_size, 0])
        self.ax1.axes.get_yaxis().set_visible(False)
        self.ax1.set_xlabel('pixels', fontsize=20)
        plt.setp(self.ax1.get_xticklabels()[-1], visible=False)
        plt.setp(self.ax1.get_xticklabels()[-2], visible=False)


    def plot_ax2(self, data, box, flake_i, **run_args):
        
        scale = run_args['scale']
        
        y1, y2, x1, x2 = flake_i['bbox']
        y1p, y2p, x1p, x2p = box
        box_size = max(abs(y1p-y2p), abs(x1p-x2p))
        
        if run_args['image_contrast2'] is None:
            run_args['image_contrast2'] = run_args['image_contrast']
        #in_range2 = ( run_args['image_contrast2'][0]*255, run_args['image_contrast2'][1]*255 )
        if run_args['image_contrast_trim'] is not None:
            in_range2 = get_in_range(data, run_args['image_contrast2'], image_contrast_trim=run_args['image_contrast_trim']*1.75)
        else:
            in_range2 = get_in_range(data, run_args['image_contrast2'])
        
        flake = data.data_rgb[y1p:y2p , x1p:x2p, :]
        flake = skimage.exposure.rescale_intensity(flake, in_range=in_range2, out_range='dtype')
        
        he, we, c = flake.shape # Actual shape of region
        if he==we:
            extent = [0, box_size*scale, box_size*scale, 0]
        elif he>we:
            extent = [0, box_size*scale*we/he, box_size*scale, 0]
        else:
            extent = [0, box_size*scale, box_size*scale*he/we, 0]
        self.ax2.imshow(flake, extent=extent, aspect='equal')
        self.ax2.axes.get_yaxis().set_visible(False)
        self.ax2.set_xlabel(r'$\mathrm{ \mu m }$', fontsize=20)
        
        xi, xf, yi, yf = self.ax2.axis()
        y, x = flake_i['center_of_mass']
        self.ax2.text(xi, yf, '({}, {})'.format(int(x),int(y)), size=12, weight='bold', horizontalalignment='left', verticalalignment='top')

        if run_args['overlays']>=1:
            c = flake_i['contour']
            xs = (c[:,0] - x1p)*scale
            ys = (c[:,1] - y1p)*scale
            self.ax2.plot(xs, ys, '-', linewidth=4.0, color='r', dashes=[4,1], alpha=0.4)
        if run_args['overlays']>=7:
            c = flake_i['convex_hull']
            xs = (c[:,1] - x1p)*scale
            ys = (c[:,0] - y1p)*scale
            self.ax2.plot(xs, ys, '-', linewidth=3.0, color='g')

        if run_args['overlays']>=2:
            rect = patches.Rectangle( ((x1-x1p)*scale, (y1-y1p)*scale), (x2-x1)*scale, (y2-y1)*scale, linewidth=2.5, edgecolor='orange', facecolor='none', alpha=0.5)
            self.ax2.add_patch(rect)
            
        if run_args['overlays']>=3:
            # Circle overlay denoting size
            y, x = flake_i['center_of_mass']
            x = (x-x1p)*scale
            y = (y-y1p)*scale
            size = flake_i['radius_um']
            circ = patches.Circle(xy=(x,y), radius=size, linewidth=2.0, edgecolor='r', facecolor='none', alpha=0.3)
            self.ax2.add_patch(circ)
            
        if run_args['overlays']>=4:
            # Cross hair
            rect = patches.Rectangle( (x-size/2, y), size, 0, linewidth=2.0, edgecolor='r', facecolor='none', alpha=0.3) # Horizontal bar
            self.ax2.add_patch(rect)
            rect = patches.Rectangle( (x, y-size/2), 0, size, linewidth=2.0, edgecolor='r', facecolor='none', alpha=0.3) # Vertical bar
            self.ax2.add_patch(rect)

        
    def plot_ax3(self, flake_i, **run_args):
        self.ax3.axis([0, 1, 0, 1])
        self.ax3.axis('off')
        
        #self.ax3.text(0.5, 0.5, 'blah')
        
        
        # Row 1 (headers)
        texts = ['', '', '$P_i \, (\mathrm{ \mu m })$', '$P_i/P_{\mathrm{circ}}$', '$P_{\mathrm{pix}}/P_i$', '$A \, (\mathrm{ \mu m ^2})$', '%', '$r \, (\mathrm{ \mu m})$']
        for icol, text in enumerate(texts):
            self.table3_cell(0, icol, text)
        
        # Column 1 (names)
        texts = ['', 'circle', 'pixel map', 'contour', 'convex hull', 'bounding box']
        for irow, text in enumerate(texts):
            self.table3_cell(irow, 0, text, horizontalalignment='left')
        fills = ['none', 'r', 'r', 'r', 'g', 'orange']
        for irow, fill in enumerate(fills):
            self.table3_cell(irow, 1, ' ', facecolor=fill, horizontalalignment='left')
            
            
            
        # Perimeter
        icol = 2
        circumference = flake_i['radius_um']*2*np.pi
        self.table3_cell(1, icol, '{:,.1f}'.format( circumference ) ) # circle
        self.table3_cell(2, icol, '{:,.1f}'.format( flake_i['perimeter_um'] ) ) # pixel map
        self.table3_cell(3, icol, '{:,.1f}'.format( flake_i['contour_perimeter_um'] ) ) # contour
        self.table3_cell(4, icol, '{:,.1f}'.format( flake_i['convex_perimeter_um'] ) ) # convex hull
        y1, y2, x1, x2 = flake_i['bbox']
        P_um = 2*(abs(y2-y1)+abs(x2-x1))*run_args['scale']
        self.table3_cell(5, icol, '{:,.1f}'.format(P_um) ) # bounding box

        # P_i/P_circ
        icol += 1
        self.table3_cell(1, icol, '{:,.1f}'.format( 1.0 ), weight='light' )
        self.table3_cell(2, icol, '{:,.1f}'.format( flake_i['perimeter_um']/circumference ) ) # pixel map
        self.table3_cell(3, icol, '{:,.1f}'.format( flake_i['contour_perimeter_um']/circumference ), weight='bold' ) # contour
        self.table3_cell(4, icol, '{:,.1f}'.format( flake_i['convex_perimeter_um']/circumference ) ) # convex hull
        self.table3_cell(5, icol, '{:,.1f}'.format(P_um/circumference) ) # bounding box            

        # P_pix/P_i
        icol += 1
        P_pix = flake_i['perimeter_um']
        self.table3_cell(1, icol, '{:,.1f}'.format( P_pix/circumference ) )
        self.table3_cell(2, icol, '{:,.1f}'.format( 1.0 ), weight='light' ) # pixel map
        self.table3_cell(3, icol, '{:,.1f}'.format( P_pix/flake_i['contour_perimeter_um'] ) ) # contour
        self.table3_cell(4, icol, '{:,.1f}'.format( P_pix/flake_i['convex_perimeter_um'] ), weight='bold' ) # convex hull
        self.table3_cell(5, icol, '{:,.1f}'.format( P_pix/P_um ) ) # bounding box            
            
            
        # Area
        icol = 5
        self.table3_cell(1, icol, '{:,.1f}'.format(flake_i['size_um']) ) # circle
        self.table3_cell(2, icol, '{:,.1f}'.format(flake_i['size_um']) ) # pixel map
        self.table3_cell(3, icol, '{:,.1f}'.format(flake_i['contour_size_um']) ) # contour
        self.table3_cell(4, icol, '{:,.1f}'.format(flake_i['convex_size_um']) ) # convex hull
        y1, y2, x1, x2 = flake_i['bbox']
        A_um = (abs(y2-y1)*abs(x2-x1))*run_args['scale2']
        self.table3_cell(5, icol, '{:,.1f}'.format(A_um) ) # bounding box
        
        # Area %
        icol += 1
        f = 100.*flake_i['size_um']
        self.table3_cell(1, icol, '{:.0f}%'.format(100), weight='light' ) # circle
        self.table3_cell(2, icol, '{:.0f}%'.format(100), weight='light' ) # pixel map
        self.table3_cell(3, icol, '{:.0f}%'.format(f/flake_i['contour_size_um']), weight='light' ) # contour
        self.table3_cell(4, icol, '{:.0f}%'.format(f/flake_i['convex_size_um']), weight='bold' ) # convex hull
        self.table3_cell(5, icol, '{:.0f}%'.format( f/A_um ) ) # bounding box

        # r
        icol += 1
        self.table3_cell(1, icol, '{:,.1f}'.format( np.sqrt(flake_i['size_um']/np.pi) ) ) # circle
        self.table3_cell(2, icol, '{:,.1f}'.format( np.sqrt(flake_i['size_um']/np.pi) ), weight='bold') # pixel map
        self.table3_cell(3, icol, '{:,.1f}'.format( np.sqrt(flake_i['contour_size_um']/np.pi) ) ) # contour
        self.table3_cell(4, icol, '{:,.1f}'.format( np.sqrt(flake_i['convex_size_um']/np.pi) ) ) # convex hull
        self.table3_cell(5, icol, '{:,.1f}'.format( np.sqrt(A_um/np.pi) ) ) # bounding box
        


        # Flake color features (34 features)
        #                                      16 els.[0..15]         1 el.[16]          16 els.[17..32]                1 el.[33]
        #flake_i['flake_color_fea'] = np.array(flake_color_fea + [flake_color_entropy] + flake_inner_color_fea + [flake_inner_color_entropy])        
            #flake_color_fea = [im_gray[flake_region>0].mean() - bk_gray[flake_region>0].mean(),            # 0
                            #im_hsv[flake_region>0, 2].mean() - bk_hsv[flake_region>0, 2].mean()] + \       # 1
                            #[im_gray[flake_region>0].mean(), im_gray[flake_region>0].std()] + \            # 2, 3
                            #list(im_hsv[flake_region>0].mean(0)) + list(im_hsv[flake_region>0].std(0)) + \ # 4,5,6        7,8,9
                            #list(im_rgb[flake_region>0].mean(0)) + list(im_rgb[flake_region>0].std(0))     # 10,11,12   13,14,15
        idx_inner = 17
        
        # Row 1 (headers)
        texts = ['', 'contour', 'inner', 'diff.']
        for icol, text in enumerate(texts):
            self.table3b_cell(0, icol, text)
        
        # Column 1 (names)
        texts = ['', 'contrast', 'gray', 'H', 'S', 'V', 'R', 'G', 'B', 'entropy']
        fills = ['none', '0.8', '0.8', 'lightblue', 'lightblue', '0.9', 'r', 'g', 'b', 'none']
        for irow, text in enumerate(texts):
            self.table3b_cell(irow, 0, text, facecolor=fills[irow])
        
        self.table3b_row(1, flake_i['flake_contrast'], flake_i['flake_color_fea'][idx_inner]/255, decimals=True ) # contrast
        self.table3b_row(2, flake_i['flake_color_fea'][2], flake_i['flake_color_fea'][2+idx_inner], flake_i['flake_color_fea'][3], flake_i['flake_color_fea'][3+idx_inner], decimals=False ) # gray
        
        #            H S V R G B
        for irow in [3,4,5,6,7,8]:
            if irow<=5:
                i = irow+1
                decimals = True
            else:
                i = irow+4
                decimals = False
                
            self.table3b_row(irow, flake_i['flake_color_fea'][i], flake_i['flake_color_fea'][i+idx_inner], flake_i['flake_color_fea'][i+3], flake_i['flake_color_fea'][i+idx_inner+3], decimals=decimals )
        
        self.table3b_row(9, flake_i['flake_color_fea'][16], flake_i['flake_color_fea'][16+idx_inner], decimals=True ) # entropy
        
        
    def table3_cell(self, irow, icol, text, ts=10, weight='normal', facecolor='none', horizontalalignment='center', verticalalignment='center'):
        xo, yo = 0, 0.95
        
        row_spacing = 0.04
        col_widths = [0.20, 0.03, 0.12, 0.12, 0.12, 0.15, 0.12, 0.1]
        
        # Center of cell
        x = xo + col_widths[icol]*0.5 + np.sum(col_widths[:icol])
        y = yo - row_spacing*(irow+0.5)

        if text!='':
            rect = patches.Rectangle( (x-col_widths[icol]/2, y-row_spacing/2), col_widths[icol], row_spacing, linewidth=0.2, edgecolor='0.1', facecolor=facecolor, alpha=0.6)
            self.ax3.add_patch(rect)

        
        if horizontalalignment=='left':
            x -= col_widths[icol]*0.5
        
        self.ax3.text(x, y, text, size=ts, weight=weight, horizontalalignment=horizontalalignment, verticalalignment=verticalalignment)
        
    def table3b_cell(self, irow, icol, text, ts=10, weight='normal', facecolor='none', horizontalalignment='center', verticalalignment='center'):
        xo, yo = 0, 0.65
        
        row_spacing = 0.04
        col_widths = [0.12, 0.17, 0.17, 0.12]
        
        # Center of cell
        x = xo + col_widths[icol]*0.5 + np.sum(col_widths[:icol])
        y = yo - row_spacing*(irow+0.5)

        if text!='':
            rect = patches.Rectangle( (x-col_widths[icol]/2, y-row_spacing/2), col_widths[icol], row_spacing, linewidth=0.2, edgecolor='0.1', facecolor=facecolor, alpha=0.6)
            self.ax3.add_patch(rect)

        
        if horizontalalignment=='left':
            x -= col_widths[icol]*0.5
        
        self.ax3.text(x, y, text, size=ts, weight=weight, horizontalalignment=horizontalalignment, verticalalignment=verticalalignment)
                
    def table3b_row(self, irow, val_a, val_b, std_a=None, std_b=None, decimals=False):
        
        # Unpack the values
        # (For some reason returned values are in one-element arrays now. Maybe due to a change in a call signature somewhwere in cv2, or a change from Python2.7 to Python3, or some other change in a library.)
        if isinstance(val_a, np.ndarray) and len(val_a)==1:
            val_a = val_a[0]
        if isinstance(val_b, np.ndarray) and len(val_b)==1:
            val_b = val_b[0]
        
        f = '{:.3f}' if decimals else '{:.0f}'
        f2 = '{:.2f}±{:.2f}' if decimals else '{:.0f}±{:.0f}'
        
        if std_a is None:
            self.table3b_cell(irow, 1, f.format( val_a ) )
        else:
            self.table3b_cell(irow, 1, f2.format( val_a, std_a ) )
        if std_b is None:
            self.table3b_cell(irow, 2, f.format( val_b ) )
        else:
            self.table3b_cell(irow, 2, f2.format( val_b, std_b ) )
            
        d = val_a-val_b
        if abs(d)<1:
            self.table3b_cell(irow, 3, '{:.3f}'.format(d) )
        else:
            self.table3b_cell(irow, 3, '{:.1f}'.format(d) )
        
        #self.table3b_cell(irow, 1, '{:,.1f} ± {}'.format( np.sqrt(A_um/np.pi) ) )
        
    
    def plot_ax5(self, flake_i, **run_args):
            
        #flake_i['flake_shape_fea'] = np.array([flake_shape_len_area_ratio] + list(flake_shape_contour_hist) + [flake_shape_fracdim])
        contour_hist = flake_i['flake_shape_fea'][1:-1]

        self.ax5.tick_params(axis='both', which='major', labelsize=10)
        self.ax5.bar(x=range(len(contour_hist)), height=contour_hist)
        self.ax5.set_xlabel('bin', fontsize=10)
        self.ax5.set_ylabel('chord distribution', fontsize=10)
    
    
    
    
    
