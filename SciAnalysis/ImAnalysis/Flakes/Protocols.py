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
import skimage.measure as measure
import cv2


try:
    from MaterialSegmentation import utils
except:
    code_PATH='/home/kyager/current/code/MaterialSegmentation/main/'
    code_PATH in sys.path or sys.path.append(code_PATH)
    import utils




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
                        'dpi' : 'auto' ,
                        'resize' : 1.0 ,
                        'bbox_expanded' : 0.1 ,
                        'overlays' : 4 , # Larger number adds more annotations
                        }
        self.run_args.update(kwargs)
        

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
        im_gray = np.array(image.convert('L', (0.2989, 0.5870, 0.1140, 0))).astype('float')
        # im_hsv = np.array(image.convert('HSV')).astype('float')
        # Alternate HSV conversion (consistent with Matlab)
        im_hsv = skimage.color.rgb2hsv(im_rgb)
        im_hsv[:,:,2] = im_hsv[:,:,2]/255.0

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
        res_map, image_labelmap, flake_centroids, flake_sizes, num_flakes = utils.perform_robustfit_multichannel(im_hsv, im_gray, run_args['image_threshold'], run_args['size_threshold'])
        if run_args['verbosity']>=4:
            print('    {} flakes'.format(num_flakes))
        
        if run_args['background']:
            if run_args['verbosity']>=4:
                print('  Segmenting background')
            _, _, bk_flake_centroids, _, n = utils.perform_robustfit(bk_gray, 10, 0)
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
            scipy.misc.toimage(res_map, cmin=0, cmax=255).save(outfile)

            arr = image_labelmap>0
            outfile = self.get_outfile(data.name, output_dir, ext='-binary.png')
            scipy.misc.toimage(arr, cmin=0, cmax=1).save(outfile)
            
            arr = np.dstack([arr]*3)
            im = np.where(arr, im_rgb, np.zeros_like(im_rgb))
            outfile = self.get_outfile(data.name, output_dir, ext='-flakes.png')
            scipy.misc.toimage(im, cmin=0, cmax=255).save(outfile)
            
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
        
        if 'image_contrast' is not None:
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
            
            flake_i = self.flake_plot(flake_i, scale, **run_args)

            flakes.append(flake_i)

        

        # Finalize
        ########################################
        
        # Save image
        #xi, xf, yi, yf = self.ax.axis()
        self.ax.axis( [0, w, h, 0] )
        
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
        
        
        # Flake contour
        _, flake_contour, _ = cv2.findContours(flake_region, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
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
        flake_shape_fracdim = utils.fractal_dimension(flake_region)
        flake_i['flake_shape_fea'] = np.array([flake_shape_len_area_ratio] + list(flake_shape_contour_hist) + [flake_shape_fracdim])


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
