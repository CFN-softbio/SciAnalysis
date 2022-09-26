#!/usr/bin/python
# -*- coding: utf-8 -*-
# vi: ts=4 sw=4

from ..Protocols import *

#import sys
#import pickle

from PIL import Image
#import matplotlib.patches as patches

from scipy import ndimage
#from scipy.stats import entropy
#from scipy.spatial.distance import cdist
import skimage
import skimage.exposure
import skimage.measure as measure
#import cv2



class find_dots(Protocol):
    
    def __init__(self, name='find_dots', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
                        'image_threshold' : 10 , # Threshold to define the foreground in the residual map
                        'size_threshold' : 0 , # Size (in pixels) to consider a contiguous region to be a                         
                        'show' : False ,
                        'dpi' : 'auto' ,
                        'resize' : 1.0 ,
                        'hist_bins' : 100 ,
                        }
        self.run_args.update(kwargs)


    #def output_exists(self, name, output_dir):
        #if 'file_extension' in self.run_args:
            #ext = self.run_args['file_extension']
        #else:
            #ext = None
        #outfile = self.get_outfile(name, output_dir, ext=ext)
        #return os.path.isfile(outfile)

    def output_exists(self, name, output_dir):
        outfile = self.get_outfile('{}/properties'.format(name), output_dir, ext='.npy')
        return os.path.isfile(outfile)

        

    @run_default
    def run(self, data, output_dir, **run_args):
        
        output_dir = os.path.join(output_dir, data.name)
        make_dir(output_dir)
        
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
        #im_hsv = skimage.color.rgb2hsv(im_rgb)
        #im_hsv[:,:,2] = im_hsv[:,:,2]/255.0 # vals in range [0..1]
        
        im_hsv = mpl.colors.rgb_to_hsv(im_rgb/255)

        h, w, c = im_rgb.shape
        if run_args['verbosity']>=5:
            print('    {}×{} = {:.2f} Mpix; {:.1f}μm × {:.1f}μm = {:.3f} mm^2'.format(w, h, w*h/1e6, w*scale, h*scale, w*h*scale2/1e6))


        if True: # abla for testing
            # Identify dots
            ########################################
            im_data = im_hsv[:,:,2] # Use brightness channel
            scale2 = data.x_scale*data.y_scale # (µm/pixel)^2
            
            results, labeled_array, properties = self._find_objects(im_data, scale2, output_dir, results, **run_args)


            # Analyze dots
            ########################################
            
            results = self._particle_stats(labeled_array, scale2, properties, output_dir, results, **run_args)

            results, properties = self._particle_color(im_rgb, im_hsv, labeled_array, scale2, properties, output_dir, results, **run_args)

            results, properties = self._particle_shape(labeled_array, scale2, properties, output_dir, results, **run_args)
            
            outfile = self.get_outfile('properties', output_dir, ext='.npy', ir=False)
            np.save(outfile, properties, allow_pickle=True)
        
        
        outfile = self.get_outfile('properties', output_dir, ext='.npy', ir=False)
        properties = np.load(outfile, allow_pickle=True).item()


        # Output image
        ########################################
        if True: # abla for testing
            self._plot_image(im_rgb, labeled_array, properties, output_dir, **run_args)
        
        self._plot_feature_space(properties, output_dir, **run_args)
        
        return results
    
    
    

        


    def _find_objects(self, im_data, scale2, output_dir, results, **run_args):
        
        im_thresh = np.where( im_data>run_args['threshold'], 1, 0 )
        if 'invert' in run_args and run_args['invert']:
            im_thresh = 1-im_thresh
            
            
        if run_args['verbosity']>=5:
            im = PIL.Image.fromarray( np.uint8(im_thresh*255.0) )
            outfile = self.get_outfile('thresholded', output_dir, ext='.png', ir=True)
            im.save(outfile)
            
            

        # Identify particles positions
        if 'diagonal_detection' in run_args and run_args['diagonal_detection']:
            #s = [[1,1,1],
            #    [1,1,1],
            #    [1,1,1]]
            s = ndimage.generate_binary_structure(2,2)
        else:
            s = [[0,1,0],
                [1,1,1],
                [0,1,0]]
        
        labeled_array, num_features = ndimage.measurements.label(im_thresh, structure=s)
        results['num_particles_raw'] = num_features
        
        properties = {
            'label': [],
            'area_pix': [],
            'area_um2': [],
            'radius_um': [],
            }
        
        # Remove objects not meeting size criteria
        print('  Computing particle sizes')
        i_max = np.max(labeled_array)+1
        n_exclude = 0
        for i in range(i_max):
            if i%100==0:
                print('    {}/{} = {:.1f}%'.format(i, i_max, 100.*i/i_max))
            particle = (labeled_array==i).astype(np.uint8)
            area_pix = np.sum(particle)
            area_um2 = area_pix*scale2 # µm^2
            radius_um = np.sqrt(area_um2/np.pi) # µm

            # Select only desired particles
            exclude = (i==0) # Don't include the background
            if 'area_min' in run_args and area_um2<run_args['area_min']:
                exclude = True
            if 'area_max' in run_args and area_um2>run_args['area_max']:
                exclude = True
            if 'radius_min' in run_args and radius_um<run_args['radius_min']:
                exclude = True
            if 'radius_max' in run_args and radius_um>run_args['radius_max']:
                exclude = True 
            
            if exclude:
                idx = np.where(labeled_array==i)
                labeled_array[idx] = 0
                n_exclude += 1
            else:
                properties['label'].append(i)
                properties['area_pix'].append(area_pix)
                properties['area_um2'].append(area_um2)
                properties['radius_um'].append(radius_um)

        properties['label'] = np.asarray(properties['label'])
        properties['area_pix'] = np.asarray(properties['area_pix'])
        properties['area_um2'] = np.asarray(properties['area_um2'])
        properties['radius_um'] = np.asarray(properties['radius_um'])

                
        results['num_particles'] = num_features - n_exclude

                
        if run_args['verbosity']>=6:
            # Colored image
            print('  Coloring dots')
            
            im = PIL.Image.fromarray( np.uint8(im_data*255.0) )
            im = im.convert('RGB')
            pix = im.load()
            
            h, w = im_data.shape
            for ix in range(w):
                for iy in range(h):
                    object_index = labeled_array[iy,ix]
                    if object_index==0:
                        c = (0, 0, 0)
                    else:
                        c = color_list[ object_index%len(color_list) ]
                        c = ( int(c[0]*255), int(c[1]*255), int(c[2]*255) )
                    pix[ix,iy] = c            
            
            
            outfile = self.get_outfile('colored', output_dir, ext='.png', ir=True)
            im.save(outfile)
            
        
        if run_args['verbosity']>=7:
            # Boundary image
            print('  Boundary image')
            
            im = PIL.Image.fromarray( np.uint8(im_data*255.0) )
            im = im.convert('RGB')
            pix = im.load()
            
            c = ( 1*255, 0*255, 0*255 )
            h, w = im_data.shape
            for ix in range(w-1):
                for iy in range(h-1):
                    
                    num_zeros = np.bincount( labeled_array[iy:iy+2,ix:ix+2].flatten() )[0]
                    if not (num_zeros==0 or num_zeros==4):
                        pix[ix,iy] = c            
            

            outfile = self.get_outfile('boundaries', output_dir, ext='.png', ir=True)
            im.save(outfile)                
                
                
        return results, labeled_array, properties



    def _particle_stats(self, labeled_array, scale2, properties, output_dir, results, **run_args):
        
        # Statistics on particles that have been found
        bins = np.bincount(labeled_array.flatten('C'))
        h, w = labeled_array.shape
        total_pixels = w*h
        background_pixels = bins[0]
        particle_pixels = total_pixels - background_pixels
        coverage = particle_pixels*1./(total_pixels*1.)
        results['coverage'] = coverage
        if run_args['verbosity']>=4:
            print('  Stats:')
            print('    {} particles'.format(results['num_particles']))
            print('    Particle coverage: {:.1f}%'.format(coverage*100.))
        
        
        # Remove 'particles' of zero size
        idx = np.nonzero(bins)
        bins = bins[idx]
        # Remove the 'surrounding field' (index 0)
        bins = bins[1:]        
        
        # Convert to physical sizes
        particle_sizes = bins*scale2 # µm^2
        particle_radii = np.sqrt(particle_sizes/np.pi) # µm
        
        results['area_average'] = np.average(particle_sizes)
        results['area_std'] = np.std(particle_sizes)
        results['area_median'] = np.median(particle_sizes)
        
        results['radius_average'] = np.average(particle_radii)
        results['radius_std'] = np.std(particle_radii)
        results['radius_median'] = np.median(particle_radii)        
        
        
        if run_args['verbosity']>=4:
            self._plot_histograms(particle_sizes, particle_radii, output_dir, results, **run_args)
        
        
        return results
    
    
    def _plot_histograms(self, particle_sizes, particle_radii, output_dir, results, **run_args):
        
        class DataHistogram_current(DataHistogram):
            def _plot_extra(self, **plot_args):
                
                xi, xf, yi, yf = self.ax.axis()
                yf *= 1.2
                self.ax.axis([xi, xf, yi, yf])
                
                self._plot_stats(**plot_args)
                
                lm_result, fit_line, fit_line_e = self.fit_peak(self)
                self.ax.axvline(lm_result.params['x_center'].value, color='r', linewidth=2.0)
                self.ax.plot(fit_line_e.x, fit_line_e.y, 'r', linewidth=2.0)
                
                
            def fit_peak(self, line, **run_args):
                import lmfit
                
                def model(v, x):
                    m = v['prefactor']*np.exp( -np.square(x-v['x_center'])/(2*(v['sigma']**2)) )
                    return m
                
                def func2minimize(params, x, data):
                    
                    v = params.valuesdict()
                    m = model(v, x)
                    
                    return m - data
                
                params = lmfit.Parameters()
                params.add('prefactor', value=np.max(line.y), min=0, max=np.max(line.y)*2.0)
                params.add('x_center', value=self.mean, min=np.min(line.x), max=np.max(line.x))
                params.add('sigma', value=self.std, min=0)
                
                lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
                
                fit_x = np.linspace(np.min(line.x), np.max(line.x), num=200)
                fit_y = model(lm_result.params.valuesdict(), fit_x)
                
                fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})
                
                x_span = abs(np.max(line.x) - np.min(line.x))
                xe = 0.5
                fit_x = np.linspace(np.min(line.x)-xe*x_span, np.max(line.x)+xe*x_span, num=2000)
                fit_y = model(lm_result.params.valuesdict(), fit_x)
                fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':2.0})
                    
                return lm_result, fit_line, fit_line_extended
                
        
        # Histogram of areas
        y, x = np.histogram(particle_sizes, bins=run_args['hist_bins'], range=[0, max(particle_sizes)*1.05])
        
        # Instead of having x be ranges for each bar, center the x on the average of each range
        xc = [ (x[i]+x[i+1])/2.0 for i in range(len(x)-1) ]
        x = x[:-1]
        
        hist = DataHistogram_current(x=x, y=y, x_label='Area', x_rlabel='$A \, (\mathrm{\mu m}^{2})$', y_label='count')
        hist.mean = results['area_average']
        hist.std = results['area_std']
        hist.median = results['area_median']
        
        outfile = self.get_outfile('particle_areas', output_dir, ext='.png')
        hist.plot(save=outfile, plot_range=[0, None, 0, None], plot_buffers=[0.15,0.05,0.18,0.05],)
        outfile = self.get_outfile('particle_areas', output_dir, ext='.dat')
        np.savetxt(outfile, np.transpose([x, y]))
        
        
        
        # Histogram of radii
        y, x = np.histogram(particle_radii, bins=run_args['hist_bins'], range=[0, max(particle_radii)*1.05])
        
        # Instead of having x be ranges for each bar, center the x on the average of each range
        xc = [ (x[i]+x[i+1])/2.0 for i in range(len(x)-1) ]
        x = x[:-1]
        
        hist = DataHistogram_current(x=x, y=y, x_label='Radius', x_rlabel='$r \, (\mathrm{\mu m})$', y_label='count')
        hist.mean = results['radius_average']
        hist.std = results['radius_std']
        hist.median = results['radius_median']
        
        outfile = self.get_outfile('particle_radii', output_dir, ext='.png')
        hist.plot(save=outfile, plot_range=[0, None, 0, None], plot_buffers=[0.15,0.05,0.18,0.05],)
        outfile = self.get_outfile('particle_radii', output_dir, ext='.dat')
        np.savetxt(outfile, np.transpose([x, y]))
        
        

    def _particle_color(self, im_rgb, im_hsv, labeled_array, scale2, properties, output_dir, results, **run_args):
        
        print('  Computing particle colors')
        properties['rgb'] = []
        properties['hsv'] = []
        properties['rgb_dominant'] = []
        properties['hsv_dominant'] = []
        properties['hue_dominant'] = []
        properties['dominant_h_only'] = []
        
        i_max = len(properties['label'])
        for i, i_label in enumerate(properties['label']):
            if i%100==0:
                print('    {}/{} = {:.1f}%'.format(i, i_max, 100.*i/i_max))
            
            particle = im_rgb[labeled_array==i_label]
            properties['rgb'].append(particle.mean(axis=0))
            
            rgb_dominant, hsv_dominant = self.dominant_color(particle)
            properties['rgb_dominant'].append(rgb_dominant)
            properties['hsv_dominant'].append(hsv_dominant)

            hsv = mpl.colors.rgb_to_hsv(rgb_dominant/255)
            properties['hue_dominant'].append(hsv[0]) # hue[0..1] (single float)
            rgb = mpl.colors.hsv_to_rgb([hsv[0], 1.0, 0.75])*255
            properties['dominant_h_only'].append(rgb) # RGB[0..255] based on dominant hue
            

            particle = im_hsv[labeled_array==i_label]
            hsv = particle.mean(axis=0)
            
            # Hue requires handling as a circular average
            theta = 2*np.pi*(particle[:,0]/255.)
            x, y = np.cos(theta), np.sin(theta)
            theta_avg = np.arctan2(np.average(y), np.average(x)) # [-pi, +pi]
            if theta_avg<0:
                theta_avg = theta_avg - 2*np.pi
            hue_avg = (theta_avg/(2*np.pi))*255 # [0, 255]
            hsv[0] = hue_avg
            
            properties['hsv'].append(hsv)

        properties['rgb'] = np.asarray(properties['rgb'])
        properties['hsv'] = np.asarray(properties['hsv'])
        properties['rgb_dominant'] = np.asarray(properties['rgb_dominant'])
        properties['hsv_dominant'] = np.asarray(properties['hsv_dominant'])
        properties['hue_dominant'] = np.asarray(properties['hue_dominant'])
        properties['dominant_h_only'] = np.asarray(properties['dominant_h_only'])


        
        return results, properties


    def dominant_color(self, rgb_list, n_colors=4):
        #https://stackoverflow.com/questions/43111029/how-to-find-the-average-colour-of-an-image-in-python-with-opencv
        
        import cv2

        # RGB
        pixels = np.float32(rgb_list.reshape(-1, 3))
        
        criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 200, .1)
        flags = cv2.KMEANS_RANDOM_CENTERS

        _, labels, palette = cv2.kmeans(pixels, n_colors, None, criteria, 10, flags)
        _, counts = np.unique(labels, return_counts=True)
        
        rgb_dominant = palette[np.argmax(counts)]
        
        
        # HSV
        hsv = mpl.colors.rgb_to_hsv(pixels/255)
        theta = 2*np.pi*hsv[:,0]
        x, y = np.cos(theta), np.sin(theta)
        pixels = np.stack( (x, y, hsv[:,1], hsv[:,2]), axis=1)

        _, labels, palette = cv2.kmeans(pixels, n_colors, None, criteria, 10, flags)
        _, counts = np.unique(labels, return_counts=True)
        
        dominant = palette[np.argmax(counts)]
        t = np.arctan2(dominant[1], dominant[0])
        if t<0:
            t = t - 2*np.pi
        hue = t/(2*np.pi)
        hsv_dominant = [hue, dominant[2], dominant[3]]
        
        return rgb_dominant, hsv_dominant
    

    def _particle_shape(self, labeled_array, scale2, properties, output_dir, results, **run_args):
        
        import skimage.measure as measure
        scale = np.sqrt(scale2) # µm/pixel
        
        print('  Computing particle shapes')
        properties['PrA'] = []
        properties['eccentricity'] = []
        properties['eccentricity_approx'] = []
        properties['width_ratio'] = []
        
        i_max = len(properties['label'])
        for i, i_label in enumerate(properties['label']):
            if i%100==0:
                print('    {}/{} = {:.1f}%'.format(i, i_max, 100.*i/i_max))
            
            particle = (labeled_array==i_label).astype(np.uint8)

            perimeter_pix = measure.perimeter(particle)
            perimeter_um = perimeter_pix*scale
            
            #area_pix = np.sum(particle)
            #area_um2 = area_pix*scale2 # µm^2
            #radius_um = np.sqrt(area_um2/np.pi) # µm
            
            #PrA = perimeter_um*radius_um/area_um2
            PrA = perimeter_um*properties['radius_um'][i]/properties['area_um2'][i]
            properties['PrA'].append(PrA)
            
            # Estimate eccentricity
            #http://gisaxs.com/index.php/PrA#Ellipse
            if PrA<=2:
                # It should be impossible to have PrA<2, but the inaccuracies of pixel-wise P and A calcs sometimes give PrA slightly less than 2.0
                d = 2
                e = 0
            else:
                d = PrA**2 - 2 + PrA*np.sqrt(PrA**2-4)
                e = np.sqrt( 1 - 4/(d**2) )
            
            properties['eccentricity_approx'].append(e)
            #properties['width_ratio'].append(0.5*d)
            
            
            
            # Fit the particle to an ellipse
            contour = measure.find_contours(particle, 0.5)[0]
            ellipse = measure.EllipseModel()
            ellipse.estimate(contour)
            
            if ellipse.params is None:
                if run_args['verbosity']>=1:
                    print("WARNING: ellipse.params is None for particle {:d} (label {:d})".format(i, i_label))
                properties['eccentricity'].append(0)
                properties['width_ratio'].append(1)
            else:
                xc, yc, a, b, theta = ellipse.params
                if a>=b:
                    eccentricity = np.sqrt(1 - b**2/a**2)
                else:
                    eccentricity = np.sqrt(1 - a**2/b**2)
                properties['eccentricity'].append(eccentricity)
                properties['width_ratio'].append( 1/(np.sqrt(1-np.square(eccentricity))) )

            
            
            
        properties['PrA'] = np.asarray(properties['PrA'])
        properties['eccentricity'] = np.asarray(properties['eccentricity'])
        properties['eccentricity_approx'] = np.asarray(properties['eccentricity_approx'])
        properties['width_ratio'] = np.asarray(properties['width_ratio'])

        
        return results, properties
    

    
        
    def _plot_image(self, im_rgb, labeled_array, properties, output_dir, **run_args):
        
        h, w, c = im_rgb.shape
        aspect = w/h
        if not isinstance(run_args['dpi'], (int, float)):
            # Try to retain the 'real' size of the image
            run_args['dpi'] = w/10
        self.fig = plt.figure( figsize=(10,10/aspect), facecolor='white' )
        self.ax = self.fig.add_axes( [0, 0, 1, 1] )
        
        im_output = im_rgb.astype(int)
        im_output[labeled_array==0] = [0, 0, 0] # Set background to black


        if run_args['verbosity']>=5:
            # Original data, with thresholded background set to black
            self.ax.imshow(im_output)
            
            self.ax.axis( [0, w, h, 0] )
            if run_args['show']:
                plt.show()
            #outfile = self.get_outfile(data.name, output_dir)
            outfile = self.get_outfile('objects', output_dir, ext='.png', ir=True)
            plt.savefig(outfile, dpi=run_args['dpi']*run_args['resize'])
        
        
        if run_args['verbosity']>=3:
            # Color-code each dot based on its average color
            
            print('  Output image')
            i_max = len(properties['label'])
            for i, i_label in enumerate(properties['label']):
                if i%100==0:
                    print('    {}/{} = {:.1f}%'.format(i, i_max, 100.*i/i_max))
                
                # Naive average (looks quite washed out)
                #rgb = [int(properties['rgb'][i][0]), int(properties['rgb'][i][1]), int(properties['rgb'][i][2])] 
                
                # HSV average (tends to look desaturated)
                #hsv = properties['hsv'][i]
                #hsv = [properties['hsv'][i][0], 1.0, 0.75]
                #rgb = mpl.colors.hsv_to_rgb(hsv)*255
                
                # Dominant colors match more closely what a human sees (especially when using only H channel)
                #rgb = properties['rgb_dominant'][i]
                #hsv = mpl.colors.rgb_to_hsv(rgb/255)
                #rgb = mpl.colors.hsv_to_rgb([hsv[0], 1.0, 0.75])*255

                #hsv = properties['hsv_dominant'][i]
                #rgb = mpl.colors.hsv_to_rgb(hsv)*255
                #rgb = mpl.colors.hsv_to_rgb([hsv[0], 1.0, 0.75])*255
                
                rgb = properties['dominant_h_only'][i]
                
                im_output[labeled_array==i_label] = rgb
                
                
            self.ax.imshow(im_output)
            
            #xi, xf, yi, yf = self.ax.axis()
            self.ax.axis( [0, w, h, 0] )
            if run_args['show']:
                plt.show()
            outfile = self.get_outfile('dominant', output_dir, ext='.png', ir=True)
            plt.savefig(outfile, dpi=run_args['dpi']*run_args['resize'])
        
        plt.close(self.fig.number)

                
    def _plot_feature_space(self, properties, output_dir, **run_args):
        
        hue = properties['hue_dominant']*360
        r = properties['radius_um']
        
        fsd = FeatureSpaceDots(x=hue, y=r, x_label='hue', y_label='size', x_rlabel='$\mathrm{hue}$', y_rlabel='$r \, (\mathrm{\mu m})$')
        fsd.properties = properties
        fsd.Protocol_run_args = run_args
        
        outfile = self.get_outfile('FeatureSpace', output_dir, ext='.png', ir=False)
        plot_range = [-180, 180, 0, np.clip(np.max(r), 2, 10)]
        fsd.plot(save=outfile, plot_range=plot_range)



class FeatureSpaceDots():

    def __init__(self, x=None, y=None, name=None, plot_args=None, **kwargs):
        
        self.x = x
        self.y = y
        self.name = name

        self.x_label = kwargs['x_label'] if 'x_label' in kwargs else 'x'
        self.y_label = kwargs['y_label'] if 'y_label' in kwargs else 'y'
        self.x_rlabel = kwargs['x_rlabel'] if 'x_rlabel' in kwargs else self.x_label
        self.y_rlabel = kwargs['y_rlabel'] if 'y_rlabel' in kwargs else self.y_label
        self.x_err = kwargs['x_err'] if 'x_err' in kwargs else None
        self.y_err = kwargs['y_err'] if 'y_err' in kwargs else None
        
        self.plot_valid_keys = ['color', 'linestyle', 'linewidth', 'marker', 'markerfacecolor', 'markersize', 'alpha', 'markeredgewidth', 'markeredgecolor', 'capsize', 'ecolor', 'elinewidth']
        self.plot_args = { 'color' : 'k',
                        'marker' : 'o',
                        'linewidth' : 3.0,
                        'rcParams': {'axes.labelsize': 35,
                                        'xtick.labelsize': 30,
                                        'ytick.labelsize': 30,
                                        },
                            }        
        if plot_args: self.plot_args.update(plot_args)
        self._kwargs = kwargs # Save incase later methods depend on these settings
        
        self.alpha_main = 0.15
        self.selection_bounds = {}
        
    
    def plot(self, save=None, show=False, plot_range=[None,None,None,None], plot_buffers=[0.18,0.02,0.15,0.04], **kwargs):
        '''Plots the data.
        
        Parameters
        ----------
        save : str
            Set to 'None' to avoid saving to disk. Provide filename to save.
        show : bool
            Set to true to open an interactive window.
        plot_range : [float, float, float, float]
            Set the range of the plotting (None scales automatically instead).
        '''  
        
        self._plot(save=save, show=show, plot_range=plot_range, plot_buffers=plot_buffers, **kwargs)
        
        
    def _plot(self, save=None, show=False, plot_range=[None,None,None,None], plot_buffers=[0.18,0.02,0.15,0.02], error=False, error_band=False, xlog=False, ylog=False, xticks=None, yticks=None, dashes=None, transparent=False, figsize=(10,10), **kwargs):
        
        # DataLine._plot()
        
        plot_args = self.plot_args.copy()
        plot_args.update(kwargs)
        self.process_plot_args(**plot_args)
        
        if not isinstance(figsize, (list,tuple,np.ndarray)):
            figsize = (figsize, figsize)
        self.fig = plt.figure( figsize=figsize, facecolor='white' )
        left_buf, right_buf, bottom_buf, top_buf = plot_buffers
        fig_width = 1.0-right_buf-left_buf
        fig_height = 1.0-top_buf-bottom_buf
        
        
        top_h_buf, right_h_buf = 0.1, 0.1
        hue_buf = 0.015
        
        self.ax = self.fig.add_axes( [left_buf, bottom_buf, fig_width-right_h_buf, fig_height-top_h_buf-hue_buf] )
        #self.ax.set_facecolor('k')
        
        self.ax_hue = self.fig.add_axes( [left_buf, bottom_buf+(fig_height-top_h_buf-hue_buf), fig_width-right_h_buf, hue_buf] )
        self.ax_top_h = self.fig.add_axes( [left_buf, bottom_buf+(fig_height-top_h_buf), fig_width-right_h_buf, top_h_buf] )
        self.ax_right_h = self.fig.add_axes( [left_buf+fig_width-right_h_buf, bottom_buf, right_h_buf, fig_height-top_h_buf-hue_buf] )

        
        
        p_args = dict([(i, plot_args[i]) for i in self.plot_valid_keys if i in plot_args])
        self._plot_main(plot_range=plot_range, error=error, error_band=error_band, dashes=dashes, **p_args)
        
        
        self.ax.set_xlabel(self.x_rlabel) #plt.xlabel(self.x_rlabel)
        self.ax.set_ylabel(self.y_rlabel) #plt.ylabel(self.y_rlabel)
        
        
        if xlog:
            plt.semilogx()
        if ylog:
            plt.semilogy()
        if xticks is not None:
            self.ax.set_xticks(xticks)
        else:
            xticks = [-180, -90, 0, +90, +180]
            self.ax.set_xticks(xticks)
            tick_labels = ['{:.0f}'.format(tick if tick>=0 else tick+360).replace('-', '−') for tick in xticks]
            self.ax.set_xticklabels(tick_labels)
        if yticks is not None:
            self.ax.set_yticks(yticks)


        self.ax.tick_params(direction='out', length=8, pad=6)


        self.ax_hue.get_xaxis().set_visible(False)
        self.ax_hue.get_yaxis().set_visible(False)
        self.ax_top_h.get_xaxis().set_visible(False)
        self.ax_top_h.get_yaxis().set_visible(False)
        self.ax_right_h.get_xaxis().set_visible(False)
        self.ax_right_h.get_yaxis().set_visible(False)



        if 'gridlines' in plot_args and plot_args['gridlines']:
            plt.grid()
        
        if 'title' in plot_args and isinstance(plot_args['title'], str):
            #size = plot_args['rcParams']['axes.labelsize']
            size = plot_args['rcParams']['xtick.labelsize']
            size *= 0.75 # Make text smaller
            plt.figtext(0, 1, plot_args['title'], size=size, weight='bold', verticalalignment='top', horizontalalignment='left')
        
        # Axis scaling
        xi, xf, yi, yf = self.ax.axis()
        if plot_range[0] != None: xi = plot_range[0]
        if plot_range[1] != None: xf = plot_range[1]
        if plot_range[2] != None: yi = plot_range[2]
        if plot_range[3] != None: yf = plot_range[3]
        self.ax.axis( [xi, xf, yi, yf] )
        
        if 'reflines' in plot_args:
            # Plot a series of vertical reference lines at the specified x-values.
            for i, xs in enumerate(plot_args['reflines']):
                color_list = ['purple', 'darkblue', 'blue', 'cyan'] # Use distinct color for first few lines
                color = 'lightblue' if i>=len(color_list) else color_list[i] # Use generic color thereafter
                if not isinstance(xs, (tuple, list, np.ndarray) ):
                    # Each refline can either be a single x-value, or a sequence of x-values that form a series
                    xs = [xs]
                for xpos in xs:
                    self.ax.axvline(xpos, color=color, dashes=[3,3])
                    self.ax.text(xpos, yf, str(xpos), size=12, color=color, verticalalignment='top', horizontalalignment='left', rotation=90)
        
        self._plot_extra(**plot_args)
        
        if save:
            if 'dpi' in plot_args:
                plt.savefig(save, dpi=plot_args['dpi'], transparent=transparent)
            else:
                plt.savefig(save, transparent=transparent)
        
        if show:
            self._plot_interact()
            plt.show()
            
        plt.close(self.fig.number)

    
    def _plot_main(self, plot_range=None, error=False, error_band=False, dashes=None, **plot_args):

        sizes = np.clip( self.properties['radius_um']*(2000/2.0), 2, 2000)
        #colors = self.properties['dominant_h_only']/255
        colors = self.properties['rgb']/255
        
        self.x = np.where(self.x>180, self.x-360, self.x)
        
        if False:
            self.ax.scatter(self.x, self.y, s=sizes, c=colors, alpha=self.alpha_main)
            
        else:
            from matplotlib.patches import Ellipse
            xi, xf, yi, yf = plot_range
            xspan, yspan = abs(xf-xi), abs(yf-yi)
            sizing = 0.03
            for i, (xc, yc) in enumerate(zip(self.x, self.y)):
                #ec = self.properties['eccentricity'][i]
                #wr = 1/(np.sqrt(1-ec**2))
                wr = self.properties['width_ratio'][i]
                
                r = yc # um
                A = np.pi*(r**2)
                a = np.sqrt(A*wr/np.pi)
                b = a/wr
                w = xspan*sizing*2*a
                h = yspan*sizing*2*b
                e = Ellipse( (xc, yc), w, h, color=colors[i], alpha=self.alpha_main)
                self.ax.add_patch(e)


    def _plot_extra(self, **plot_args):
        '''This internal function can be over-ridden in order to force additional
        plotting behavior.'''
        
        ra = self.Protocol_run_args
        xi, xf, yi, yf = self.ax.axis()
        

        # Bounds used in data analysis
        size, lw, color = 12, 1, 'b'
        if 'radius_min' in ra:
            lim = ra['radius_min']
            self.ax.axhline(lim, color=color, alpha=0.1, linewidth=lw)
            s = '$r>{:.2f} \, \mathrm{{\mu m}}$'.format(lim)
            self.ax.text(xf, lim, s, size=size, color=color, alpha=0.5, verticalalignment='bottom', horizontalalignment='right')
        if 'area_min' in ra:
            lim = np.sqrt(ra['area_min'])
            self.ax.axhline(lim, color=color, alpha=0.1, linewidth=lw)
            s = '$A>{:.2f} \, \mathrm{{\mu m}}^2$'.format(np.square(lim))
            self.ax.text(xf, lim, s, size=size, color=color, alpha=0.5, verticalalignment='bottom', horizontalalignment='right')
        if 'radius_max' in ra:
            lim = ra['radius_max']
            self.ax.axhline(lim, color=color, alpha=0.1, linewidth=lw)
            s = '$r<{:.2f} \, \mathrm{{\mu m}}$'.format(lim)
            self.ax.text(xf, lim, s, size=size, color=color, alpha=0.5, verticalalignment='top', horizontalalignment='right')
        if 'area_max' in ra:
            lim = np.sqrt(ra['area_max'])
            self.ax.axhline(lim, color=color, alpha=0.1, linewidth=lw)
            s = '$A<{:.2f} \, \mathrm{{\mu m}}^2$'.format(np.square(lim))
            self.ax.text(xf, lim, s, size=size, color=color, alpha=0.5, verticalalignment='top', horizontalalignment='right')
            
            
        # The hue reference
        cmap = mpl.cm.get_cmap('hsv')
        gradient = np.linspace(0, 1, 256*4)
        shift = int( -1*(xi/(xf-xi))*len(gradient) )
        gradient = np.roll(gradient, shift)
        gradient = np.vstack((gradient, gradient))
        self.ax_hue.imshow(gradient, aspect='auto', cmap=cmap)


        # Top histogram
        # TODO: Apply a hard limit to the histogram ranges, so that the bins are always equally spaced
        nbins = 200
        hvals, hbins, hpatches = self.ax_top_h.hist(self.x, bins=nbins, range=[xi, xf])
        self.ax_top_h.axis([xi, xf, 0, np.max(hvals)*0.75])
        
        xspan = abs(xf-xi)
        for x, patch in zip(hbins, hpatches):
            if x<0:
                x = x+360 
            color = cmap(x/xspan)
            patch.set_facecolor(color)
            
            
        # Right histogram
        nbins = 200
        hvals, hbins, hpatches = self.ax_right_h.hist(self.y, bins=nbins, range=[yi, yf], color='0.5', orientation='horizontal')
        self.ax_right_h.axis([0, np.max(hvals)*0.75, yi, yf])

        from matplotlib.patches import Rectangle
        
        # Selection boxes
        x_name, y_name = 'hue_dominant', 'radius_um'
        if self.selection_bounds is not None and len(self.selection_bounds)>0:
            for metric_name, bounds in self.selection_bounds.items():
                xi, xf = bounds[x_name]
                yi, yf = bounds[y_name]
                
                if metric_name=='small':
                    color = 'g'
                elif metric_name=='med':
                    color = 'gold'
                elif metric_name=='big':
                    color = 'r'
                else:
                    color = '0.5'
                
                rect = Rectangle((xi, yi), xf-xi, yf-yi, fill=False, color=color, linewidth=2, linestyle=(0, (3, 6)))
                self.ax.add_patch(rect)
                
                self.ax.text((xi+xf)*0.5, yf, metric_name, size=14, verticalalignment='bottom', horizontalalignment='center', color=color)
                
    

    def process_plot_args(self, **plot_args):
        
        if 'rcParams' in plot_args:
            for param, value in plot_args['rcParams'].items():
                plt.rcParams[param] = value    
