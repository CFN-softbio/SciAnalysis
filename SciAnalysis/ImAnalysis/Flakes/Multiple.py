#!/usr/bin/python
# -*- coding: utf-8 -*-
# vi: ts=4 sw=4

import pickle
from ..Protocols import *

import matplotlib.patches as patches
import skimage



def re_name_convention(name_convention=None, **kwargs):
    # Naming convention for raw data files, returned as
    # a string suitable for use in RE (regular expressions).
    
    if 'filename_re' in kwargs and kwargs['filename_re'] is not None:
        return kwargs['filename_re']
    
    if name_convention is None or name_convention=='indexed':
        filename_re = '^.+_x(\d+)_y(\d+).*$' # Default (Indexed)
        
    elif name_convention=='ixiy':
        filename_re = '^.+_ix(-?\d+)_iy(-?\d+).*$'
        
    elif name_convention=='xy':
        filename_re = '^.+_x(-?\d+\.\d+)_y(-?\d+\.\d+).*$'
        
    else:
        # If none of the above conditions are met, we simply return the 
        # name_convention, which allows for a custom definition.
        filename_re = name_convention
            
    return filename_re



class mean_image(ProtocolMultiple):
    
    def __init__(self, name='mean_image', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.tiff'
        self.run_args = {
                        'file_extension' : '.tif',
                        'force' : False,
                        'verbosity' : 3,
                        'outname' : None,
                        }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, datas, output_dir, basename, **run_args):
        
        results = {}


        average = None
        for i, data in enumerate(datas):
            
            if run_args['verbosity']>=5:
                print('  Adding image {} ({}); {:.1f}%'.format(i, data.name, 100.*i/len(datas)))
            
            try:
                data_rgb = data.data_rgb
            except:
                data_rgb = plt.imread(data.infile) # Deferred load
                
            
            if average is None:
                average = data_rgb.astype(int)
            else:
                average += data_rgb
                
                
        average = average/len(datas)
        average = average.astype(np.uint8)
        
        if run_args['outname'] is None:
            run_args['outname'] = basename
        outfile = self.get_outfile(run_args['outname'], output_dir, ext=run_args['file_extension'])
        plt.imsave(outfile, average)
        
        return results
        


class average_image(ProtocolMultiple):
    
    def __init__(self, name='average_image', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.tiff'
        self.run_args = {
                        'file_extension' : '.tif',
                        'force' : False,
                        'verbosity' : 3,
                        'average_type' : 'median',
                        'outname' : None,
                        'average_type_in_filename' : True,
                        'proportiontocut' : 0.1,
                        }
        self.run_args.update(kwargs)
        

    @run_default
    def run(self, datas, output_dir, basename, **run_args):
        
        results = {}

        l = len(datas)
        try:
            h, w, c = datas[0].data_rgb.shape
        except:
            data_rgb = plt.imread(datas[0].infile)
            h, w, c = data_rgb.shape
        
        full_dataset = np.zeros( (l, h, w, c) )

        for i, data in enumerate(datas):
            
            if run_args['verbosity']>=5:
                print('  Loading image {} ({}); {:.1f}%'.format(i, data.name, 100.*i/len(datas)))
            
            try:
                data_rgb = data.data_rgb
            except:
                data_rgb = plt.imread(data.infile) # Deferred load
                
            full_dataset[i,:,:] = data_rgb
                

        if run_args['verbosity']>=4:
            print('Computing average of type "{}"'.format(run_args['average_type']))
            
                
        if run_args['average_type']=='mean':
            average = np.mean(full_dataset, axis=0).astype(np.uint8)
        elif run_args['average_type']=='median':
            average = np.median(full_dataset, axis=0).astype(np.uint8)
        elif run_args['average_type']=='mode':
            from scipy.stats import mode
            modes, counts = mode(full_dataset, axis=0)
            average = modes[0].astype(np.uint8)
        elif run_args['average_type']=='trim_mean':
            from scipy.stats import trim_mean
            average = trim_mean(full_dataset, run_args['proportiontocut'], axis=0).astype(np.uint8)
        else:
            print('ERROR: average_type "{}" not recognized.'.format(run_args['average_type']))

        if run_args['verbosity']>=5:
            avg = np.mean(average, axis=(0,1))
            print('  Image average rgb = {}'.format(avg))

        if run_args['outname'] is None:
            run_args['outname'] = basename
        if run_args['average_type_in_filename']:
            run_args['outname'] = '{}_{}'.format(run_args['outname'], run_args['average_type'])
        outfile = self.get_outfile(run_args['outname'], output_dir, ext=run_args['file_extension'])
        plt.imsave(outfile, average)
        
        return results


        



class tile_img(ProtocolMultiple):
    
    def __init__(self, name='tile_img', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.jpg'
        self.run_args = {
                        'blur' : None,
                        'verbosity' : 3,
                        'image_contrast' : (0, 1),
                        'file_extension' : '.png',
                        'filename_re' : None,
                        'name_convention' : None,
                        'dpi' : 200,
                        'spacing_x' : +1.0,
                        'spacing_y' : +1.0,
                        'overlap' : 0.0,
                        }
        self.run_args.update(kwargs)
    
    
    def _preliminaries(self, datas, **run_args):
        
        # Single image size
        try:
            h, w = datas[0].data.shape
            #h, w, c = datas[0].data_rgb.shape
        except:
            img = plt.imread(datas[0].infile) # Deferred load
            h, w, c = img.shape
        if run_args['verbosity']>=4:
            print('  Each image: {}×{} = {:,d} pixels'.format(w, h, w*h))
        
        # Determine total image size
        import re
        filename_re = re_name_convention(**run_args)
        filename_re = re.compile(filename_re)
        
        nrows = 0
        ncols = 0
        for data in datas:
            m = filename_re.match(data.name)
            if m:
                nrows = max(nrows, int(m.groups()[1]) )
                ncols = max(ncols, int(m.groups()[0]) )
                    
            else:
                if run_args['verbosity']>=1:
                    print("ERROR: filename_re didn't match for: {}".format(data.name))


        aspect = (w*ncols)/(h*nrows)
        # direction x == image width == columns
        # direction y == image height == rows
        spacing_x = w*(1-run_args['overlap'])*run_args['spacing_x']
        spacing_y = h*(1-run_args['overlap'])*run_args['spacing_y']
        

        
        if run_args['verbosity']>=4:
            print('  Grid: {}×{} = {} grid positions; {} images total ({:.1f}% coverage)'.format(ncols, nrows, ncols*nrows, len(datas), 100.*len(datas)/(ncols*nrows)))
            print('  Aspect ratio: {:.2f}'.format(aspect))
            
        return h, w, aspect, ncols, nrows, spacing_x, spacing_y, filename_re
    
                
    @run_default
    def run(self, datas, output_dir, basename, **run_args):
        
        results = {}
        
        h, w, aspect, ncols, nrows, spacing_x, spacing_y, filename_re = self._preliminaries(datas, **run_args)
        results['nrows'] = nrows
        results['ncols'] = ncols
        results['aspect'] = aspect
        results['num_images'] = len(datas)
        

        
        self.fig = plt.figure( figsize=(10,10/aspect), facecolor='white' )
        self.ax = self.fig.add_axes( [0, 0, 1, 1] )
        xi, xf, yi, yf = 0, 0, 0, 0
        
        for i, data in enumerate(datas):
            if run_args['verbosity']>=5:
                print('  Adding image {} ({}); {:.1f}%'.format(i, data.name, 100.*i/len(datas)))
                
            m = filename_re.match(data.name)
            if m:
                col = int(m.groups()[0])
                row = int(m.groups()[1])
                
                try:
                    data_rgb = data.data_rgb
                except:
                    data_rgb = plt.imread(data.infile) # Deferred load
                    
                    
                # WARNING: Kludge. Certain PNG files are loaded as float arrays (values from 0 to 1) instead of integers (0 to 255)
                if data_rgb.dtype==np.float32:
                    data_rgb *= 255
                    
                if run_args['image_contrast'] is not None:
                    in_range = ( run_args['image_contrast'][0]*255, run_args['image_contrast'][1]*255 )
                    import skimage.exposure
                    data_rgb = skimage.exposure.rescale_intensity(data_rgb, in_range=in_range, out_range='dtype')
                    
                left = (col-1)*spacing_x
                right = left+w
                bottom = (row-1)*spacing_y
                top = bottom+h
                extent = (left, right, bottom, top)
                self.ax.imshow(data_rgb, extent=extent)
                
                xi = min(xi, left)
                xf = max(xf, right)
                yi = min(yi, bottom)
                yf = max(yf, top)
                
                    
            else:
                if run_args['verbosity']>=1:
                    print("ERROR: filename_re didn't match for: {}".format(data.name))


        #xi, xf, yi, yf = self.ax.axis()
        self.ax.axis( [xi, xf, yi, yf] )
        results['xi'] = xi
        results['xf'] = xf
        results['yi'] = yi
        results['yf'] = yf


        outfile = self.get_outfile(basename, output_dir, ext=run_args['file_extension'])
        results['files_saved'] = [
            { 'filename': '{}'.format(outfile) ,
             'description' : 'combined image of tile scan' ,
             'type' : 'image'
            } ,
            ]
            
        
        plt.savefig(outfile, dpi=run_args['dpi'])
        #plt.show()
        plt.close(self.fig.number)
        
        return results    
    
    
    
    

class tile_svg(tile_img):
    
    
    def __init__(self, name='tile_svg', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.svg'
        self.run_args = {
                        'blur' : None,
                        'verbosity' : 3,
                        'subdir' : None,
                        'subdir_ext' : '.jpg',
                        'file_extension' : '.svg',
                        'filename_re' : None,
                        'name_convention' : None,
                        'spacing_x' : +1.0,
                        'spacing_y' : +1.0,
                        'overlap' : 0.0,
                        'canvas_scale' : 0.1,
                        }
        self.run_args.update(kwargs)
    
    
    @run_default
    def run(self, datas, output_dir, basename, **run_args):
        
        results = {}
        
        h, w, aspect, ncols, nrows, spacing_x, spacing_y, filename_re = self._preliminaries(datas, **run_args)
        scale = run_args['canvas_scale']

        results['nrows'] = nrows
        results['ncols'] = ncols
        results['aspect'] = aspect
        results['num_images'] = len(datas)
        
        
        self.layer_names = ['images', 'labels']
        self.im_id = 4000
        self.txt_id = 4500
        self.layers = ['', '']
        
        
        # Image grid
        for i, data in enumerate(datas):
            if run_args['verbosity']>=5:
                print('  Adding image {} ({}); {:.1f}%'.format(i, data.name, 100.*i/len(datas)))
                
            m = filename_re.match(data.name)
            if m:
                col = int(m.groups()[0])
                row = int(m.groups()[1])
                
                left = (col-1)*spacing_x
                bottom = (row-1)*spacing_y
                
                if run_args['subdir'] is None:
                    path = ''.join( ('../', data.infile) )
                else:
                    path = ''.join( (run_args['subdir'], data.name, run_args['subdir_ext']) )
                self._svg_add_image(path, left*scale, bottom*scale, width=w*scale, aspect=aspect)
                
                    
            else:
                if run_args['verbosity']>=1:
                    print("ERROR: filename_re didn't match for: {}".format(data.name))


        # Row and column headers
        for ic in range(ncols):
            left = (ic)*spacing_x
            self._svg_add_text('x{:03d}'.format(ic+1), left*scale, -spacing_y*scale, layer=1, fontsize=w*scale*0.2)
        for ir in range(nrows):
            bottom = (ir)*spacing_y
            self._svg_add_text('y{:03d}'.format(ir+1), -spacing_x*scale, bottom*scale, layer=1, fontsize=w*scale*0.2)
            



        
        outfile = self.get_outfile(basename, output_dir, ext=run_args['file_extension'])
        results['files_saved'] = [
            { 'filename': '{}'.format(outfile) ,
             'description' : 'SVG layout file of tile scan' ,
             'type' : 'image'
            } ,
            ]
        
        svg = ''.join( (self._svg_header(), self._svg_layers(), self._svg_footer()) )
        with open(outfile, 'w') as fout:
            fout.write(svg)
            
        
        return results    
    
    
    
    
    def _svg_header(self):
        svg = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<!-- Created with Inkscape (http://www.inkscape.org/) -->

<svg
   xmlns:dc="http://purl.org/dc/elements/1.1/"
   xmlns:cc="http://creativecommons.org/ns#"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
   xmlns:svg="http://www.w3.org/2000/svg"
   xmlns="http://www.w3.org/2000/svg"
   xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
   xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
   width="210mm"
   height="297mm"
   viewBox="0 0 210 297"
   version="1.1"
   id="svg832"
   inkscape:version="0.92.3 (2405546, 2018-03-11)"
   sodipodi:docname="layout.svg">
  <defs
     id="defs826" />
  <sodipodi:namedview
     id="base"
     pagecolor="#ffffff"
     bordercolor="#666666"
     borderopacity="1.0"
     inkscape:pageopacity="0.0"
     inkscape:pageshadow="2"
     inkscape:zoom="0.24748738"
     inkscape:cx="-314.36137"
     inkscape:cy="136.71534"
     inkscape:document-units="mm"
     inkscape:current-layer="layer2"
     showgrid="false"
     inkscape:window-width="2560"
     inkscape:window-height="1576"
     inkscape:window-x="0"
     inkscape:window-y="1200"
     inkscape:window-maximized="1" />
  <metadata
     id="metadata829">
    <rdf:RDF>
      <cc:Work
         rdf:about="">
        <dc:format>image/svg+xml</dc:format>
        <dc:type
           rdf:resource="http://purl.org/dc/dcmitype/StillImage" />
        <dc:title />
      </cc:Work>
    </rdf:RDF>
  </metadata>"""
        return svg
    
    def _svg_layers(self):
        svg = ''
        for i, (name, layer_svg) in enumerate(zip(self.layer_names, self.layers)):
            svg += """  <g
     inkscape:label="{}"
     inkscape:groupmode="layer"
     id="layer{}" >\n""".format(name, i+1)
     
            svg += layer_svg
            
            svg += """
  </g>\n"""
     
        return svg
    
    def _svg_footer(self):
        svg = """</svg>"""
        return svg
    
    
    def _svg_add_image(self, filename, xc, yc, width=200, layer=0, aspect=None):
        
        if aspect is None:
            # Determine actual image size (aspect ratio)
            im = ndimage.imread(filename)
            h, w, d = im.shape
            aspect = float(w)/float(h)
        height = width/aspect
        
        svg = """    <image
       width="{:.4f}"
       height="{:.4f}"
       xlink:href="{}"
       id="image{:d}"
       x="{:.4f}"
       y="{:.4f}" />
""".format( width, height, filename, self.im_id, xc-width/2, yc-height/2 )
        
        self.im_id += 1
        self.layers[layer] += svg    
    
    def _svg_add_text(self, text, xc, yc, layer=0, fontsize=30):
        
        svg = """    <text
       sodipodi:linespacing="100%"
       id="text{:d}"
       y="{:.4f}"
       x="{:.4f}"
       style="font-size:{:.6f}px;font-style:normal;font-weight:normal;line-height:100%;letter-spacing:0px;word-spacing:0px;fill:#000000;fill-opacity:1;stroke:none;font-family:Bitstream Vera Sans"
       xml:space="preserve"><tspan
         style="text-align:center;line-height:100%;writing-mode:lr-tb;text-anchor:middle"
         y="{:.4f}"
         x="{:.4f}"
         id="tspan{:d}"
         sodipodi:role="line">{}</tspan></text>
""".format(self.txt_id, yc, xc, fontsize, yc, xc, self.txt_id+1, text)
        
        self.txt_id += 2
        self.layers[layer] += svg
        
    
    
    
    
    


class Data2D_flake_histogram(Data2D):
    
    def _xy_axes(self):
        #dim_y,dim_x = self.data.shape
        
        return self.x_axis, self.y_axis

    def _plot(self, save=None, show=False, ztrim=[0.01, 0.01], size=10.0, plot_buffers=[0.1,0.1,0.1,0.1], xlog=False, ylog=False, side_histograms=None, **kwargs):
        
        # Data2D._plot()
        
        plot_args = self.plot_args.copy()
        plot_args.update(kwargs)
        self.process_plot_args(**plot_args)
        
        
        self.fig = plt.figure( figsize=(size,size), facecolor='white' )
        left_buf, right_buf, bottom_buf, top_buf = plot_buffers
        if 'colorbar' in plot_args and plot_args['colorbar']:
            right_buf += 0.03
        fig_width = 1.0-right_buf-left_buf
        fig_height = 1.0-top_buf-bottom_buf
        
        
        if side_histograms is not None:
            sh = side_histograms
            self.ax = self.fig.add_axes( [left_buf, bottom_buf, fig_width-sh, fig_height-sh] )
            
            # Extra graphs
            self.ax_top = self.fig.add_axes( [left_buf, bottom_buf+fig_height-sh, fig_width-sh, sh] )
            self.ax_right = self.fig.add_axes( [left_buf+fig_width-sh, bottom_buf, sh, fig_height-sh] )
            
        else:
            sh = 0
            self.ax = self.fig.add_axes( [left_buf, bottom_buf, fig_width, fig_height] )
        
        plt.sca(self.ax)
        
        
        # Set zmin and zmax. Top priority is given to a kwarg to this plot function.
        # If that is not set, the value set for this object is used. If neither are
        # specified, a value is auto-selected using ztrim.
        
        values = np.sort( self.data.flatten() )
        if 'zmin' in plot_args and plot_args['zmin'] is not None:
            zmin = plot_args['zmin']
        elif self.z_display[0] is not None:
            zmin = self.z_display[0]
        else:
            zmin = values[ +int( len(values)*ztrim[0] ) ]
            
        if 'zmax' in plot_args and plot_args['zmax'] is not None:
            zmax = plot_args['zmax']
        elif self.z_display[1] is not None:
            zmax = self.z_display[1]
        else:
            idx = -int( len(values)*ztrim[1] )
            if idx>=0:
                idx = -1
            zmax = values[idx]
            
        if zmax==zmin:
            zmax = max(values)
            
        print( '        data: %.2f to %.2f\n        z-scaling: %.2f to %.2f\n' % (np.min(self.data), np.max(self.data), zmin, zmax) )
        
        self.z_display[0] = zmin
        self.z_display[1] = zmax
        self._plot_z_transform()
            
        
        shading = 'flat'
        #shading = 'gouraud'
        
        if 'cmap' in plot_args:
            cmap = plot_args['cmap']
            
        else:
            # http://matplotlib.org/examples/color/colormaps_reference.html
            #cmap = mpl.cm.RdBu
            #cmap = mpl.cm.RdBu_r
            #cmap = mpl.cm.hot
            #cmap = mpl.cm.gist_heat
            cmap = mpl.cm.jet
         
        x_axis, y_axis = self.xy_axes()
        extent = [x_axis[0], x_axis[-1], y_axis[0], y_axis[-1]]

        # NOTE: This will not correctly handle if z_display is not set to 'linear'
        
        #self.im = plt.imshow(self.Z, vmin=0, vmax=1, cmap=cmap, interpolation='nearest', extent=extent, origin='lower')
        mesh = self.ax.pcolormesh( self.x_axis, self.y_axis, self.Z, cmap=cmap, vmin=zmin, vmax=zmax, shading=shading )
        
        if xlog:
            self.ax.semilogx()
        if ylog:
            self.ax.semilogy()
        
        
        if self.regions is not None:
            for region in self.regions:
                plt.imshow(region, cmap=mpl.cm.spring, interpolation='nearest', alpha=0.75)
                #plt.imshow(np.flipud(region), cmap=mpl.cm.spring, interpolation='nearest', alpha=0.75, origin='lower')

        x_label = self.x_rlabel if self.x_rlabel is not None else self.x_label
        y_label = self.y_rlabel if self.y_rlabel is not None else self.y_label
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        
        if 'xticks' in kwargs and kwargs['xticks'] is not None:
            self.ax.set_xticks(kwargs['xticks'])
        if 'yticks' in kwargs and kwargs['yticks'] is not None:
            self.ax.set_yticks(kwargs['yticks'])
        
        if 'colorbar' in plot_args and plot_args['colorbar']:
            cbaxes = self.fig.add_axes([left_buf+fig_width+0.01, bottom_buf, 0.02, fig_height-sh]) 
            cbaxes.tick_params(labelsize=10, width=0.5)
            
            if self.z_display[2]=='linear':
                cbar = plt.colorbar(mesh, ax=self.ax, cax=cbaxes, fraction=0.04, pad=0.03, aspect=30)
            else:
                # Handling this correctly would require vmin=0, vmax=1
                pass
                #n = 4
                #colorbar_labels = [ zmin + i*(zmax-zmin)/n for i in range(n+1) ]
                #tick_positions = self._plot_z_transform(data=colorbar_labels, set_Z=False)
                #cbar = plt.colorbar(mesh, ax=self.ax, cax=cbaxes, ticks=tick_positions, fraction=0.04, pad=0.03, aspect=30)
                #colorbar_labels = ["{:.2g}".format(c) for c in colorbar_labels]
                #cbar.ax.set_yticklabels(colorbar_labels, size=18)
        
                
        
        if 'plot_range' in plot_args:
            plot_range = plot_args['plot_range']
            # Axis scaling
            xi, xf, yi, yf = self.ax.axis()
            if plot_range[0] != None: xi = plot_range[0]
            if plot_range[1] != None: xf = plot_range[1]
            if plot_range[2] != None: yi = plot_range[2]
            if plot_range[3] != None: yf = plot_range[3]
            self.ax.axis( [xi, xf, yi, yf] )
        
        if 'title' in plot_args:
            #size = plot_args['rcParams']['axes.labelsize']
            size = plot_args['rcParams']['xtick.labelsize']
            plt.figtext(0, 1, plot_args['title'], size=size, weight='bold', verticalalignment='top', horizontalalignment='left')
        
        self._plot_extra(side_histograms=side_histograms, xlog=xlog, ylog=ylog, **plot_args)
        
        if save:
            if 'transparent' not in plot_args:
                plot_args['transparent'] = True
            if 'dpi' in plot_args:
                plt.savefig(save, dpi=plot_args['dpi'], transparent=plot_args['transparent'])
            else:
                plt.savefig(save, transparent=plot_args['transparent'])
        
        if show:
            self._plot_interact()
            plt.show()
            
        plt.close(self.fig.number)



    def _plot_extra(self, xlog=False, ylog=False, side_histograms=None, **plot_args):
                
        xi, xf, yi, yf = self.ax.axis()
        
        target_color = plot_args['target_color']


        size = plot_args['rcParams']['xtick.labelsize']
        s = r'$\rho_{{\mathrm{{overall}} }} = {:.1f} \, \mathrm{{flakes/mm^{{2}} }} \,\, ({:,d} \, \mathrm{{total\,flakes)}}$'.format(self.overall_density, int(self.num_flakes))
        plt.figtext(0, 0, s, size=size, color='0.5', verticalalignment='bottom', horizontalalignment='left')
            
        
        if side_histograms is not None:
            
            # Top plot
            if True:
                
                plt.sca(self.ax_top)
                
                x, y, w = self.bins_to_staircase(self.hist_r, self.bin_edges_r)
                yplotf = max(y)*0.75
                yplotf = 20
                
                self.ax_top.plot(x, y, '-', linewidth=2, color=self.color_e)
                self.ax_top.fill_between(x, y, y2=0, color=self.color_e, alpha=0.5)
                
                if xlog:
                    self.ax_top.semilogx()
                self.ax_top.axis( [xi, xf, 0, yplotf] )
                
                # Ticks off
                if True:
                    self.ax_top.get_xaxis().set_ticks([])
                    self.ax_top.tick_params(axis='x',which='both',bottom='off')
                
                self.ax_top.set_ylabel(r'$\rho \, (\mathrm{mm^{-2}})$', size=15)
                #self.ax_top.get_yaxis().set_ticks([10, 20, 30, 40, 50])
                self.ax_top.get_yaxis().set_ticks(range(5, yplotf+5, 5))
                self.ax_top.tick_params(axis='y', which='major', labelsize=10, width=1)


                lm_result, fit_line, fit_line_extended = self.fit_power(line=DataLine(x=x, y=y), verbosity=3, fit_range=[1e0, 2e1])
                self.ax_top.plot(fit_line_extended.x, fit_line_extended.y, '-', linewidth=1, color=self.color_e)
                
                prefactor = lm_result.params['prefactor'].value
                alpha = lm_result.params['alpha'].value
                s = r'$\rho = {:.0f} r^{{ {:.1f} }}$'.format(prefactor, alpha)
                self.ax_top.text(xf, yplotf, s, size=20, color=self.color_e, horizontalalignment='right', verticalalignment='top')
                
                if plot_args['target']:
                    x, y, w = self.bins_to_staircase(self.hist_rtarget, self.bin_edges_r)
                    self.ax_top.plot(x, y, '-', linewidth=1, color=target_color)
                    idx = np.where( (x>plot_args['target'][0]) & (x<plot_args['target'][1]) )
                    self.ax_top.plot(x[idx], y[idx], '-', linewidth=2, color=target_color)
                    
                    num_target_top = np.sum(y[idx])
                    
                    lm_result, fit_line, fit_line_extended = self.fit_power(line=DataLine(x=x, y=y), verbosity=3, fit_range=[plot_args['target'][0], plot_args['target'][1]])

                    self.ax_top.plot(fit_line_extended.x, fit_line_extended.y, '-', linewidth=1, color=target_color)
                    
                    prefactor = lm_result.params['prefactor'].value
                    alpha = lm_result.params['alpha'].value
                    s = r'$\rho = {:.0f} r^{{ {:.1f} }}$'.format(prefactor, alpha)
                    self.ax_top.text(xf, yplotf*0.75, s, size=20, color=target_color, horizontalalignment='right', verticalalignment='top')
                    
            
            
            # Right plot
            if True:
                plt.sca(self.ax_right)
                
                x, y, w = self.bins_to_staircase(self.hist_v, self.bin_edges_v)
                self.ax_right.plot(y, x, '-', linewidth=2, color=self.color_e)
                self.ax_right.fill_betweenx(x, y, 0, color=self.color_e, alpha=0.5)
                
                
                if ylog:
                    self.ax_right.semilogx()
                self.ax_right.axis( [0, yplotf, yi, yf] )

                self.ax_right.set_xlabel(r'$\rho \, (\mathrm{mm^{-2}})$', size=15)
                self.ax_right.get_xaxis().set_ticks(range(5, yplotf+5, 5))
                self.ax_right.tick_params(axis='x', which='major', labelsize=10, width=1)
                
                self.ax_right.axhline(0, linewidth=1, color='k', alpha=0.5)

                # Ticks off
                if True:
                    self.ax_right.get_yaxis().set_ticks([])
                    self.ax_right.tick_params(axis='y',which='both',bottom='off')

                if plot_args['target']:
                    x, y, w = self.bins_to_staircase(self.hist_vtarget, self.bin_edges_v)
                    self.ax_right.plot(y, x, '-', linewidth=1, color=target_color)
                    idx = np.where( (x>plot_args['target'][2]) & (x<plot_args['target'][3]) )
                    self.ax_right.plot(y[idx], x[idx], '-', linewidth=2, color=target_color)
                    
                    num_target_right = np.sum(y[idx])
                    
        if plot_args['target']:
            xti, xtf, yti, ytf = plot_args['target']
            rect = patches.Rectangle((xti,yti),xtf-xti,ytf-yti,linewidth=3,edgecolor=target_color,facecolor='none')
            self.ax.add_patch(rect)
            
            self.ax.text(xtf, ytf, 'target flakes', size=25, color=target_color, horizontalalignment='right', verticalalignment='bottom')
            s = r'${:.1f} \,  \mathrm{{ flakes/mm^{{2}} }}$'.format( (num_target_top+num_target_right)*0.5 )
            self.ax.text(xtf, yti, s, size=25, color=target_color, horizontalalignment='right', verticalalignment='top')
            

    def fit_power(self, line=None, **run_args):
        
        # Usage: lm_result, fit_line, fit_line_extended = self.fit_peaks(**run_args)
        
        if line is None:
            line = self

        line_full = line
        if 'fit_range' in run_args:
            line = line.sub_range(run_args['fit_range'][0], run_args['fit_range'][1])
        
        import lmfit
        
        def model(v, x):
            '''Gaussians with constant background.'''
            m = v['prefactor']*np.power(x, v['alpha'])
            return m
        
        def func2minimize(params, x, data):
            
            v = params.valuesdict()
            m = model(v, x)
            
            return m - data
        
        params = lmfit.Parameters()
        params.add('prefactor', value=40, min=0, vary=True)
        params.add('alpha', value=-2.0, min=-3.0, max=-0.1, vary=True)
        
        
        
        lm_result = lmfit.minimize(func2minimize, params, args=(line.x, line.y))
        
        if run_args['verbosity']>=5:
            print('Fit results (lmfit):')
            lmfit.report_fit(lm_result.params)
            
        fit_x = line.x
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':4.0})
        
        #fit_x = np.linspace(np.min(line_full.x)*0.5, np.max(line_full.x)*2.0, num=600)
        fit_x = np.logspace( np.log(np.min(line_full.x)), np.log(np.max(line_full.x)*2.0), num=1000)
        fit_y = model(lm_result.params.valuesdict(), fit_x)
        fit_line_extended = DataLine(x=fit_x, y=fit_y, plot_args={'linestyle':'-', 'color':'r', 'marker':None, 'linewidth':3.0})        

        return lm_result, fit_line, fit_line_extended

    def bins_to_centers(self, hist, bin_edges):
        
        x = (bin_edges[1:] + bin_edges[:-1])/2
        y = hist
        w = bin_edges[1:] - bin_edges[:-1]
        
        return x, y, w


    def bins_to_staircase(self, hist, bin_edges):
        
        xn = np.repeat(bin_edges, 2)
        x = xn[1:-1]
        y = np.repeat(hist, 2)
        w = np.repeat(bin_edges[1:] - bin_edges[:-1], 2)
        
        return x, y, w

        
    def _plot_interact(self):
        super()._plot_interact()
        self.fig.canvas.mpl_connect('button_press_event', self._button_press_event )
        self._last_xy = None
        
    def _button_press_event(self, event, verbosity=3):
        
        if verbosity>=5:
            print('button_press_event: {}'.format(event))
        
        if event.inaxes!=self.ax:
            return

        current_plot_limits = self.ax.axis()
        x = event.xdata
        y = event.ydata
        
        if self._last_xy is not None and self._last_xy==(x,y):
            if verbosity>=4:
                print('    Double-click at (x,y) = ({:.3g}, {:.3g})'.format(x,y))
            self.launch_subwindow(x,y, verbosity=verbosity)
            self._last_xy = None
        else:
            if verbosity>=4:
                print('    Click at (x,y) = ({:.3g}, {:.3g})'.format(x,y))
            self._last_xy = (x,y)
                
        
            
            
        self.fig.canvas.draw() # Update
        
        
    def launch_subwindow(self, x, y, span=1, max_matches=8*8, prepend='    launch_subwindow: ', verbosity=3):
        
        if verbosity>=5:
            print('{} (x, y) = ({:.4g}, {:.4g})'.format(prepend, x, y))
            
            
        # Determine the data range associated with the (x,y) click
        idx_xstart = np.where( self.x_axis<x )[0][-1] - span
        idx_xend = np.where( self.x_axis>x )[0][0] + span
        
        idx_xstart = np.clip(idx_xstart, 0, len(self.x_axis) )
        idx_xend = np.clip(idx_xend, 0, len(self.x_axis) )

        idx_ystart = np.where( self.y_axis<y )[0][-1] - span
        idx_yend = np.where( self.y_axis>y )[0][0] + span
        
        idx_ystart = np.clip(idx_ystart, 0, len(self.y_axis) )
        idx_yend = np.clip(idx_yend, 0, len(self.y_axis) )
        
        target_zone = [self.x_axis[idx_xstart], self.x_axis[idx_xend], self.y_axis[idx_ystart], self.y_axis[idx_yend]]

        if verbosity>=5:
            print('{} xidx {:d} to {:d} (x = {:.3g} to {:.3g})'.format(prepend, idx_xstart, idx_xend, target_zone[0], target_zone[1]))
            print('{} yidx {:d} to {:d} (y = {:.3g} to {:.3g})'.format(prepend, idx_ystart, idx_yend, target_zone[2], target_zone[3]))
        

        # Extract data matching the range
        idx = np.where( (self.r_sizes>target_zone[0]) & (self.r_sizes<target_zone[1]) & (self.values>target_zone[2]) & (self.values<target_zone[3]) )
        idx_sort = np.argsort( -1*self.r_sizes[idx] )
        matches = len(self.r_sizes[idx])
        
        if verbosity>=3:
            print('  {} flakes identified ({} maximum matches)'.format(matches, max_matches))
        
        flakes_sorted = np.asarray(self.flakes)[idx][idx_sort]
        
        name_convention = self._kwargs['name_convention'] if 'name_convention' in self._kwargs else None
        g = ImageGrid(name_convention=name_convention)
        icount = 0
        for i, flake in enumerate(flakes_sorted):
            
            if icount<max_matches:
                icount += 1
                if verbosity>=4:
                    xf, yf = g.file_xy(flake['infile'])
                    y, x = flake['center_of_mass']
                    print('            [{:d}, {:d}, {:.0f}, {:.0f}], # r = {:.3f}, v = {:.3f}'.format(xf, yf, x, y, flake['radius_um'], flake['flake_contrast']) )
                
                g.flakes.append(flake)
                
            
            
        if verbosity>=3:
            print('    {} flakes added to ImageGrid'.format(len(g.flakes)))
        
        image_contrast = (0, 1) # Unmodified
        #image_contrast = (0.3, 0.78) # 80/255, 200/255
        #image_contrast = (0.5, 1)
        image_contrast = (0.35, 0.6)
        
        g.scale = self.scale
        g.target_zone = target_zone
        g.plot(save=None, show=True, image_contrast=image_contrast, dpi=300)
                    
        
class ImageGrid(object):
    
    def __init__(self, plot_args=None, **kwargs):
        
        self.flakes = []
        self.axes = []
        self.ax_ims = []
        self.flakes_raw = []


        # Naming convention for raw data files
        filename_re = re_name_convention(**kwargs)
        self.re_files = re.compile(filename_re)
        
        self.plot_args = { 'color' : 'k',
                        'marker' : 'o',
                        'linewidth' : 3.0,
                        'rcParams': {'axes.labelsize': 35,
                                        'xtick.labelsize': 30,
                                        'ytick.labelsize': 30,
                                        },
                        'image_contrast': (0,1),
                        'crop_pad': 0.2,
                            }        
        if plot_args: self.plot_args.update(plot_args)        
    
    
    def file_xy(self, filename):
        
        m = self.re_files.match(filename)
        if m:
            xf = int(m.groups()[0])
            yf = int(m.groups()[1])
        else:
            print('RE fail for filename: {}'.format(filename))
        
        return xf, yf

    def process_plot_args(self, **plot_args):
        
        if 'rcParams' in plot_args:
            for param, value in plot_args['rcParams'].items():
                plt.rcParams[param] = value    
    
    def plot(self, save=None, show=False, size=10, textsize=2.0, linewidth=2.0, transparent=True, verbosity=3, **kwargs):
        
        plot_args = self.plot_args.copy()
        plot_args.update(kwargs)
        self.process_plot_args(**plot_args)
        
        
        self.fig = plt.figure( figsize=(size,size), facecolor='white' )


        num = len(self.flakes)
        ncols = np.ceil(np.sqrt(num))
        nrows = np.ceil(num/ncols)
        in_range = ( plot_args['image_contrast'][0]*255, plot_args['image_contrast'][1]*255 )
        
        for i, flake in enumerate(self.flakes):
            
            filename = flake['infile'].replace('\\', '/') # String replace in case files were saved on another platform.
            box = flake ['bbox']
            center = flake['center_of_mass']
            radius = flake['radius_um']
            value = flake['flake_contrast']
            
            if verbosity>=6:
                print('ImageGrid.plot for flake {}, infile: {}'.format(i, filename))
                d = os.path.dirname(os.path.realpath(filename))
                n = len([name for name in os.listdir(d) if '.tif' in name])
                if os.path.exists(filename):
                    print('    exists in {} ({} files)'.format(d, n))
                else:
                    print('    not found in {} ({} files)'.format(d, n))
                print('        cwd: {}'.format(os.getcwd()))
                #print('        cwd: {}'.format(os.path.realpath(os.getcwd())))
                #print('        radius_um = {:.1f}; center = [{:.1f}, {:.1f}]'.format(radius, center[0], center[1]))
            
            
            # Load image
            img = plt.imread(filename)
            h, w, c = img.shape
            y1, y2, x1, x2 = box
            
            # Make the crop border a bit bigger than the flake
            pad = plot_args['crop_pad']*max( abs(x2-x1), abs(y2-y1) )
            x1p = int(np.clip(x1 - pad, 0, w))
            x2p = int(np.clip(x2 + pad, 0, w))
            y1p = int(np.clip(y1 - pad, 0, h))
            y2p = int(np.clip(y2 + pad, 0, h))
            
            # Display flake
            flake = img[y1p:y2p , x1p:x2p, :]
            self.flakes_raw.append(flake)
            flake = skimage.exposure.rescale_intensity(flake, in_range=in_range, out_range='dtype')
            
            
            ax = self.fig.add_subplot(nrows, ncols, i+1)
            self.axes.append(ax)
            self.ax_ims.append( ax.imshow(flake) )

            # Bounding box
            rect = patches.Rectangle( (x1-x1p, y1-y1p), x2-x1, y2-y1, linewidth=0.5*linewidth, edgecolor='orange', facecolor='none', alpha=0.5)
            ax.add_patch(rect)
            ax.text(x2-x1p, y2-y1p, '{:.3f}'.format(value), size=3*textsize, color='orange', horizontalalignment='left', verticalalignment='top', alpha=0.5)
            

            # Size circle
            y, x = center
            xm = x - x1p
            ym = y - y1p
            rpix = radius/self.scale
            circ = patches.Circle(xy=(xm, ym), radius=rpix, linewidth=0.5*linewidth, edgecolor='r', facecolor='none', alpha=0.3)
            ax.add_patch(circ)
            rect = patches.Rectangle( (xm-0.5*rpix, ym), rpix, 0, linewidth=0.5*linewidth, edgecolor='r', facecolor='none', alpha=0.3)
            ax.add_patch(rect)
            rect = patches.Rectangle( (xm, ym-0.5*rpix), 0, rpix, linewidth=0.5*linewidth, edgecolor='r', facecolor='none', alpha=0.3)
            ax.add_patch(rect)
            s = r'${:.1f} \, \mathrm{{\mu m}}$'.format(radius)
            ax.text(xm, ym+rpix, s, size=6*textsize, color='r', horizontalalignment='center', verticalalignment='top', alpha=0.5)
            

            xf, yf = self.file_xy(filename)
                
            ax.text(0, 0, 'x{}, y{}'.format(xf, yf), size=4*textsize, horizontalalignment='left', verticalalignment='bottom')
            ax.text(x2p-x1p, 0, r'$({:.0f}, {:.0f})$'.format(x, y), size=3*textsize, horizontalalignment='right', verticalalignment='bottom')
            ax.set_xticks([])
            ax.set_yticks([])
            ax.axis('off')
            
            if True:
                print('        radius_um = {:.1f}; center = [{:.1f}, {:.1f}]'.format(radius, center[0], center[1]))
                      
            
            
        target_zone = self.target_zone
        s = 'Flakes of size [{:.1f}, {:.1f}] μm and contrast [{:.3f}, {:.3f}]'.format(target_zone[0], target_zone[1], target_zone[2], target_zone[3])
        self.fig.suptitle(s, fontsize=12)
        plt.tight_layout(pad=3.0, w_pad=0.05, h_pad=0.05)

        if save:
            if 'dpi' in plot_args:
                plt.savefig(save, dpi=plot_args['dpi'], transparent=transparent)
            else:
                plt.savefig(save, transparent=transparent)

        if show:
            self._plot_interact()
            plt.show()
            
        plt.close(self.fig.number)



    def _plot_interact(self):
        
        self.fig.canvas.set_window_title('ImageGrid')
        #plt.get_current_fig_manager().toolbar.pan()
        self.fig.canvas.toolbar.pan()
        #self.fig.canvas.mpl_connect('scroll_event', self._scroll_event )
        #self.fig.canvas.mpl_connect('motion_notify_event', self._move_event )
        self.fig.canvas.mpl_connect('key_press_event', self._key_press_event)
        self.fig.canvas.mpl_connect('button_press_event', self._button_press_event)
        
        #self.ax.format_coord = self._format_coord       
        
    def _key_press_event(self, event):
        
        if event.key in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']:
            
            avg = 0
            for im, flake in zip(self.ax_ims, self.flakes_raw):
                avg += np.average(flake)
            avg /= len(self.flakes_raw)*255
            amt = int(event.key)/10.0
            image_contrast = ( avg*amt , 1.0-(1.0-avg)*amt )
            
            in_range = ( image_contrast[0]*255, image_contrast[1]*255 )
            
            for im, flake in zip(self.ax_ims, self.flakes_raw):
                
                flake = skimage.exposure.rescale_intensity(flake, in_range=in_range, out_range='dtype')
                im.set_data(flake)
                
            self.fig.canvas.draw() # Update
            
        
    def _button_press_event(self, event, verbosity=3):
        
        if verbosity>=5:
            print('button_press_event: {}'.format(event))
           
        if event.dblclick:
        
            for ax, flake in zip(self.axes, self.flakes):
                if event.inaxes==ax:
                    if verbosity>=4:
                        print('    Loading flake index {} ({})'.format(flake['index'], flake['name']))
                    from ..Data import Data2DImageRGB
                    from .Protocols import flake_images
                    data = Data2DImageRGB(infile=flake['infile'])
                    f = flake_images()
                    f.run_args['scale'] = flake['scale']
                    f.run_args['scale2'] = flake['scale2']
                    f.run_args['image_contrast_trim'] = 0.4
                    f.plot_flake(data, flake, output_dir='./flake_images', show=True, **f.run_args)
                    return
                

                        
    

    
class histogram(tile_img):
    
    
    def __init__(self, name='histogram', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.png'
        self.run_args = {
                        'scale' : None,
                        'verbosity' : 3,
                        'file_extension' : '.png',
                        'min_area_pixels' : 50,
                        'min_radius_um' : None,
                        'max_radius_um' : 50,
                        'value_range' : (-0.08, 0.08),
                        'value_key' : 'flake_contrast',
                        'target' : [2, 45, -0.015, 0], # 'boxed-in' region for the 'target' (desired) flakes
                        'interact' : False,
                        'zmin' : 0,
                        'ztrim' : [0, 0.005],
                        'colorbar' : True,
                        }
        self.run_args.update(kwargs)
    
    
    def print_structure(self, dictionary, prepend=' '):
        for key, value in dictionary.items():

            if isinstance(value, str):
                print('{}{} : {}'.format(prepend, key, value))
            elif isinstance(value, (int, np.int) ):
                print('{}{} : {:d}'.format(prepend, key, value))
            elif isinstance(value, (float, np.float) ):
                print('{}{} : {:.3f}'.format(prepend, key, value))
            elif hasattr(value, 'shape'):
                print('{}{} : array {}'.format(prepend, key, value.shape))
            elif isinstance(value, (list, tuple)):
                print('{}{} : list len {}'.format(prepend, key, len(value) ))
            else:
                print('{}{} : -'.format(prepend, key))
                
                    
    
    @run_default
    def run(self, datas, output_dir, basename, **run_args):
        
        results = {}
        
        # Aggregate results
        ########################################
        flakes = []
        r_sizes = []
        values = []
        for data in datas:
            with open(data.infile, 'rb') as fin:
                saved = pickle.load(fin) # 'res_map', 'image_labelmap', 'flakes'
                if len(flakes)==0:
                    h, w = saved['res_map'].shape
                
                for flake in saved['flakes']:
                    flakes.append(flake)
                    #self.print_structure(flake)
                    r_sizes.append(flake['radius_um'])
                    values.append(flake[run_args['value_key']])
                    
                if run_args['verbosity']>=5:
                    print('      {} flakes added from image {}'.format(len(saved['flakes']), data.infile))


        r_sizes = np.asarray(r_sizes)
        values = np.asarray(values)
                
        if run_args['verbosity']>=4:
            print('  {:,d} flakes identified in {:d} images'.format(len(flakes), len(datas)))
            print('    radius ranging from {:.2f} μm to {:,.1f} μm (avg = {:.2f} μm, stdev = {:.2f} μm)'.format(min(r_sizes), max(r_sizes), np.average(r_sizes), np.std(r_sizes)))
            print('    values ("{}") ranging from {:.3f} to {:.3f} (avg = {:.3f}, stdev = {:.3f})'.format(run_args['value_key'], min(values), max(values), np.average(values), np.std(values) ))
            
            

        # Set values
        ########################################
        # Set scale for conversion from pixels to um
        if run_args['scale'] is None:
            data = datas[0]
            run_args['scale'] = np.average((data.x_scale, data.y_scale)) # um/pixel
        run_args['scale2'] = run_args['scale']*run_args['scale'] # um^2/pixel^2

        results['num_images'] = len(datas)
        results['num_flakes'] = len(flakes)
        results['total_area_pix'] = results['num_images']*h*w
        results['total_area_um2'] = results['total_area_pix']*run_args['scale2']
        results['total_area_mm2'] = results['total_area_um2']/1e6
        
        results['overall_density'] = results['num_flakes']/results['total_area_mm2']
    
        if run_args['min_radius_um'] is None:
            run_args['min_radius_um'] = np.sqrt(run_args['min_area_pixels']/np.pi)*run_args['scale']
            
            
            
        # Generate histograms
        ########################################
        bin_edges_x = np.logspace(np.log(run_args['min_radius_um']), np.log(run_args['max_radius_um']), num=150, endpoint=True)
        bin_edges_y = np.linspace(run_args['value_range'][0], run_args['value_range'][1], num=200, endpoint=True)
        
        H, xedges, yedges = np.histogram2d(r_sizes, values, bins=[bin_edges_x,bin_edges_y])
        H /= results['total_area_mm2']

        if run_args['verbosity']>=4:
            print('  {:,d} histogram bins'.format(H.size))
            print('    H histogram ranging from {:.3g} to {:.3g} (avg = {:.3g}, stdev = {:.3g})'.format(min(H.flatten()), max(H.flatten()), np.average(H.flatten()), np.std(H.flatten()) ))
        
        bin_edges_rx = np.logspace(np.log(run_args['min_radius_um']), np.log(100), num=400, endpoint=True)
        hist_r, bin_edges_r = np.histogram(r_sizes, bins=bin_edges_rx, range=None, normed=False, weights=None, density=None)
        hist_r = hist_r/results['total_area_mm2']
        
        if run_args['verbosity']>=4:
            print('    hist_r ranging from {:.2f} to {:,.1f} (avg = {:.2f}, stdev = {:.2f})'.format(min(hist_r), max(hist_r), np.average(hist_r), np.std(hist_r) ))
        
        if run_args['target']:
            idx = np.where( ( values>(run_args['target'][2]) ) & ( values<(run_args['target'][3]) ) )
            hist_rtarget, bin_edges_rtarget = np.histogram(r_sizes[idx], bins=bin_edges_rx, range=None, normed=False, weights=None, density=None)
            hist_rtarget = hist_rtarget/results['total_area_mm2']
        
        bin_edges_vy = np.linspace(run_args['value_range'][0], run_args['value_range'][1], num=200, endpoint=True)
        hist_v, bin_edges_v = np.histogram(values, bins=bin_edges_vy, range=None, normed=False, weights=None, density=None)
        hist_v = hist_v/results['total_area_mm2']

        if run_args['verbosity']>=4:
            print('    hist_v ranging from {:.2f} to {:,.1f} (avg = {:.2f}, stdev = {:.2f})'.format(min(hist_v), max(hist_v), np.average(hist_v), np.std(hist_v) ))

        if run_args['target']:
            idx = np.where( ( r_sizes>(run_args['target'][0]) ) & ( r_sizes<(run_args['target'][1]) ) )
            hist_vtarget, bin_edges_v = np.histogram(values[idx], bins=bin_edges_vy, range=None, normed=False, weights=None, density=None)
            hist_vtarget = hist_vtarget/results['total_area_mm2']
        
        
        
        # Generate plot
        ########################################
        d = Data2D_flake_histogram(**run_args)
        d.data = H.transpose()
        d.x_axis = (bin_edges_x[1:] + bin_edges_x[:-1])/2
        d.y_axis = (bin_edges_y[1:] + bin_edges_y[:-1])/2
        
        d.hist_r = hist_r
        d.bin_edges_r = bin_edges_r
        d.hist_v = hist_v
        d.bin_edges_v = bin_edges_v
        
        if run_args['target']:
            d.hist_rtarget = hist_rtarget
            d.hist_vtarget = hist_vtarget
        
        d.color_e = 'darkblue'
        
        
        plot_args = { 'rcParams': {'axes.labelsize': 40,
                                        'xtick.labelsize': 20,
                                        'ytick.labelsize': 20,
                                        'xtick.major.size': 12,
                                        'xtick.major.width': 2,
                                        'xtick.minor.size': 6,
                                        'xtick.minor.width': 1,
                                        'ytick.major.size': 10,
                                        'ytick.major.width': 2,
                                        'ytick.minor.size': 5,
                                        'ytick.minor.width': 1,

                                        },
                            }    
        plot_args['target'] = run_args['target']
        plot_args['target_color'] = 'orange'
        
        d.x_label = 'flake radius'
        d.x_rlabel = r'$\mathrm{flake \,\, radius,} \,\, r \, (\mathrm{\mu m})$'
        d.y_label = 'image contrast'
        d.y_rlabel = r'$\mathrm{image \,\, contrast}$'
        d.plot_args = plot_args
        
        d.num_flakes = results['num_flakes']
        d.overall_density = results['overall_density']
        
        # Optionally massage the data?
        #d.blur(1.0)
        
        # Make connections required for interactive aspects
        d.scale = run_args['scale']
        d.r_sizes = r_sizes
        d.values = values
        d.flakes = flakes
        
        
        outfile = self.get_outfile(basename, output_dir, ext=run_args['file_extension'])
        results['files_saved'] = [
            { 'filename': '{}'.format(outfile) ,
             'description' : '2D histogram of flakes' ,
             'type' : 'image'
            } ,
            ]
            
        plot_range = run_args['min_radius_um'], run_args['max_radius_um'], run_args['value_range'][0], run_args['value_range'][1]
        d.plot(outfile, show=run_args['interact'], xlog=True, ylog=False, plot_buffers=[0.18, 0.05, 0.18, 0.05], plot_range=plot_range, yticks=[-0.08, -0.04, 0, 0.04, 0.08], side_histograms=0.15, dpi=300, plot_args=plot_args, **run_args)
        
    
        return results
    
    
    
