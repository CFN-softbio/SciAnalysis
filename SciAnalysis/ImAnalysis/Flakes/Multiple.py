#!/usr/bin/python
# -*- coding: utf-8 -*-
# vi: ts=4 sw=4

from ..Protocols import *



class average_image(ProtocolMultiple):
    
    def __init__(self, name='average_image', **kwargs):
        
        self.name = self.__class__.__name__ if name is None else name
        
        self.default_ext = '.tiff'
        self.run_args = {
                        'file_extension' : '.tif',
                        'force' : False,
                        'verbosity' : 3,
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
        
        outfile = self.get_outfile(basename, output_dir, ext=run_args['file_extension'])
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
                        'file_extension' : '.jpg',
                        'filename_re' : 'tile_x(\d\d\d)_y(\d\d\d)',
                        'dpi' : 200,
                        'spacing_x' : +1.0,
                        'spacing_y' : -1.0,
                        'overlap' : 0.0,
                        }
        self.run_args.update(kwargs)
    
    
    def _preliminaries(self, datas, **run_args):
        
        # Single image size
        try:
            h, w = datas[0].data.shape
        except:
            img = plt.imread(datas[0].infile) # Deferred load
            h, w, c = img.shape
        if run_args['verbosity']>=4:
            print('  Each image: {}×{} = {:,d} pixels'.format(w, h, w*h))
        
        # Determine total image size
        import re
        filename_re = re.compile(run_args['filename_re'])
        
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
                    
                if 'image_contrast' is not None:
                    in_range = ( run_args['image_contrast'][0]*255, run_args['image_contrast'][1]*255 )
                    import skimage
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
                        'filename_re' : 'tile_x(\d\d\d)_y(\d\d\d)',
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
        
    
    

    
