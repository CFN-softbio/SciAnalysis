#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os, sys
import re
import glob
from scipy import ndimage


HEADER ="""<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<!-- Created with Inkscape (http://www.inkscape.org/) -->

<svg
   xmlns:dc="http://purl.org/dc/elements/1.1/"
   xmlns:cc="http://creativecommons.org/ns#"
   xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
   xmlns:svg="http://www.w3.org/2000/svg"
   xmlns="http://www.w3.org/2000/svg"
   xmlns:xlink="http://www.w3.org/1999/xlink"
   xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"
   xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"
   width="744.09448819"
   height="1052.3622047"
   id="svg2"
   version="1.1"
   inkscape:version="0.48.4 r9939"
   sodipodi:docname="Mauter.svg">
  <defs
     id="defs4" />
  <sodipodi:namedview
     id="base"
     pagecolor="#ffffff"
     bordercolor="#666666"
     borderopacity="1.0"
     inkscape:pageopacity="0.0"
     inkscape:pageshadow="2"
     inkscape:zoom="0.12374369"
     inkscape:cx="7481.1052"
     inkscape:cy="1884.3233"
     inkscape:document-units="px"
     inkscape:current-layer="layer1"
     showgrid="false"
     inkscape:window-width="2560"
     inkscape:window-height="1515"
     inkscape:window-x="-4"
     inkscape:window-y="-4"
     inkscape:window-maximized="1" />
  <metadata
     id="metadata7">
    <rdf:RDF>
      <cc:Work
         rdf:about="">
        <dc:format>image/svg+xml</dc:format>
        <dc:type
           rdf:resource="http://purl.org/dc/dcmitype/StillImage" />
        <dc:title></dc:title>
      </cc:Work>
    </rdf:RDF>
  </metadata>
  <g
     inkscape:label="Layer 1"
     inkscape:groupmode="layer"
     id="layer1">
"""

FOOTER = """
  </g>
</svg>
"""

# Filename
################################################################################
class Filename(object):
    '''Parses a filename into pieces following the desired pattern.'''
    
    def __init__(self, filepath):
        '''Creates a new Filename object, to make it easy to separate the filename
        into its pieces (path, file, extension).'''
        
        self.full_filepath = filepath
        self._update()
        
        
    def _update(self):
    
        path, filename, filebase, ext = self.file_split(self.full_filepath)
        
        self.path = path            # filesystem path to file
        self.filename = filename    # filename (including extension)
        self.filebase = filebase    # filename (without extension)
        self.ext = ext              # extension
        

    def file_split(self, filepath):
        
        
        filepath, filename = os.path.split(filepath)
        filebase, ext = os.path.splitext(filename)
        
        return filepath, filename, filebase, ext


    def split(self):
        # filepath, filename, filebase, ext = f.split()
        return self.path, self.filename, self.filebase, self.ext

    def get_filepath(self):
        return self.full_filepath

    def get_path(self):
        return self.path+'/'
        
    def get_filename(self):
        return self.filename
        
    def get_filebase(self):
        return self.filebase
        
    def get_ext(self):
        return self.ext

    def timestamp(self):
        statinfo = os.stat(self.full_filepath)
        filetimestamp = statinfo.st_mtime
        return filetimestamp

    def matches_basename(self, filepath):
        path, filename, filebase, ext = self.file_split(filepath)
        
        return self.filebase==filebase
    
    def append(self, text):
        self.full_filepath = os.path.join(self.path, self.filebase + text + self.ext)
        self._update()

    # End class Filename(object)
    ########################################



class LayoutSVG():
    
    def __init__(self):
        
        self.svg = ''
        self.im_id = 4000
        self.txt_id = 4500
        
        
    def add_image(self, filename, xc, yc, width=200):
        
        # Determine actual image size (aspect ratio)
        im = ndimage.imread(filename)
        h, w, d = im.shape
        aspect = float(w)/float(h)
        height = width/aspect

        svg = """    <image
       width="{:.4f}"
       height="{:.4f}"
       xlink:href="file://{}"
       id="image{:d}"
       x="{:.4f}"
       y="{:.4f}" />
""".format( width, height, filename, self.im_id, xc-width/2, yc-height/2 )
        
        self.im_id += 1
        self.svg += svg
        
        
    def add_text(self, text, xc, yc, fontsize=30):
        
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
        self.svg += svg
        
        
    def add_vtext(self, text, xc, yc, fontsize=30):
        
        svg = u"""    <text
       xml:space="preserve"
       style="font-size:{:.6f}px;font-style:normal;font-weight:normal;line-height:100%;letter-spacing:0px;word-spacing:0px;fill:#000000;fill-opacity:1;stroke:none;font-family:Bitstream Vera Sans"
       x="{:.4f}"
       y="{:.4f}"
       id="text{:d}"
       sodipodi:linespacing="100%"
       transform="matrix(0,-1,1,0,0,0)"><tspan
         sodipodi:role="line"
         id="tspan{:d}"
         x="{:.4f}"
         y="{:.4f}"
         style="text-align:center;line-height:100%;writing-mode:lr-tb;text-anchor:middle">{}</tspan></text>        
""".format(fontsize, yc, xc, self.txt_id, self.txt_id+1, yc, xc, text)

        self.txt_id += 2
        self.svg += svg
        
        
        
    def pattern_image_grid(self, file_re, source_dir='./', col_spacing=200, row_spacing=None, x0=0, y0=0):
        
        if row_spacing is None:
            row_spacing = col_spacing
        
        match_re = re.compile(file_re)
        
        infiles = glob.glob(source_dir+'*')
        infiles.sort()
        
        # Figure out the list of valid files
        filenames = []
        row_strs = []
        col_strs = []
        for infile in infiles:
            
            f = Filename(infile)
            filepath, filename, filebase, ext = f.split()
            
            m = match_re.match(filename)
            if m:
                
                if filename in filenames:
                    print('ERROR: File already added.')
                filenames.append(filename)
                
                col_str = m.groups()[0]
                row_str = m.groups()[1]
                
                if col_str not in col_strs:
                    col_strs.append(col_str)
                if row_str not in row_strs:
                    row_strs.append(row_str)
                    
                    
        # Col headers
        for ic, col_str in enumerate( sorted(col_strs) ):
            
            x = x0 + ic*col_spacing
            y = y0 - float(row_spacing)*( len(row_strs) - 1.0 + 0.5 )
            self.add_text(col_str, x, y, fontsize=14)
        
        # Row headers
        for ir, row_str in enumerate( sorted(row_strs) ):
            
            x = x0 - float(col_spacing)/2
            y = y0 + ir*row_spacing
            
            row_str = 'θ = '+row_str+'°'

            self.add_vtext(row_str, x, y, fontsize=25)
            
            
        # Images
        for ic, col_str in enumerate( sorted(col_strs) ):
            for ir, row_str in enumerate( sorted(row_strs) ):
                
                matching = [s for s in filenames if col_str in s]
                matching = [s for s in matching if row_str in s]
                
                if len(matching)<1:
                    print("  WARNING: No file for col='{}' row='{}'".format(col_str, row_str))
                    
                else:
                
                    filename = matching[0]
                    
                    x = x0 + ic*col_spacing
                    y = y0 - ir*row_spacing
                    print('Image: {}'.format(filename))
                    self.add_image(source_dir+filename, x, y, width=col_spacing)
                
                
        
        
        
    
    def save(self, outfile='layout.svg'):
        
        svg = HEADER + self.svg + FOOTER
        
        fout = open(outfile, 'w')
        fout.write(svg)
        fout.close()



layout = LayoutSVG()


source_dir = '/media/extend2/CMS/data/2017_04Apr_01/waxs/analysis/qr_image/'
#file_re = '(.+)_th(\d+\.\d+)_60.00s_.+_saxs\.jpg'
file_re = 'Name-(.+)_\d+\.\ds_T.+C_th(\d+\.\d+)_10.00s_\d+_waxs\.png'
layout.pattern_image_grid(file_re, source_dir=source_dir, col_spacing=200, row_spacing=200)


outfile = sys.argv[0][:-3]+'.svg'
layout.save(outfile)




