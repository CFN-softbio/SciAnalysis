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
        
        
    def add_text_c(self, text, xc, yc, color='FFFFFF', fontsize=30):
        
        svg = """    <text
       sodipodi:linespacing="100%"
       id="text{:d}"
       y="{:.4f}"
       x="{:.4f}"
       style="font-size:{:.6f}px;font-style:normal;font-weight:normal;line-height:100%;letter-spacing:0px;word-spacing:0px;fill:#{};fill-opacity:1;stroke:none;font-family:Bitstream Vera Sans"
       xml:space="preserve"><tspan
         style="text-align:center;line-height:100%;writing-mode:lr-tb;text-anchor:middle"
         y="{:.4f}"
         x="{:.4f}"
         id="tspan{:d}"
         sodipodi:role="line">{}</tspan></text>
""".format(self.txt_id, yc, xc, fontsize, color, yc, xc, self.txt_id+1, text)
        
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
        
        
        
    def pattern_image_grid(self, file_re, source_dir='./', col_spacing=200, row_spacing=None, x0=0, y0=0, verbosity=3):
        
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
                    
                    
        #col_strs = sorted(col_strs, key=int)
        #row_strs = sorted(row_strs, key=int)
        col_strs = sorted(col_strs)
        row_strs = sorted(row_strs)
                    
        # Col headers
        for ic, col_str in enumerate( col_strs ):
            
            x = x0 + ic*col_spacing
            y = y0 - float(row_spacing)*( len(row_strs) - 1.0 + 0.5 )
            
            col_str = 'ab '+col_str+''
            
            self.add_text(col_str, x, y, fontsize=14)
        
        # Row headers
        for ir, row_str in enumerate( row_strs ):
            
            x = x0 - float(col_spacing)/2
            y = y0 + ir*row_spacing
            
            row_str = 'cc '+row_str+''

            self.add_vtext(row_str, x, y, fontsize=25)
            
            
        # Images
        for ic, col_str in enumerate( col_strs ):
            for ir, row_str in enumerate( row_strs ):
                
                if verbosity>=5:
                    print("  Checking col='{}' row='{}'".format(col_str, row_str))
                
                matching = [s for s in filenames if col_str in s]
                matching = [s for s in matching if row_str in s]
                
                if len(matching)<1:
                    if verbosity>=1:
                        print("  WARNING: No file for col='{}' row='{}'".format(col_str, row_str))
                    
                else:
                
                    if len(matching)>1:
                        
                        filename = None
                        if verbosity>=4:
                            print("  NOTE: {} matches for col='{}' row='{}'".format(len(matching), col_str, row_str))

                        for ifile, filename_cur in enumerate(matching):
                            if verbosity>=5:
                                print("    Checking {}: {}".format(ifile, filename_cur))
                            m = match_re.match(filename_cur)
                            if m:
                                if col_str==m.groups()[0] and row_str==m.groups()[1]:
                                    if verbosity>=5:
                                        print("      Confirmed match {}: {}".format(ifile, filename_cur))
                                    filename = filename_cur
                                    break
                            else:
                                if verbosity>=1:
                                    print("  ERROR: File {} didn't match RE (that should be impossible)".format(filename_cur))
                            
                        True
                        
                    else:
                        filename = matching[0]
                        
                    #filename = matching[0]
                    
                    
                    x = x0 + ic*col_spacing
                    y = y0 - ir*row_spacing
                    if verbosity>=3:
                        print("  Image: {}".format(filename))
                    if verbosity>=4:
                        print("    col='{}', row='{}' yields '{}'".format(col_str, row_str, filename))
                        
                    if filename is None:
                        if verbosity>=1:
                            print("  ERROR: filename is None for col='{}' row='{}'".format(len(matching), col_str, row_str))
                    else:
                        self.add_image(source_dir+filename, x, y, width=col_spacing)
                
                
        


    def match(self, filename, match_names=None, verbosity=3):
        '''Used to internally filter what files will actually be considered.
        One can tweak this method to enforce different matching behavior.'''
        
        
        if match_names is None:
            return True
        
        else:
            for match_name in match_names:
                if match_name in filename:
                    return True
            
        return False
        

        
    def pattern_multi_series(self, infiles, parse_re, match_names=None, source_dir='./', col_spacing=200, row_spacing=None, x0=0, y0=0, parse_row_idx=0, parse_col_idx=1, overlay_text_color='000000', verbosity=3):
        
        if row_spacing is None:
            row_spacing = col_spacing
        
        match_re = re.compile(parse_re)
        
        # Figure out the list of rows
        filenames = []
        row_strs = []
        col_strs = []
        for infile in infiles:
            
            f = Filename(infile)
            filepath, filename, filebase, ext = f.split()
            
            
            m = match_re.match(filename)
            if self.match(filename, match_names=match_names) and m:
                #print(filename)
                
                if filename in filenames:
                    print('ERROR: File already added.')
                filenames.append(filename)
                
                row_str = m.groups()[parse_row_idx]
                col_str = m.groups()[parse_col_idx]
                
                if row_str not in row_strs:
                    row_strs.append(row_str)

                if col_str not in col_strs:
                    col_strs.append(col_str)
                    
            else:
                if verbosity>=3:
                    print('  Note: File does not match RE: {}'.format(infile))
                    
                    
        if verbosity>=2:
            print('================================================================================')
            print('{} filenames matched the RE (out of {} infiles); {:.1f}%'.format(len(filenames), len(infiles), 100.*len(filenames)/len(infiles)))
            print('================================================================================')
                  
                  
        
        row_strs = sorted(row_strs)
        #row_strs = sorted(row_strs, key=int)
        col_strs = sorted(col_strs)

        if verbosity>=2:
            print('{} rows identified'.format(len(row_strs)))

            
        # Images
        max_cols = 0
        for ir, row_str in enumerate( row_strs ):
            
            
            matching_files = []
            col_strs = []
            for filename in filenames:
                m = match_re.match(filename)
                if m:
                    if m.groups()[parse_row_idx] == row_str:
                        matching_files.append(filename)
                        col_strs.append(m.groups()[parse_col_idx])
                
            col_strs = sorted(col_strs, key=float)
            if len(col_strs)>max_cols:
                max_cols = len(col_strs)
            
            for ic, col_str in enumerate( col_strs ):
            
                if verbosity>=5:
                    print("  Checking col='{}' row='{}'".format(col_str, row_str))
                
                
                filename = None
                
                for match_file in matching_files:
                    m = match_re.match(match_file)
                    if m:
                        if m.groups()[parse_col_idx] == col_str:
                            filename = match_file
                            break
                    
                if filename is None:
                    if verbosity>=1:
                        print("  ERROR: filename is None for col='{}' row='{}'".format(col_str, row_str))
                
                else:
                    x = x0 + ic*col_spacing
                    y = y0 - ir*row_spacing
                    if verbosity>=3:
                        print("  Image: {}".format(filename))
                    if verbosity>=4:
                        print("    col='{}', row='{}' yields '{}'".format(col_str, row_str, filename))
                        
                    self.add_image(source_dir+filename, x, y, width=col_spacing)
                    
                    s = filename
                    self.add_text_c(s, x, y-row_spacing*0.5*0.9, color=overlay_text_color, fontsize=14)
                

        
        # Row headers
        for ir, row_str in enumerate( row_strs ):
            
            x = x0 - float(col_spacing)/2
            y = y0 + ir*row_spacing
            
            #row_str = 'cc '+row_str+''

            self.add_vtext(row_str, x, y, fontsize=80)
            

                        
        # Col headers
        if verbosity>=2:
            print('{} columns'.format(max_cols))
        # Note the col_strs here have been reset within the Image loop above.
        # So this will print the col_strs associated with the last processed row.
        for ic, col_str in enumerate( col_strs ):
            
            x = x0 + ic*col_spacing
            y = y0 - float(row_spacing)*( len(row_strs) - 1.0 + 0.5 )
            
            #col_str = 't '+str(ic)+''
            
            self.add_text(col_str, x, y, fontsize=80)                                    
        
        
        
    
    def save(self, outfile='layout.svg'):
        
        svg = HEADER + self.svg + FOOTER
        
        fout = open(outfile, 'w')
        fout.write(svg)
        fout.close()






layout = LayoutSVG()


if True:
    # The prefix that all the samples match
    sample_names = 'runT'
    protocol = 'thumbnails'
    
    match_names = None
    #match_names = ['sample1', 'sample2'] # Turn this on to select specific samples to match
    
    source_dir = '../../analysis/{}'.format(protocol)
    pattern = '{}*_th*_saxs.*'.format(sample_names) 
    parse_re = '^(.+)_th\d\.\d+_.+_(\d+)_saxs\.[a-zA-Z]+$' # In this RE, the bracketed regions denote the parts that will be used for row and column sorting.
    
    col_spacing, row_spacing = 487, 619+15 
    outfile = '{}-{}.svg'.format(sample_names, protocol)
    
    



source_dir = os.path.abspath(source_dir) + '/'
infiles = glob.glob(os.path.join(source_dir, pattern))
infiles.sort()
print('Considering {} files...'.format(len(infiles)))

layout.pattern_multi_series(infiles, parse_re, match_names=match_names, source_dir=source_dir, col_spacing=col_spacing, row_spacing=row_spacing, verbosity=2)

layout.save(outfile)



