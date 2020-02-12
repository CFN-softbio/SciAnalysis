#!/usr/bin/python3

import glob
import re
import random

import numpy as np
#import scipy.misc
import imageio
from PIL import Image



#width, height, desig = 3840, 2160, '4k' # 4k UHD
width, height, desig = 7680, 2160, '-4kW' # 4k UHD wide


if False:
    source_dir = '../thumbnails2/'
    ext = 'jpg'
    infiles = glob.glob(source_dir+'PWM*x*yy*20.00s*_saxs.{}'.format(ext))
    print('{} infiles'.format(len(infiles)))

    filename_re = re.compile('^.+_x(-?\d+\.\d+)_yy(-?\d+\.\d+)_.+_(\d+)_saxs\.{}$'.format(ext))

    positions = []
    sIDs = []
    for infile in infiles:
        
        m = filename_re.match(infile)
        if m:
            x = float(m.groups()[0])
            y = float(m.groups()[1])
            positions.append( [x,y] )
            
            sIDs.append( int(m.groups()[2]) )


    print('{} matched files'.format(len(sIDs)))
    positions = np.asarray(positions)
    sIDs = np.asarray(sIDs)
    infiles = np.asarray(infiles)

    # Sort
    idx = np.argsort(sIDs)
    sIDs = sIDs[idx]
    infiles = infiles[idx]
    positions = positions[idx]
    


if True:
    source_dir = '/media/extend2/CMS/'
    infiles = glob.glob(source_dir+'**/thumbnails/*.png', recursive=True)
    #infiles = glob.glob(source_dir+'**/thumbnails/*.jpg', recursive=True)
    print('{} infiles'.format(len(infiles)))

    random.shuffle(infiles)
    
    infiles = infiles[:3000]


    aspect = width/height # w/h, x/y, col/row
    num_col = int(np.floor( np.sqrt(len(infiles)*aspect) ))
    num_row = len(infiles)//num_col
    print('{} columns by {} rows = {:,} elements ({:.1f}% of {:,} total)'.format(num_col, num_row, num_col*num_row, 100.*num_col*num_row/len(infiles), len(infiles)))

    positions = []
    for icol in range(num_col):
        for irow in range(num_row):
            positions.append( [icol, irow] )

    positions = np.asarray(positions)



# For testing/debugging
#print(infiles[0])
#infiles = infiles[:100]
#exit()


def sm_imresize(image, size, handle_float=False):
    '''Replacement for deprecated scipy.misc.imresize function.'''
    #image = scipy.misc.imresize(image, size) # Deprecated
    
    h, w, c = image.shape
    if isinstance(size, (int, float)):
        hn = int(h*size)
        wn = int(w*size)
    elif len(size)==2:
        hn, wn = size
    else:
        print('Error in sm_imresize.')
    
    if handle_float:
        image = np.copy(image)*255
        image = np.array( Image.fromarray( image.astype(np.uint8) ).resize((wn,hn)) )
        image = image/255
    else:
        image = np.array( Image.fromarray( image.astype(np.uint8) ).resize((wn,hn)) )
        #image = resize(image, output_shape=(hn,wn), preserve_range=True) # Doesn't work
                        
    return image


class Canvas():
    
    def __init__(self, width=3840, height=2160, fill=0):
        self.width = width
        self.height = height
        
        self.image = np.ones( (self.height, self.width, 3), dtype=np.uint8 )*fill
        
        self.im_scaling, self.canvas_scaling = None, None


    def compute_canvas_transform(self, positions, im_width=500, fit=True):
        
        xs = np.asarray(positions)[:,0]
        ys = np.asarray(positions)[:,1]
        
        xspan = np.max(xs)-np.min(xs)
        yspan = np.max(ys)-np.min(ys)
        
        xscale = self.width/xspan
        yscale = self.height/yspan
        
        self.origin = [ np.average(xs), np.average(ys) ]
        
        if fit:
            self.canvas_scaling = min(xscale, yscale) # Fit all data into space
        else:
            self.canvas_scaling = max(xscale, yscale) # Leave no gaps in image
        
        
        aspect = self.width/self.height
        num_col = int(np.floor( np.sqrt(len(positions)*aspect) ))
        num_row = len(positions)//num_col
        
        self.im_scaling = ( self.width/(num_col) )/im_width
        
        
        
    def add_images(self, infiles, positions, canvas_scaling=None, im_scaling=None, im_scaling_tweak=1.0, force_square=True, force_size=None, verbosity=3):
        
        for i, (infile, (x,y)) in enumerate(zip(infiles, positions)):
            
            if verbosity>=3 and i%100==0:
                print('    Adding image {:d}/{:d} ({:.1f}% complete)'.format(i, len(infiles), 100*i/len(infiles)))
            
            if verbosity>=4:
                print('        Adding {}'.format(infile))
            self.add_image(infile, x, y, canvas_scaling=canvas_scaling, im_scaling=im_scaling, im_scaling_tweak=im_scaling_tweak, force_square=force_square, force_size=force_size, verbosity=verbosity)
            
        
        
    def add_image(self, infile, x, y, canvas_scaling=None, im_scaling=None, im_scaling_tweak=1.0, force_square=True, force_size=None, verbosity=3):
        
        im_scaling = im_scaling if im_scaling is not None else self.im_scaling
        im_scaling *= im_scaling_tweak
        canvas_scaling = canvas_scaling if canvas_scaling is not None else self.canvas_scaling
        
        #im = scipy.misc.imread(infile) # Deprecated
        im = imageio.imread(infile)
        
        h, w, c = im.shape
        

        if force_square:
            if w>h:
                mid = w//2
                im = im[:,mid-h//2:mid+h//2,:]
            elif h>w:
                mid = h//2
                im = im[mid-w//2:mid+w//2,:,:]
                
        if force_size is not None:
            #im = scipy.misc.imresize(im, [force_size, force_size]) # Deprecated
            im = sm_imresize(im, [force_size, force_size])
            
        if abs(im_scaling-1.0)>0.01:
            #im = scipy.misc.imresize(im, im_scaling) # Deprecated
            im = sm_imresize(im, im_scaling)
            h, w, c = im.shape
                
        
        xc = (x-self.origin[0])*canvas_scaling + self.width/2
        yc = (y-self.origin[1])*canvas_scaling + self.height/2
        
        xi = int(xc-w/2)
        xf = xi+w
        yi = int(yc-h/2)
        yf = yi+h
        
        # Target (canvas)
        xit = np.clip(xi, 0, self.width)
        xft = np.clip(xf, 0, self.width)
        yit = np.clip(yi, 0, self.height)
        yft = np.clip(yf, 0, self.height)
        
        # Source (image)
        #xis, xfs, yis, yfs = 0, w, 0, h # Assuming target is strictly within canvas (no clipping)
        xis = xit-xi
        xfs = w-(xf-xft)
        yis = yit-yi
        yfs = h-(yf-yft)
        
        if verbosity>=5:
            print('            Copying from im{} to canvas{}'.format( (xis, xfs, yis, yfs), (xit, xft, yit, yft) ))
        
        self.image[yit:yft,xit:xft,:] = im[yis:yfs,xis:xfs,:3]
        
        
        
        
    def save(self, outfile='output.png'):
        #scipy.misc.imsave(outfile, self.image) # Deprecated
        imageio.imsave(outfile, self.image)
        
        
        




canvas = Canvas(width=width, height=height, fill=0)

canvas.compute_canvas_transform(positions, im_width=500, fit=True)
canvas.add_images(infiles, positions, force_square=True, force_size=500, im_scaling_tweak=1, verbosity=3)

canvas.save('Xray_layout{}.png'.format(desig))






