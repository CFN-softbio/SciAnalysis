#!/usr/bin/python3

import numpy as np
import PIL.Image

#infile = 'CHX_Eiger1M_blemish2-mask.npy'
infile = 'CHX_Eiger1M_flatfield.npy'

outfile = infile[:-4] + '.png'

pixmask = np.load(infile)
#pixmask = pixmask < 1

#img = np.where( pixmask<1, 255, 0)
#img = np.where( pixmask<1, 0, 255)
img = np.where( pixmask<0.1, 0, 255)


img = PIL.Image.fromarray(np.uint8(img))

img.save(outfile)