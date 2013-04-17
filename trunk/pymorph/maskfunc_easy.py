import os, sys, pyfits
import numpy as np
import config as c
import pymconvolve
from mask_or_fit import *

class MaskFunc:
    """The class for making mask for GALFIT. It uses the masking conditions 
       from config.py. The output mask image will have the name
       M_string(galid).fits """
    def __init__(self, cutimage, xcntr, ycntr, NXPTS, NYPTS, line_s):
        self.cutimage = cutimage
        self.xcntr = xcntr
        self.ycntr = ycntr
        self.NXPTS = NXPTS
        self.NYPTS = NYPTS
        self.line_s  = line_s
        self.mask    = gmask(cutimage, xcntr, ycntr, NXPTS, NYPTS, line_s)

def gmask(cutimage, xcntr, ycntr, NXPTS, NYPTS, line_s):
    imagefile = c.imagefile
    sex_cata = c.sex_cata
    threshold = c.threshold
    thresh_area = c.thresh_area
    mask_reg = c.mask_reg
    avoidme = c.avoidme
    
    target = SEx_obj(NXPTS, NYPTS, line_s)
    
    mask_file = 'M_' + c.fstring + '.fits'
    mag_zero = c.mag_zero #magnitude zero point

    f_tmp_mask =pyfits.open('seg.fits')
    tmp_mask = f_tmp_mask[0].data
    f_tmp_mask.close()
    for line_j in open('SegCat.cat','r'):
        if line_j[0] != '#': #line is not a comment
        #try:
            neighbor = SEx_obj(NXPTS, NYPTS, line_j)
            
            if target.mask_or_fit(neighbor,threshold,thresh_area,avoidme)!=1:
                tmp_mask[np.where(tmp_mask == neighbor.sex_num)] = 0

        #except:
        #    pass
    if c.NoMask:
        tmp_mask[np.where(tmp_mask > 0)] = 0
    elif c.NormMask:
        pass
    boxcar = np.reshape(np.ones(3 * 3), (3, 3))
    tmp_mask = pymconvolve.Convolve(tmp_mask, boxcar)
    tmp_mask[tmp_mask > 1e-5] = 1
    tmp_mask[tmp_mask != 1] = 0
    hdu = pyfits.PrimaryHDU(tmp_mask.astype(np.float32))
    hdu.writeto(mask_file)
