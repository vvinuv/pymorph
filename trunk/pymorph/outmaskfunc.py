import os, sys, pyfits
import numpy as n
import config as c
from mask_or_fit import *

class OutMaskFunc:
    """The class for making mask for output image from GALFIT to run ellipse.
       The file name of the mask is OEM_string(galid).fits"""
    def __init__(self, outimage, xcntr, ycntr, NXPTS, NYPTS, line_s):
        self.outimage = outimage
        self.xcntr = xcntr
        self.ycntr = ycntr
        self.NXPTS = NXPTS
        self.NYPTS = NYPTS
        self.line_s  = line_s
        self.mask    = mask(outimage, xcntr, ycntr, NXPTS, NYPTS, line_s)

def mask(outimage, xcntr, ycntr, NXPTS, NYPTS, line_s):
    imagefile = c.imagefile
    sex_cata = c.sex_cata
    threshold = c.threshold
    thresh_area = c.thresh_area
    mask_reg = c.mask_reg
    avoidme = c.avoidme
    
    x = n.reshape(n.arange(NXPTS * NYPTS),(NXPTS, NYPTS)) / NYPTS
    x = x.astype(n.float32)
    y = n.reshape(n.arange(NXPTS * NYPTS),(NXPTS, NYPTS)) % NYPTS
    y = y.astype(n.float32)


    target = SEx_obj(NXPTS, NYPTS, line_s)

    mask_file = 'OEM_' + c.fstring + '.fits'
    z = n.zeros((NXPTS, NYPTS))
    
    for line_j in open(sex_cata,'r'):
        if line_j[0] != '#': #line is not a comment
            neighbor = SEx_obj(NXPTS, NYPTS, line_j)
            
            if target.mask_or_fit(neighbor,threshold,thresh_area,avoidme)>=0:
                #neighbor.set_axis_rat(1.0) #make masks circular

                # recenter in chip coordinates
                xn = xcntr - target.xcntr + neighbor.xcntr
                yn = ycntr - target.ycntr + neighbor.ycntr
                neighbor.set_center(xn, yn)

                R = neighbor.calc_rad(x,y)
                z[n.where(R<=mask_reg*neighbor.maj_axis)] = 1
	
    hdu = pyfits.PrimaryHDU(n.swapaxes(z, 0, 1).astype(n.float32))
    hdu.writeto(mask_file)
