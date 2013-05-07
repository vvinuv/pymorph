import os
from os.path import exists
import sys
import pyfits
import numpy as n
import config as c
import ndimage as im
from mask_or_fit import *

class ElliMaskFunc:
    """The class for making mask for ellipse task. This will mask all the 
       neibour objects using the masking conditions in config.py"""
    def __init__(self, cutimage, xcntr, ycntr, NXPTS, NYPTS, line_s, galflag):
        self.cutimage = cutimage
        self.xcntr = xcntr
        self.ycntr = ycntr
        self.NXPTS = NXPTS
        self.NYPTS = NYPTS
        self.line_s  = line_s
        self.galflag = galflag
        self.mask    = emask(cutimage, xcntr, ycntr, NXPTS, NYPTS, line_s, galflag)

def emask(cutimage, xcntr, ycntr, NXPTS, NYPTS, line_s, galflag):
    imagefile = c.imagefile
    sex_cata = c.sex_cata
    clus_cata = c.out_cata
    threshold = c.threshold
    thresh_area = c.thresh_area
    mask_reg = c.mask_reg
    avoidme = c.avoidme
    mag_zero = c.mag_zero #magnitude zero point

    x = n.reshape(n.arange(NXPTS * NYPTS),(NXPTS, NYPTS)) / NYPTS
    x = x.astype(n.float32)
    y = n.reshape(n.arange(NXPTS * NYPTS),(NXPTS, NYPTS)) % NYPTS
    y = y.astype(n.float32)

    target = SEx_obj(NXPTS, NYPTS, line_s)
    R =target.calc_rad(x,y)
    
    mask_file = 'EM_' + c.fstring + '.fits'

    tmp_mask = n.zeros((NXPTS, NYPTS))
    f = pyfits.open(c.datadir + cutimage)
    galaxy = f[0].data
    f.close()
    galaxy = n.swapaxes(galaxy, 0, 1)

    startR = 3.0
    while startR < 4*target.radius: #max(NXPTS / 2.0, NYPTS / 2.0):
        galaxymax = galaxy[n.where(R < startR)].max()
        tmp_mask[n.where(galaxy > galaxymax)] = 1
        galaxy[n.where(tmp_mask == 1)] = 0.0
        galaxy[n.where(R < startR)] = 0.0
        startR += 1.0
        
    tmp_mask = im.binary_fill_holes(tmp_mask)
    tmp_mask = im.binary_erosion(tmp_mask)
    z = n.zeros((NXPTS, NYPTS))		
    z1 = n.zeros((NXPTS, NYPTS))


    for line_j in open(sex_cata,'r'):
        if line_j[0] != '#': #line is not a comment
            
            neighbor = SEx_obj(NXPTS, NYPTS, line_j)
            if target.mask_or_fit(neighbor,threshold,thresh_area,avoidme)==1:
                # neighbor.set_axis_rat(1.0) #make masks circular

                # recenter in chip coordinates
                xn = xcntr - target.xcntr + neighbor.xcntr
                yn = ycntr - target.ycntr + neighbor.ycntr
                neighbor.set_center(xn, yn)
                
                R = neighbor.calc_rad(x,y)
                if galflag ==1:
                    z[n.where(R<=mask_reg*neighbor.maj_axis)] = 1
                    z1[n.where(R <= 2 * mask_reg * neighbor.maj_axis)] = 1
                elif galflag == 0:  
                    z[n.where(R <= neighbor.maj_axis * 2.5)] = 1
        
    if(galflag):
        if exists("TmpElliMask1.fits"):
            os.remove("TmpElliMask1.fits")
        hdu = pyfits.PrimaryHDU(n.swapaxes(z1, 0, 1).astype(n.float32))
        hdu.writeto("TmpElliMask1.fits")
	if c.NoMask:
	    z[n.where(z > 0)] = 0
        elif c.NormMask:
            pass
        else:
            z = z + tmp_mask
            z[n.where(z > 0)] = 1
            
        if exists("TmpElliMask.fits"):
            os.remove("TmpElliMask.fits")
        hdu = pyfits.PrimaryHDU(n.swapaxes(z, 0, 1).astype(n.float32))
        hdu.writeto("TmpElliMask.fits")

        hdu = pyfits.PrimaryHDU(n.swapaxes(z, 0, 1).astype(n.float32))
        hdu.writeto(mask_file)
    else:
        try:
            os.remove("BMask.fits")
        except:
            pass
 	if c.NoMask:
	    z[n.where(z > 0)] = 0
        elif c.NormMask:
            pass
        else:
            z = z + tmp_mask
            z[n.where(z > 0)] = 1
        z = im.binary_dilation(z, iterations=15)
        hdu = pyfits.PrimaryHDU(n.swapaxes(z, 0, 1).astype(n.float32))
        hdu.writeto("BMask.fits")
        os.system("cp BMask.fits B.fits")
