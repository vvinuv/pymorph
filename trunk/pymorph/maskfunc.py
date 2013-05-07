import os, sys, pyfits
import numpy as n
import config as c
import ndimage as im
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
    avoidme = c.avoidme
    mask_reg = c.mask_reg
    x = n.reshape(n.arange(NXPTS * NYPTS),(NXPTS, NYPTS)) / NYPTS
    x = x.astype(n.float32)
    y = n.reshape(n.arange(NXPTS * NYPTS),(NXPTS, NYPTS)) % NYPTS
    y = y.astype(n.float32)
    
    target = SEx_obj(NXPTS, NYPTS, line_s)
    R = target.calc_rad(x,y)
    
    mask_file = 'M_' + c.fstring + '.fits'
    tmp_mask = n.zeros((NXPTS, NYPTS))
    f = pyfits.open(c.datadir + cutimage)
    galaxy = f[0].data
    f.close()
    galaxy = n.swapaxes(galaxy, 0, 1)

    # search for and mask unusually bright pixels, assuming a smooth profile
    startR = 3.0 # the pixel radius that we start at
    
    while startR < 4*target.radius: #max(NXPTS / 2.0, NYPTS / 2.0): 
        galaxymax = galaxy[n.where(R < startR)].max()
        tmp_mask[n.where(galaxy > galaxymax)] = 1
        galaxy[n.where(tmp_mask == 1)] = 0.0
        galaxy[n.where(R < startR)] = 0.0
        startR += 1.0 # step through radii in single pixel annuli

    tmp_mask = im.binary_fill_holes(tmp_mask)
    tmp_mask = im.binary_erosion(tmp_mask)
    f = pyfits.open('TmpElliMask1.fits')
    ellip_mask = f[0].data
    f.close()
    ellip_mask = n.swapaxes(ellip_mask, 0, 1)

    tmp_mask[n.where(ellip_mask == 1)] = 0
    z = n.zeros((NXPTS, NYPTS))

    for line_j in open(sex_cata,'r'):
        if line_j[0] != '#': #line is not a comment
            
            neighbor = SEx_obj(NXPTS, NYPTS, line_j)
            if target.mask_or_fit(neighbor,threshold,thresh_area,avoidme)==1:
                #neighbor.set_axis_rat(1.0) #make masks circular

                # recenter in chip coordinates
                #print 'masking ', line_j
                xn = xcntr - target.xcntr + neighbor.xcntr
                yn = ycntr - target.ycntr + neighbor.ycntr
                #print neighbor.xcntr, neighbor.ycntr
                #print xn, yn
                
                neighbor.set_center(xn, yn)

                R = neighbor.calc_rad(x,y)
                z[n.where(R<=mask_reg*neighbor.maj_axis)] = 1
            #elif target.mask_or_fit(neighbor,threshold,thresh_area,avoidme)==0:
                #print "source ", line_j, " should be fit"
                
    if c.NoMask:
        z[n.where(z > 0)] = 0
    elif c.NormMask:
        pass
    else:
        z = z + tmp_mask
        z[n.where(z > 0)] = 1
        z = im.binary_dilation(z, iterations=2)
    z = im.binary_fill_holes(z)
    hdu = pyfits.PrimaryHDU(n.swapaxes(z, 0, 1).astype(n.float32))
    hdu.writeto(mask_file)
    try:
        os.remove('TmpElliMask.fits')
    except:
        pass
