import os, sys, pyfits
import numpy as n
import config as c
import ndimage as im

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
    x = n.reshape(n.arange(NXPTS * NYPTS),(NXPTS, NYPTS)) / NYPTS
    x = x.astype(n.float32)
    y = n.reshape(n.arange(NXPTS * NYPTS),(NXPTS, NYPTS)) % NYPTS
    y = y.astype(n.float32)
    values = line_s.split()
    mask_file = 'M_' + str(cutimage)[:-5] + '.fits'
    xcntr_o  = float(values[1]) #x center of the object
    ycntr_o  = float(values[2]) #y center of the object
    mag    = float(values[7]) #Magnitude
    radius = float(values[9]) #Half light radius
    mag_zero = c.mag_zero #magnitude zero point
    sky	 = float(values[10]) #sky 
    pos_ang = float(values[11]) #position angle
    axis_rat = 1.0 / float(values[12]) #axis ration b/a
    area_o = float(values[13]) # object area
    major_axis = float(values[14])	#major axis of the object
    eg = 1.0 - axis_rat
    one_minus_eg_sq    = (1.0 - eg)**2.0
    si = n.sin(pos_ang * n.pi / 180.0)
    co = n.cos(pos_ang * n.pi / 180.0)
    tx = (x - xcntr + 1.0) * co + (y - ycntr + 1.0) * si
    ty = (xcntr - 1.0 -x) * si + (y - ycntr + 1.0) * co
    R = n.sqrt(tx**2.0 + ty**2.0 / one_minus_eg_sq)
    tmp_mask = n.zeros((NXPTS, NYPTS))
    f = pyfits.open(cutimage)
    galaxy = f[0].data
    f.close()
    galaxy = n.swapaxes(galaxy, 0, 1)
    startR = 3.0
    while startR < max(NXPTS / 2.0, NYPTS / 2.0):
        galaxymax = galaxy[n.where(R < startR)].max()
        tmp_mask[n.where(galaxy > galaxymax)] = 1
        galaxy[n.where(tmp_mask == 1)] = 0.0
        galaxy[n.where(R < startR)] = 0.0
        startR += 1.0
    tmp_mask = im.binary_fill_holes(tmp_mask)
    tmp_mask = im.binary_erosion(tmp_mask)
    f = pyfits.open('TmpElliMask1.fits')
    ellip_mask = f[0].data
    f.close()
    ellip_mask = n.swapaxes(ellip_mask, 0, 1)
    tmp_mask[n.where(ellip_mask == 1)] = 0
    z = n.zeros((NXPTS, NYPTS))
    for line_j in open(sex_cata,'r'):
        try:
            values = line_j.split()
            xcntr_n  = float(values[1]) #x center of the neighbour
            ycntr_n  = float(values[2]) #y center of the neighbour
            mag    = float(values[7]) #Magnitude
            radius = float(values[9]) #Half light radius
            sky      = float(values[10]) #sky
            pos_ang = float(values[11])  #position angle
            axis_rat = 1.0/float(values[12]) #axis ration b/a
            si = n.sin(pos_ang * n.pi / 180.0)
            co = n.cos(pos_ang * n.pi / 180.0)
            area_n = float(values[13]) #neighbour area
            maj_axis = float(values[14])#major axis of neighbour
            eg = 1.0 - axis_rat
            one_minus_eg_sq    = (1.0-eg)**2.0
            if(abs(xcntr_n - xcntr_o) < NXPTS / 2.0 + 30.0 and \
               abs(ycntr_n - ycntr_o) < NYPTS / 2.0 + 30.0 and \
               xcntr_n != xcntr_o and ycntr_n != ycntr_o):
                if(abs(xcntr_n - xcntr_o) > threshold * (major_axis + \
                   maj_axis) or abs(ycntr_n - ycntr_o) > threshold * \
                   (major_axis + maj_axis) or area_n < thresh_area * area_o):
                    if((xcntr_o - xcntr_n) < 0):
                        xn = xcntr + abs(xcntr_n - xcntr_o)
                    if((ycntr_o - ycntr_n) < 0):
                        yn = ycntr + abs(ycntr_n - ycntr_o)
                    if((xcntr_o - xcntr_n) > 0):
                        xn = xcntr - (xcntr_o - xcntr_n)
                    if((ycntr_o - ycntr_n) > 0):
                        yn = ycntr - (ycntr_o - ycntr_n)
                    tx = x - xn + 1.0 
                    ty = y - yn + 1.0
                    R = n.sqrt(tx**2.0 + ty**2.0)
#                    tx = (x - xn + 0.5) * co + (y - yn + 0.5) * si
#                    ty = (xn - 0.5 -x) * si + (y - yn + 0.5) * co
#                    R = n.sqrt(tx**2.0 + ty**2.0 / one_minus_eg_sq)
                    z[n.where(R<=mask_reg*maj_axis)] = 1
        except:
            pass
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
