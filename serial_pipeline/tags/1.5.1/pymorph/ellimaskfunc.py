import os
from os.path import exists
import sys
import pyfits
import numpy as n
import config as c
import ndimage as im

class ElliMaskFunc:
    """The class for making mask for ellipse task. This will mask all the 
       neibour objects using the masking conditions in config.py"""
    def __init__(self, cutimage, size, line_s, galflag):
        self.cutimage = cutimage
        self.size = size
        self.line_s  = line_s
        self.galflag = galflag
        self.mask    = mask(cutimage, size, line_s, galflag)

def mask(cutimage, size, line_s, galflag):
    imagefile = c.imagefile
    sex_cata = c.sex_cata
    clus_cata = c.out_cata
    threshold = c.threshold
    thresh_area = c.thresh_area
    mask_reg = c.mask_reg
    x = n.reshape(n.arange(size*size),(size,size)) % size
    x = x.astype(n.float32)
    y = n.reshape(n.arange(size*size),(size,size)) / size
    y = y.astype(n.float32)
    values = line_s.split()
    mask_file = 'EM_' + str(cutimage)[:-5] + '.fits'
    xcntr_o  = float(values[1]) #x center of the object
    ycntr_o  = float(values[2]) #y center of the object
    if c.galcut:
        xcntr = xcntr_o
        ycntr = ycntr_o
    else:
        xcntr = size / 2.0 + 1.0 + xcntr_o - int(xcntr_o)
        ycntr = size / 2.0 + 1.0 + ycntr_o - int(ycntr_o)
    
    mag    = float(values[7]) #Magnitude
    radius = float(values[9]) #Half light radius
    mag_zero = c.mag_zero #magnitude zero point
    sky	 = float(values[10]) #sky 
    pos_ang = float(values[11]) #position angle
    axis_rat = 1.0 / float(values[12]) #axis ration b/a
    area_o = float(values[13]) # object area
    major_axis 	= float(values[14])	#major axis of the object
    eg = 1.0 - axis_rat
    one_minus_eg_sq    = (1.0 - eg)**2.0
    si = n.sin(pos_ang * n.pi / 180.0)
    co = n.cos(pos_ang * n.pi / 180.0)
    tx = (x - xcntr + 0.5) * co + (y - ycntr + 0.5) * si
    ty = (xcntr - 0.5 -x) * si + (y - ycntr + 0.5) * co
    R = n.sqrt(tx**2.0 + ty**2.0 / one_minus_eg_sq)
    tmp_mask = n.zeros((size, size))
    f = pyfits.open(cutimage)
    galaxy = f[0].data
    f.close()
    startR = 3.0
    while startR < size/2.0:
        galaxymax = galaxy[n.where(R < startR)].max()
        tmp_mask[n.where(galaxy > galaxymax)] = 1
        galaxy[n.where(tmp_mask == 1)] = 0.0
        galaxy[n.where(R < startR)] = 0.0
        startR += 1.0
    tmp_mask = im.binary_fill_holes(tmp_mask)
    tmp_mask = im.binary_erosion(tmp_mask)
    z = n.zeros((size,size))		
    z1 = n.zeros((size,size))
#    print z
    for line_j in open(sex_cata,'r'):
        try:
            values = line_j.split()
            xcntr_n  = float(values[1]) #x center of the neighbour
            ycntr_n  = float(values[2]) #y center of the neighbour
            mag    = float(values[7]) #Magnitude
            radius = float(values[9]) #Half light radius
            sky      = float(values[10]) #sky
            pos_ang = float(values[11]) #position angle
            axis_rat = 1.0/float(values[12]) #axis ration b/a
            si = n.sin(pos_ang * n.pi / 180.0)
            co = n.cos(pos_ang * n.pi / 180.0)
            area = float(values[13])
            maj_axis = float(values[14])#major axis of neighbour
            eg = 1.0 - axis_rat
            one_minus_eg_sq    = (1.0-eg)**2.0
            if(abs(xcntr_n - xcntr_o) < size/2.0 and \
               abs(ycntr_n - ycntr_o) < size/2.0 and \
               xcntr_n != xcntr_o and ycntr_n != ycntr_o and galflag == 1):
                if((xcntr_o - xcntr_n) < 0):
                    xn = xcntr + abs(xcntr_n - xcntr_o)
                if((ycntr_o - ycntr_n) < 0):
                    yn = ycntr + abs(ycntr_n - ycntr_o)
                if((xcntr_o - xcntr_n) > 0):
                    xn = xcntr - (xcntr_o -xcntr_n)
                if((ycntr_o - ycntr_n) > 0):
                    yn = ycntr - (ycntr_o -ycntr_n)
                tx = (x - xn + 0.5) * co + (y - yn + 0.5) * si
                ty = (xn - 0.5 -x) * si + (y - yn + 0.5) * co
                R = n.sqrt(tx**2.0 + ty**2.0 / one_minus_eg_sq)
                z[n.where(R<=mask_reg*maj_axis)] = 1
                z1[n.where(R<=2*mask_reg*maj_axis)] = 1
            if(abs(xcntr_n - xcntr_o) < size/2.0 and \
               abs(ycntr_n - ycntr_o) < size/2.0 and galflag == 0):  
                if((xcntr_o - xntr_n) < 0):
                    xn = xcntr + abs(xcntr_n - xcntr_o)
                if((ycntr_o - ycntr_n) < 0):
                    yn = ycntr + abs(ycntr_n - ycntr_o)
                if((xcntr_o - xcntr_n) > 0):
                    xn = xcntr - (xcntr_o -xcntr_n)
                if((ycntr_o - ycntr_n) > 0):
                    yn = ycntr - (ycntr_o -ycntr_n)
                if(xcntr_n == xcntr_o and ycntr_n == ycntr_o):
                    xn = xcntr 
                    yn = ycntr
                tx = (x - xn + 0.5) * co + (y - yn + 0.5) * si
                ty = (xn - 0.5 -x) * si + (y - yn + 0.5) * co
                R = n.sqrt(tx**2.0 + ty**2.0 / one_minus_eg_sq)
                z[n.where(R<=mask_reg*maj_axis*2.0)] = 1
                z1[n.where(R<=2*mask_reg*maj_axis)] = 1
        except:
            i=1
    if(galflag):
        if exists("TmpElliMask1.fits"):
            os.remove("TmpElliMask1.fits")
        else:
            pass
        hdu = pyfits.PrimaryHDU(z1.astype(n.float32))
        hdu.writeto("TmpElliMask1.fits")
        hdu = pyfits.PrimaryHDU(z.astype(n.float32))
        hdu.writeto("TmpElliMask.fits")
        z = z + tmp_mask
        z[n.where(z > 0)] = 1
        hdu = pyfits.PrimaryHDU(z.astype(n.float32))
        hdu.writeto(mask_file)
    else:
        try:
            os.remove("BMask.fits")
        except:
            pass
        hdu = pyfits.PrimaryHDU(z.astype(n.float32))
        hdu.writeto("BMask.fits")
#line = '1    193.378    158.284 214.8573569 +56.7789966      3555934     3804.634   8.8786   0.0012     36.075     1433.745 -54.4    1.668    19672    38.968  16  0.00'
#ElliMaskFunc('n5585_lR.fits', 313, line, 0)

