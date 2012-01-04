import os
from os.path import exists
import sys
import pyfits
import numpy as n
import config as c
import pymconvolve 


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
    values = line_s.split()
    mask_file = 'EM_' + c.fstring + '.fits'
    xcntr_o  = xcntr * 1.0 #x center of the object
    ycntr_o  = ycntr * 1.0 #y center of the object
    mag    = float(values[7]) #Magnitude
    radius = float(values[9]) #Half light radius
    mag_zero = c.mag_zero #magnitude zero point
    sky	 = float(values[10]) #sky 
    pos_ang = float(values[11]) #position angle
    axis_rat = 1.0 / float(values[12]) #axis ration b/a
    area_o = float(values[13]) # object area
    major_axis 	= float(values[14])	#major axis of the object
    eg = 1.0 - axis_rat
    f_tmp_mask =pyfits.open('seg.fits')
    tmp_mask = f_tmp_mask[0].data
    f_tmp_mask.close()
    for line_j in open('SegCat.cat','r'):
        try:
            values = line_j.split()
            id_n = float(values[0])
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
            if n.abs(xcntr_n - xcntr_o) < 5.0 and n.abs(ycntr_n - ycntr_o) < 5.0 and galflag == 1:
                tmp_mask[n.where(tmp_mask == id_n)] = 0
        except:
            pass
    boxcar = np.reshape(np.ones(3 * 3), (3, 3))
    tmp_mask = pymconvolve.Convolve(tmp_mask, boxcar)
    tmp_mask[n.where(tmp_mask > 0)] = 1
    if(galflag):
        hdu = pyfits.PrimaryHDU(tmp_mask.astype(n.float32))
        hdu.writeto(mask_file)
    else:
        try:
            os.remove("BMask.fits")
        except:
            pass
        hdu = pyfits.PrimaryHDU(tmp_mask.astype(n.float32))
        hdu.writeto("BMask.fits")
#line = '1    193.378    158.284 214.8573569 +56.7789966      3555934     3804.634   8.8786   0.0012     36.075     1433.745 -54.4    1.668    19672    38.968  16  0.00'
#ElliMaskFunc('n5585_lR.fits', 313, line, 0)

