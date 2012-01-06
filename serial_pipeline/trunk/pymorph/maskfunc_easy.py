import os, sys, pyfits
import numpy as np
import config as c
import pymconvolve

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
    values = line_s.split()
    mask_file = 'M_' + c.fstring + '.fits'
    xcntr_o  = xcntr * 1.0 #x center of the object
    ycntr_o  = ycntr * 1.0 #y center of the object
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
            pos_ang = float(values[11])  #position angle
            axis_rat = 1.0/float(values[12]) #axis ration b/a
            si = np.sin(pos_ang * np.pi / 180.0)
            co = np.cos(pos_ang * np.pi / 180.0)
            area_n = float(values[13]) #neighbour area
            maj_axis = float(values[14])#major axis of neighbour
            eg = 1.0 - axis_rat
            one_minus_eg_sq    = (1.0-eg)**2.0
            if(abs(xcntr_n - xcntr_o) > threshold * (major_axis + \
               maj_axis) or abs(ycntr_n - ycntr_o) > threshold * \
              (major_axis + maj_axis) or area_n < thresh_area * area_o):
                if abs(xcntr_n - xcntr_o) < 5.0 and \
                  abs(ycntr_n - ycntr_o) < 5.0:    
                    tmp_mask[np.where(tmp_mask == id_n)] = 0
                pass
            else:
                tmp_mask[np.where(tmp_mask == id_n)] = 0
        except:
            pass
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
