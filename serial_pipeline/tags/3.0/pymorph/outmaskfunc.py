import os, sys, pyfits
import numpy as n
import config as c

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
    x = n.reshape(n.arange(NXPTS * NYPTS),(NXPTS, NYPTS)) / NYPTS
    x = x.astype(n.float32)
    y = n.reshape(n.arange(NXPTS * NYPTS),(NXPTS, NYPTS)) % NYPTS
    y = y.astype(n.float32)
    values = line_s.split()
    mask_file = 'OEM_' + c.fstring + '.fits'
    xcntr_o  = float(values[1]) #x center of the object
    ycntr_o  = float(values[2]) #y center of the object
    mag    = float(values[7]) #Magnitude
    radius = float(values[9]) #Half light radius
    mag_zero = c.mag_zero #magnitude zero point
    sky	 = float(values[10]) #sky 
    pos_ang = float(values[11]) - 90.0 #position angle
    axis_rat = 1.0 / float(values[12]) #axis ration b/a
    area_o = float(values[13])   # object's area
    major_axis = float(values[14])	#major axis of the object
    z = n.zeros((NXPTS, NYPTS))
    for line_j in open(sex_cata,'r'):
        try:
            values = line_j.split()
            xcntr_n  = float(values[1]) #x center of the neighbour
            ycntr_n  = float(values[2]) #y center of the neighbour
            mag    = float(values[7]) #Magnitude
            radius = float(values[9]) #Half light radius
            sky      = float(values[10]) #sky
            pos_ang = float(values[11]) - 90.0 #position angle
            axis_rat = 1.0/float(values[12]) #axis ration b/a
            area_n = float(values[13])
            maj_axis = float(values[14])#major axis of neighbour
            if(abs(xcntr_n - xcntr_o) <= (major_axis + maj_axis) * \
               threshold and \
               abs(ycntr_n - ycntr_o) <= (major_axis  + maj_axis) * \
               threshold and area_n >= thresh_area * area_o and \
               xcntr_n != xcntr_o and ycntr_n != ycntr_o):
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
                z[n.where(R<=mask_reg*maj_axis)] = 1
        except:
            pass	
    hdu = pyfits.PrimaryHDU(n.swapaxes(z, 0, 1).astype(n.float32))
    hdu.writeto(mask_file)
