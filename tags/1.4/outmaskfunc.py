import os, sys, pyfits
import numpy as n
import config as c

class OutMaskFunc:
    """The class for making mask for GALFIT """
    def __init__(self, outimage, size, line_s):
        self.outimage = outimage
        self.size = size
        self.line_s  = line_s
        self.mask    = mask(outimage, size, line_s)

def mask(outimage, size, line_s):
    imagefile = c.imagefile
    sex_cata = c.sex_cata
    threshold = c.threshold
    thresh_area = c.thresh_area
    mask_reg = c.mask_reg
    x = n.reshape(n.arange(size*size),(size,size)) % size
    x = x.astype(n.float32)
    y = n.reshape(n.arange(size*size),(size,size)) / size
    y = y.astype(n.float32)
    values = line_s.split()
    mask_file = 'OEM_' + str(outimage)[:-5] + '.fits'
    xcntr_o  = float(values[1]) #x center of the object
    ycntr_o  = float(values[2]) #y center of the object
    xcntr = size / 2.0 + 1.0 + xcntr_o - int(xcntr_o)
    ycntr = size / 2.0 + 1.0 + ycntr_o - int(ycntr_o)
    mag    = float(values[7]) #Magnitude
    radius = float(values[9]) #Half light radius
    mag_zero = c.mag_zero #magnitude zero point
    sky	 = float(values[10]) #sky 
    pos_ang = float(values[11]) - 90.0 #position angle
    axis_rat = 1.0 / float(values[12]) #axis ration b/a
    major_axis = float(values[14])	#major axis of the object
    z = n.zeros((size,size))
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
            area = float(values[13])
            maj_axis = float(values[14])#major axis of neighbour
            if(abs(xcntr_n - xcntr_o) < major_axis * threshold and \
               abs(ycntr_n - ycntr_o) < major_axis * threshold and \
               xcntr_n != xcntr_o and ycntr_n != ycntr_o):
                if((xcntr_o - xcntr_n) < 0):
                    xn = xcntr + abs(xcntr_n - xcntr_o)
                if((ycntr_o - ycntr_n) < 0):
                    yn = ycntr + abs(ycntr_n - ycntr_o)
                if((xcntr_o - xcntr_n) > 0):
                    xn = xcntr - (xcntr_o - xcntr_n)
                if((ycntr_o - ycntr_n) > 0):
                    yn = ycntr - (ycntr_o - ycntr_n)
                tx = x - xn + 0.5 
                ty = y - yn + 0.5
                R = n.sqrt(tx**2.0 + ty**2.0)
                z[n.where(R<=mask_reg*maj_axis)] = 1
                print xn, yn

            elif(abs(xcntr_n - xcntr_o) < size/2.0 - 1.0 and \
                 abs(ycntr_n- ycntr_o) < size/2.0 - 1.0 and \
                 area >= thresh_area):
                if(abs(xcntr_n - xcntr_o) >= major_axis * threshold or \
                   abs(ycntr_n - ycntr_o) >= major_axis):
                    if((xcntr_o - xcntr_n) < 0):
                        xn = xcntr + abs(xcntr_n - xcntr_o)
                    if((ycntr_o - ycntr_n) < 0):
                        yn = ycntr + abs(ycntr_n - ycntr_o)
                    if((xcntr_o - xcntr_n) > 0):
                        xn = xcntr - (xcntr_o - xcntr_n)
                    if((ycntr_o - ycntr_n) > 0):
                        yn = ycntr - (ycntr_o - ycntr_n)
                    print xn, yn
                    tx = x - xn + 0.5 
                    ty = y - yn + 0.5
                    R = n.sqrt(tx**2.0 + ty**2.0)
                    z[n.where(R<=mask_reg*maj_axis)] = 1
        except:
            pass	
    hdu = pyfits.PrimaryHDU(z.astype(n.float32))
    hdu.writeto(mask_file)
