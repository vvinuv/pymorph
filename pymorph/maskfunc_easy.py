import os
import fitsio
import numpy as np
import pymconvolve
from mask_or_fit import *



class tMaskFunc:
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

    def gmask(self, seg_img, fstring, center_limit=5, seg_limit=1e-5):
        try:
            fseg = fitsio.FITS('/Users/vinu/github/pymorph_refactoring/pymorph/' + seg_img, 'r')
            seg_mask = fseg[0].read()
            fseg.close()
        except:
            print('Problem fits')

        mask_file = 'EM_{}.fits'.format(fstring)

        for line_j in open('/Users/vinu/github/pymorph_refactoring/pymorph/SegCat.cat','r'):
            try:
                values = line_j.split()
                id_n = float(values[0])
                xcntr_n  = float(values[1]) #x center of the neighbour
                ycntr_n  = float(values[2]) #y center of the neighbour

                if np.abs(xcntr_n - self.xcntr_o) < center_limit and np.abs(ycntr_n - self.ycntr_o) < center_limit and self.galflag == 1:
                    seg_mask[np.where(seg_mask == id_n)] = 0
                    pos_ang = float(values[11]) #position angle
                    axis_rat = 1.0/float(values[12]) #axis ration b/a
                    si = np.sin(pos_ang * np.pi / 180.0)
                    co = np.cos(pos_ang * np.pi / 180.0)
                    area = float(values[13])
                    maj_axis = float(values[14])#major axis of neighbour
                    eg = 1.0 - axis_rat
                    one_minus_eg_sq    = (1.0-eg)**2.0
                    break
            except:
                pass

        boxcar = np.reshape(np.ones(3 * 3), (3, 3))
        seg_mask = pymconvolve.Convolve(seg_mask, boxcar)
        seg_mask = np.where(seg_mask > seg_limit, 1., 0.)

        if self.galflag:
            fits = fitsio.FITS(mask_file, 'rw')
            fits.write(seg_mask)
        else:
            fits = fitsio.FITS(mask_file, 'rw')
            fits.write(seg_mask, clobber=True)





class MaskFunc:
    """

    The class for making mask for GALFIT. It uses the masking conditions 
    from config.py. The output mask image will have the name
    M_string(galid).fits 

    """
    def __init__(self, mimg, xcntr, ycntr, NXPTS, NYPTS, line_s):
        self.mimg = mimg
        self.xcntr = xcntr
        self.ycntr = ycntr
        self.NXPTS = NXPTS
        self.NYPTS = NYPTS
        self.line_s  = line_s

    def gmask(self, threshold, thresh_area, 
              avoidme=0, NoMask=False, seg_limit=1e-5):
        
        target = GetSExObj(self.NXPTS, self.NYPTS, line_s)
        
        mask_file = 'M_{}.fits'.format(fstring)

        fseg = fitsio.FITS('seg.fits', 'r')
        seg_mask = fseg.read()
        fseg.close()

        for line_j in open('SegCat.cat','r'):
            if line_j[0] != '#': #line is not a comment
                neighbor = GetSExObj(self.NXPTS, self.NYPTS, line_j)
                
                if target.get_mask(neighbor, threshold, thresh_area, avoidme) != 1:
                    seg_mask[np.where(seg_mask == neighbor.sex_num)] = 0

        if NoMask:
            seg_mask[np.where(seg_mask > 0)] = 0
        else:
            pass

        boxcar = np.reshape(np.ones(3 * 3), (3, 3))
        seg_mask = pymconvolve.Convolve(seg_mask, boxcar)
        seg_mask = np.where(seg_mask > seg_limit, 1, 0)

        fseg = fitsio.FITS(mask_file, 'rw')
        fseg.write(mask_file)
        fseg.close()


