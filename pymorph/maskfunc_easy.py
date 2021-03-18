import os
import fitsio
import numpy as np
import pymconvolve
from mask_or_fit import *




class MaskFunc:
    """

    The class for making mask for GALFIT. It uses the masking conditions 
    from config.py. The output mask image will have the name
    M_string(galid).fits 

    """
    def __init__(self, mimg, xcntr, ycntr, NXPTS, NYPTS, values):

        self.mimg = mimg
        self.xcntr = xcntr
        self.ycntr = ycntr
        self.NXPTS = NXPTS
        self.NYPTS = NYPTS
        self.values  = values

    def gmask(self, threshold, thresh_area, fstring, seg_fits, seg_cat, 
              avoidme=0, NoMask=False, seg_limit=1e-5, verbose=False):
        
        target = GetSExObj(NXPTS=self.NXPTS, NYPTS=self.NYPTS, 
                           values=self.values)
        
        if verbose:
            print('Maskfunc: mimg ', mimg)

        fseg = fitsio.FITS(seg_fits, 'r')
        seg_mask = fseg[0].read()
        fseg.close()

        #print('1, seg_mask')
        #print(seg_mask)
        for line_j in open(seg_cat, 'r'):
            if line_j[0] != '#': #line is not a comment
                obj = []
                for line in line_j.split():
                    obj.append(float(line))
                neighbor = GetSExObj(NXPTS=self.NXPTS, NYPTS=self.NYPTS, 
                                     values=obj)
                
                if target.get_mask(neighbor, threshold, thresh_area, avoidme) != 1:
                    #print('Yes', neighbor.sex_num)
                    seg_mask[np.where(seg_mask == neighbor.sex_num)] = 0

        if NoMask:
            seg_mask[np.where(seg_mask > 0)] = 0
        else:
            pass

        boxcar = np.reshape(np.ones(3 * 3), (3, 3))
        seg_mask = pymconvolve.Convolve(seg_mask, boxcar)
        #print('2, seg_mask')
        #print(seg_mask)
        #sys.exit()
        seg_mask[seg_mask > seg_limit] = 1
        seg_mask[seg_mask != 1] = 0

        #seg_mask = np.where(seg_mask > seg_limit, 1, 0)

        fseg = fitsio.FITS(self.mimg, 'rw')
        fseg.write(seg_mask, clobber=True)
        fseg.close()


