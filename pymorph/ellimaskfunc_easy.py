import os
import sys
import fitsio
import numpy as np
from .pymconvolve import pConvolve

class ElliMaskFunc:

    """
    
    The class for making mask for ellipse task. This will mask all the 
    neibour objects using the masking conditions in config.py

    """

    def __init__(self, DATADIR, fstring):

        mimg = 'EM_{}.fits'.format(fstring)
        self.mimg = os.path.join(DATADIR, mimg) 

    def emask(self, seg_file, seg_cata, xcntr_o, ycntr_o, 
              center_limit=5, seg_limit=1e-5):

        #from astropy.io import fits
        print('seg file', seg_file)
        #print('seg_cata', seg_cata)
       
        
        #XXX fitsio has some issues. I reported. Till that time I am using
        #astropy.io 
        seg = fitsio.read(seg_file)
        #print('seg', seg)               
        #print(type(seg))
        
        #fseg = fits.open(seg_file)
        #seg = fseg[1].data
        #fseg.close()
        #print('seg', seg)

        
        #print('seg_cata', seg_cata)
        #print(np.unique(seg))

        #Find values of objects        
        values = np.genfromtxt(seg_cata)
        #print(values)
        if values.ndim > 1:
            id_n = values[:, 0].astype(int)
            xcntr_n = values[:, 1]
            ycntr_n = values[:, 2]
        else:
            id_n = values[0].astype(int)
            xcntr_n = values[1]
            ycntr_n = values[2]
            
        #Find distance to all objects wrt target and get the minimum distance
        #to find the target 
        dist_neigh_obj = np.sqrt((xcntr_n - xcntr_o)**2. + (ycntr_n - ycntr_o)**2)
        #print(dist_neigh_obj)
        if values.ndim > 1:
            id_n = id_n[np.argmin(dist_neigh_obj)]
            
        #print('id_n', id_n)
        #print('seg_file', seg_file)
        #print(seg)
        #print('id_n', id_n)    
        
        seg[np.where(seg == id_n)] = 0
        
        #sys.exit()
        boxcar = np.reshape(np.ones(3 * 3), (3, 3))
        seg = pConvolve(seg, boxcar)
        mask = np.where(seg > seg_limit, 1., 0.)

        #print(mask.shape, type(mask), self.mimg)

        
        fitsio.write(self.mimg, mask, clobber=True)
        
    def emask_tmp(self, seg_file, seg_cata, fstring, 
              center_limit=5, seg_limit=1e-5):

        #from astropy.io import fits
        mask_file = 'EM_{}.fits'.format(fstring)
        
        print('seg file', seg_file)
        print('seg_cata', seg_cata)
        
        f = fitsio.FITS(seg_file, 'r')
        seg = f[0].read()
        f.close()
               
        #print(type(seg))
        
        #fseg = fits.open(seg_file)
        #seg = fseg[1].data
        #fseg.close()
        #print('seg', seg)

        
        print('seg_cata', seg_cata)
        print(np.unique(seg))
        for line_j in open(seg_cata, 'r'):
            print(line_j)
            values = line_j.split()
            id_n = float(values[0])
            xcntr_n  = float(values[1]) #x center of the neighbour
            ycntr_n  = float(values[2]) #y center of the neighbour

            if np.abs(xcntr_n - xcntr_o) < center_limit and \
               np.abs(ycntr_n - ycntr_o) < center_limit and \
               self.galflag == 1:
                seg[np.where(seg == id_n)] = 0
                print("Mask")
                break

        boxcar = np.reshape(np.ones(3 * 3), (3, 3))
        seg = pConvolve(seg, boxcar)
        mask = np.where(seg > seg_limit, 1., 0.)

        print(mask.shape, type(mask), mask_file)
        f = fitsio.FITS(mask_file, 'rw')
        if self.galflag:
            f.write(mask)
        else:
            f.write(mask, clobber=True)

        #f.writeto(mask_file, mask, overwrite=True)

        

    #line = '1    193.378    158.284 214.8573569 +56.7789966      3555934     3804.634   8.8786   0.0012     36.075     1433.745 -54.4    1.668    19672    38.968  16  0.00'
    #ElliMaskFunc('n5585_lR.fits', 313, line, 0)

