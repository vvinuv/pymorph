import fitsio
import numpy as np
import pymconvolve 
#from astropy.io import fits

class ElliMaskFunc:

    """
    
    The class for making mask for ellipse task. This will mask all the 
    neibour objects using the masking conditions in config.py

    """

    def __init__(self, xcntr_o, ycntr_o, 
                 center_limit=5., seg_limit=1e-5):

        self.xcntr_o = xcntr_o
        self.ycntr_o = ycntr_o

    def emask(self, seg_file, seg_cata, fstring, 
              center_limit=5, seg_limit=1e-5):

        #from astropy.io import fits
        mask_file = 'EM_{}.fits'.format(fstring)
        
        #print('seg file', seg_file)
        #print('seg_cata', seg_cata)
        
        f = fitsio.FITS(seg_file, 'r')
        seg = f[0].read()
        f.close()
               
        #print(type(seg))
        
        #fseg = fits.open(seg_file)
        #seg = fseg[1].data
        #fseg.close()
        #print('seg', seg)

        
        #print('seg_cata', seg_cata)
        print(np.unique(seg))
        
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
            
        dist_neigh_obj = np.sqrt((xcntr_n - self.xcntr_o)**2. + (ycntr_n - self.ycntr_o)**2)
        #print(dist_neigh_obj)
        if values.ndim > 1:
            id_n = id_n[np.argmin(dist_neigh_obj)]
            
        #print('id_n', id_n)
        seg[np.where(seg == id_n)] = 0

        boxcar = np.reshape(np.ones(3 * 3), (3, 3))
        seg = pymconvolve.Convolve(seg, boxcar)
        mask = np.where(seg > seg_limit, 1., 0.)

        #print(mask.shape, type(mask), mask_file)

        
        f = fitsio.FITS(mask_file, 'rw')
        f.write(mask, overwrite=True)
        
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

            if np.abs(xcntr_n - self.xcntr_o) < center_limit and \
               np.abs(ycntr_n - self.ycntr_o) < center_limit and \
               self.galflag == 1:
                seg[np.where(seg == id_n)] = 0
                print("Mask")
                break

        boxcar = np.reshape(np.ones(3 * 3), (3, 3))
        seg = pymconvolve.Convolve(seg, boxcar)
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

