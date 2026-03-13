import os
import fitsio
import numpy as np
from .pymconvolve import pConvolve
from .mask_or_fit import GetSExObj

class MaskFunc:
    """

    The class for making mask for GALFIT. It uses the masking conditions 
    from config.py. The output mask image will have the name
    M_string(galid).fits 

    """
    def __init__(self, DATADIR, fstring):
        mimg = 'M_{}.fits'.format(fstring)
        self.mimg = os.path.join(DATADIR, mimg)

    def gmask(self, threshold, thresh_area, seg_fits, fit_neighbor_cutimage, 
              avoidme=0, NoMask=False, seg_limit=1e-5, verbose=False):
        
        if 1:#verbose:
            print('Maskfunc: mimg ', self.mimg)
        seg_mask = fitsio.read(seg_fits)

        #fit_neighbor_cutimage from configfunc.py. It has coordinates of neighbors and the target to be fitted. Then we need to remove the corresponind IDs of those objects 
        neighbor_id = seg_mask[fit_neighbor_cutimage[:, 1], fit_neighbor_cutimage[:, 0]]
        #print(neighbor_id)
        for nid in neighbor_id:
            seg_mask[seg_mask == nid] = 0
        #print('1, seg_mask')
        #print(seg_mask)

        if NoMask:
            seg_mask[np.where(seg_mask > 0)] = 0
        else:
            pass

        #print('MMMMMM')
        boxcar = np.reshape(np.ones(3 * 3), (3, 3))
        seg_mask = pConvolve(seg_mask, boxcar)
        #print('2, seg_mask')
        #print(seg_mask)
        #sys.exit()
        seg_mask[seg_mask > seg_limit] = 1
        seg_mask[seg_mask != 1] = 0

        #seg_mask = np.where(seg_mask > seg_limit, 1, 0)

        fitsio.write(self.mimg, seg_mask, clobber=True)



class MaskFunc_tmp:
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

    def gmask(self, threshold, thresh_area, seg_fits, seg_cat, 
              avoidme=0, NoMask=False, seg_limit=1e-5, verbose=False):
        
        target = GetSExObj(NXPTS=self.NXPTS, NYPTS=self.NYPTS, 
                           values=self.values)
        print('MMMMMM')
        print(self.values)
        if verbose:
            print('Maskfunc: mimg ', mimg)
        print('seg_fits', seg_fits, seg_cat)
        seg_mask = fitsio.read(seg_fits)

        #print('1, seg_mask')
        #print(seg_mask)
        for line_j in open(seg_cat, 'r'):
            if line_j[0] != '#': #line is not a comment
                obj = []
                for line in line_j.split():
                    obj.append(float(line))
                neighbor = GetSExObj(NXPTS=self.NXPTS, NYPTS=self.NYPTS, 
                                     values=obj)
                
                print(obj)
                if target.get_mask(neighbor, threshold, thresh_area, avoidme) != 1:
                    #print('Yes', neighbor.sex_num)
                    seg_mask[np.where(seg_mask == neighbor.sex_num)] = 0

        if NoMask:
            seg_mask[np.where(seg_mask > 0)] = 0
        else:
            pass

        print('MMMMMM')
        boxcar = np.reshape(np.ones(3 * 3), (3, 3))
        seg_mask = pConvolve(seg_mask, boxcar)
        #print('2, seg_mask')
        #print(seg_mask)
        #sys.exit()
        seg_mask[seg_mask > seg_limit] = 1
        seg_mask[seg_mask != 1] = 0

        #seg_mask = np.where(seg_mask > seg_limit, 1, 0)

        fitsio.write(self.mimg, seg_mask, clobber=True)


