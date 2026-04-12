import pyfits, os
import numpy as n
import config as c
import numpy.ma as ma
from flagfunc import *

class BkgdFunc:
    "The class which will provide the blank sky region and sky deviation to the casgm class"
    def __init__(self, cutimage, xcntr, ycntr, bxcntr, bycntr, eg, pa, sky):
        self.cutimage = cutimage
        self.xcntr = xcntr
        self.ycntr = ycntr
        self.bxcntr = bxcntr
        self.bycntr = bycntr
        self.eg = eg
        self.pa = pa
        self.sky = sky
        self.bkgd = bkgd(cutimage, xcntr, ycntr, bxcntr, bycntr, eg, pa, sky)
        return    

def bkgd(cutimage, xcntr, ycntr, bxcntr, bycntr, eg, pa, sky):
    xcntr = xcntr - 1 # this is because python index starts from 0
    ycntr = ycntr - 1
    angle = c.angle
    back_extraction_radius = c.back_extraction_radius
    f = pyfits.open(c.datadir + cutimage)
    z = f[0].data
    header = f[0].header
    if ('sky' in header):
        sky = header['sky']
    f.close()
    z = n.swapaxes(z, 0, 1)

    nxpts = z.shape[0]
    nypts = z.shape[1]

    if (bxcntr != 9999 and bycntr != 9999):
        back_ini_xcntr = bxcntr 
        back_ini_ycntr = bycntr 
        back_region = z[int(bycntr - back_extraction_radius):int(bycntr + \
                      back_extraction_radius), int(bxcntr -\
                      back_extraction_radius):int(bxcntr + \
                      back_extraction_radius)]
        skysig = back_region.std()
        # Define BackMask here so it exists for the write later
        BackMask = n.where(abs(z - sky) < skysig * 3.0, 0, z)
    else:
        f = pyfits.open("BMask.fits")
        bgmask = f[0].data
        f.close()
        bgmask = n.swapaxes(bgmask, 0, 1)
        bgmaskedgalaxy = ma.masked_array(z, bgmask)
        bgmaskedgalaxy1d = bgmaskedgalaxy.compressed()
        skysig = ma.std(bgmaskedgalaxy1d)
        sky_iter = ma.average(bgmaskedgalaxy1d)
        skysig_iter = skysig * 1.5
        
        x = n.reshape(n.arange(nxpts * nypts), (nxpts, nypts)) / nypts
        x = x.astype(n.float32)
        y = n.reshape(n.arange(nxpts * nypts), (nxpts, nypts)) % nypts
        y = y.astype(n.float32)
        
        countback = 1 
        FLAG_BACK = FLAG_BACK1 = 0
        while FLAG_BACK1 == 0:
            bxcntr = back_extraction_radius
            # Use // for integer division to satisfy range() in Python 3
            for i in range(int((nxpts - int(2 * back_extraction_radius)) // 4)):
                bycntr = back_extraction_radius
                for j in range(int((nypts - int(2 * back_extraction_radius)) // 4)):
                    if(FLAG_BACK1 == 0):
                        tx = x - bxcntr
                        ty = y - bycntr
                        R = n.sqrt(tx**2.0 + ty**2.0)
                        
                        ValueOfRegion = n.abs(z[n.where(R <= back_extraction_radius)] - sky_iter)
                        SizeOfRegion = ValueOfRegion.size
                        SizeOfIgnRegion = ValueOfRegion[n.where(ValueOfRegion > skysig_iter)].size
                        
                        BackMask = n.where(abs(z - sky_iter) < skysig * 3.0, 0, z) 
                        BackMask = n.where(abs(BackMask) > 0, 1, BackMask) 
                        
                        if(SizeOfIgnRegion < 0.2 * SizeOfRegion and FLAG_BACK == 0):
                            FLAG_BACK = 1
                            FLAG_BACK1 = 1
                            back_ini_xcntr = bxcntr
                            back_ini_ycntr = bycntr

                    bycntr += 4.0
                bxcntr += 4.0
            
            skysig_iter *= 1.3
            if countback == 3 and FLAG_BACK1 == 0:
                FLAG_BACK1 = 1
                # Ensure these functions exist in your flagfunc module
                c.Flag = SetFlag(c.Flag, GetFlag('BACK_FAILED'))
            countback += 1

    os.system('rm -f BackMask.fits')
    BackMask = n.swapaxes(BackMask, 0, 1)
    hdu = pyfits.PrimaryHDU(BackMask)
    hdu.writeto('BackMask.fits')
    
    return back_ini_xcntr, back_ini_ycntr, skysig