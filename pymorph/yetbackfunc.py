import os
import numpy as np
import pymconvolve
import numpy.ma as ma
import fitsio
from mask_or_fit import GetSExObj
from runsexfunc import PySex

def QuarterMask(z, zm, xcntr, ycntr, bbya, pa, quarter):

    nxpts, nypts = z.shape
    zmm = np.ones_like(z)
    co = np.cos(pa * np.pi / 180.0)
    si = np.sin(pa * np.pi / 180.0)
    one_minus_eg_sq = (bbya)**2.0

    x, y = np.meshgrid(np.arange(nxpts), np.arange(nypts))
    xrot = (x - xcntr) * co + (y - ycntr) * si
    xrot = xrot.T
    xsq = xrot**2.0
    yrot = (xcntr - x) * si + (y - ycntr) * co
    yrot = yrot.T
    ysq = yrot**2.0 
    r = np.sqrt(xsq + ysq / one_minus_eg_sq)

    if quarter == 0: 
        condition = xrot > -1e5
    if quarter == 1:
        condition = (xrot - 0 >= 0) & (yrot - 0 >= 0)
    if quarter == 2:
        condition = (xrot - 0 < 0) & (yrot - 0 >= 0)
    if quarter == 3:
        condition = (xrot - 0 < 0) & (yrot - 0 < 0)
    if quarter == 4:
        condition = (xrot - 0 >= 0) & (yrot - 0 < 0)

    zmm[condition] = 0
    zmm = zm + zmm
    zmm[np.where(zmm > 0)] = 1

    return np.median(ma.masked_array(z, zmm).compressed())

#XXX gimg was here instead of cutimage
def FindYetSky(sex_params, SEX_PATH, cutimage, wimg, seg_file, 
               X0, Y0, scat, SEx_GAIN,
               center_err=5., median_std=1.3, sconfig='seg'):

    #from astropy.io import fits

    ##XXX
    #print('scat', scat)
    PS = PySex(SEX_PATH)
    #PS.RunSex(sex_params, gimg, wimg, scat, SEx_GAIN, sconfig='seg')
    PS.RunSex(sex_params, cutimage, wimg, scat, SEx_GAIN, sconfig='seg')

    #f = fitsio.FITS(gimg)
    f = fitsio.FITS(cutimage)
    z = f[0].read()
    f.close()

    #print(z.shape)
    #print(gimg)
    #print(cutimage)

    fseg = fitsio.FITS(seg_file)
    zm = fseg[0].read()
    fseg.close()

    #f = fits.open(gimg)
    #z = f[0].data
    #f.close()

    #fseg = fits.open(seg_file)
    #zm = fseg[0].data
    #fseg.close()

    #print(zm.shape)

    SexSky, SkyYet = 9999, 9999
    SkyMed, SkyMin = 9999, 9999
    SkyQua, SkySig = 9999, 9999

    
    sex_values = np.genfromtxt(scat, skip_header=0)
    #print(sex_values)
    if sex_values.ndim == 1:
        sex_values = np.expand_dims(sex_values, axis=0)
        
    for values in sex_values:
        #print(values.shape)
        obj = GetSExObj(NXPTS=None, NYPTS=None, values=values)
        SexId = obj.sex_num
        xcntr = obj.xcntr
        ycntr = obj.ycntr
        pa = obj.pos_ang
        bbya = obj.bbya
        a = obj.maj_axis
        b = a * bbya
        hr = obj.radius
        sky = obj.sky

        if np.abs(X0 - obj.xcntr) < center_err and np.abs(Y0 - obj.ycntr) < center_err:
            
            boxcar = np.reshape(np.ones(3 * 3), (3, 3))
            zm = pymconvolve.Convolve(zm, boxcar)
            zm[np.where(zm > 0)] = 1
            
            #print(np.unique(zm, return_counts=True))

            SkyQua = []

            for ii in np.arange(1, 5): 
                SkyQua.append(QuarterMask(z, zm, 
                                         obj.xcntr - 1.0, obj.ycntr - 1.0, 
                                         bbya, pa, ii))

            SkyQua = np.array(SkyQua)
            SexSky = obj.sky

            tmpstd = np.std(ma.masked_array(z, zm).compressed())
            tmpmed = np.median(ma.masked_array(z, zm).compressed())
            zm[np.where((z - tmpmed) > median_std * tmpstd)] = 1
            SkyYet = np.median(ma.masked_array(z, zm).compressed())
            SkyMed = np.median(SkyQua)
            SkyMin = np.min(SkyQua)
            SkySig = np.std(ma.masked_array(z, zm).compressed())
            #               os.system('rm -f SegCat.cat default_seg.sex seg.fits') 
            break

    #print(SexSky, SkyYet, SkyMed, SkyMin, SkyQua, SkySig)
    return SexSky, SkyYet, SkyMed, SkyMin, SkyQua, SkySig
