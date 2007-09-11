import pyfits,os,time,sys
from numarray import *
import numarray.ma as ma
import config as c
from concfunc import *
from asymfunc import *
from clumfunc import *
from ginifunc import *
from momentfunc import *
from bgndfunc import *
from pyraf import iraf
import numarray.nd_image as im

class CasGm:
    """The class which will find CASGM parameters"""
    def __init__(self, cutimage, maskimage, xcntr, ycntr, eg, pa, sky):
        self.cutimage   = cutimage
        self.maskimage  = maskimage
        self.xcntr = xcntr
        self.ycntr = ycntr
        self.eg = eg
        self.pa = pa
        self.sky = sky
        self.casgm = casgm(cutimage, maskimage, xcntr, ycntr, eg, pa, sky)
        return

def casgm(cutimage, maskimage, xcntr, ycntr, eg, pa, sky):
    xcntr = xcntr-1 #this is because python index statrs from 0
    ycntr = ycntr-1
    angle = c.angle
    back_extraction_radius = c.back_extraction_radius
    tstart = time.clock() # Start clock to measure time
    f = pyfits.open(cutimage)
    z = f[0].data
    f.close()
    print cutimage
    nxpts = z.shape[0]
    nypts = z.shape[1]

    #The function which will find the blank region 
    try:
        back_ini_xcntr = c.back_ini_xcntr #initial x coordinate of
        back_ini_ycntr = c.back_ini_ycntr #blank portion center
    except:
        f = pyfits.open("BMask.fits")
        bgmask = f[0].data
        f.close()
        bgmaskedgalaxy = ma.masked_array(z, bgmask)
        bgmaskedgalaxy1d = bgmaskedgalaxy.compressed()
        skysig = im.standard_deviation(bgmaskedgalaxy1d)
        sky_iter = ma.average(bgmaskedgalaxy1d)
        skysig_iter = skysig * 0.8
        x = reshape(arange(nxpts * nypts),(nxpts, nypts)) % nypts
        x = x.astype(Float32)
        y = reshape(arange(nxpts * nypts),(nxpts, nypts)) / nypts
        y = y.astype(Float32)
        countback = 1 #After some iteration the following loop quits
        FLAG_BACK = FLAG_BACK1 = 0
        while FLAG_BACK1 == 0:
            bxcntr = back_extraction_radius
            for i in range((nxpts - int(2 * back_extraction_radius)) / 4):
                bycntr = back_extraction_radius
                for j in range((nypts - int(2 * back_extraction_radius)) / 4):
#                print bxcntr, bycntr, FLAG_BACK1, FLAG_BACK
                    if(FLAG_BACK1 == 0):
                        tx = x - bxcntr
                        ty = y - bycntr
                        R = sqrt(tx**2.0 + ty**2.0)
#                        B_AVE = z[where(R <= back_extraction_radius)].mean()
#                        B_STD = z[where(R <= back_extraction_radius)].stddev()
#                        if(B_AVE < sky_iter * 2.5 and B_AVE > sky_iter / 2.5 and sky_iter > 0 and FLAG_BACK == 1 and B_STD < skysig_iter):
#                            FLAG_BACK1 = 1
#                            back_ini_xcntr1 = bxcntr
#                            back_ini_ycntr1 = bycntr
#                        if(B_AVE > sky_iter * 2.5 and B_AVE < sky_iter / 2.5 and sky_iter < 0 and FLAG_BACK == 1 and B_STD < skysig_iter):
#                            FLAG_BACK1 = 1
#                            back_ini_xcntr1 = bxcntr
#                            back_ini_ycntr1 = bycntr
#                        if(B_AVE < sky_iter * 2.5 and B_AVE > sky_iter / 2.5 and sky_iter > 0 and FLAG_BACK == 0 and B_STD < skysig_iter):
#                            FLAG_BACK = 1 
#                            back_ini_xcntr = bxcntr
#                            back_ini_ycntr = bycntr
#                        if(B_AVE > sky_iter * 2.5 and B_AVE < sky_iter / 2.5 and sky_iter < 0 and FLAG_BACK == 0 and B_STD < skysig_iter):
#                            FLAG_BACK = 1
#                            back_ini_xcntr = bxcntr
#                            back_ini_ycntr = bycntr
                        bmax = abs(z[where(R <= back_extraction_radius)] \
                               - sky_iter).max()
                        if(bmax < skysig_iter and FLAG_BACK == 1):
                            FLAG_BACK1 = 1
                            back_ini_xcntr1 = bxcntr
                            back_ini_ycntr1 = bycntr  
                        if(bmax < skysig_iter and FLAG_BACK == 0):
                            FLAG_BACK = 1
                            back_ini_xcntr = bxcntr
                            back_ini_ycntr = bycntr
                    bycntr += 4.0
                bxcntr += 4.0
            skysig_iter *= 1.2
            countback += 1
    z = z - sky
    f=pyfits.open(maskimage)
    mask = f[0].data
    f.close()
    maskedgalaxy = ma.masked_array(z, mask)
    z = ma.filled(maskedgalaxy, value=0.0)
    ########################
    #   CONCENTRATION      #
    ########################
    con=concentration(z, xcntr, ycntr, nxpts, nypts, 0, 0, sky)
    extraction_radius=con.total_rad
    if(extraction_radius == 9999):
        return 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999
    else:
        sigma=0.25*extraction_radius/1.5

        print "EXTRACTION RADIUS ",con.total_rad
        print "CONCENTRATIN AND ERROR ",con.concen,con.error_con

        ########################
        #   ASYMMETRY          #
        ########################

        asy = asymmetry(cutimage, maskimage, xcntr, ycntr, 0, 0, \
              extraction_radius, sky, angle, 1, 0)
        extraction_radius = asy.image_asymm[8]
        ABS_ZSUM = asy.image_asymm[6] * (back_extraction_radius * \
                   back_extraction_radius) / (extraction_radius * \
                   extraction_radius * 1.0)
        back_asy = asymmetry(cutimage, maskimage, back_ini_xcntr, \
                             back_ini_ycntr, 0, 0, back_extraction_radius, \
                             sky, angle, 0, ABS_ZSUM)
        try:
            back_asy1 = asymmetry(cutimage, maskimage, back_ini_xcntr1, \
                       back_ini_ycntr1,\
                       0, 0, back_extraction_radius, sky, angle, 0, ABS_ZSUM)
            ASY = asy.image_asymm[0] - (back_asy.image_asymm[0] +\
                                    back_asy1.image_asymm[0]) / 2.0
            ASY_ERROR = 2 * sqrt(asy.image_asymm[1]**2 + \
                      back_asy.image_asymm[1]**2 + back_asy1.image_asymm[1]**2)
        except:
            ASY = asy.image_asymm[0] - back_asy.image_asymm[0]
            ASY_ERROR = 2 * sqrt(asy.image_asymm[1]**2 \
                        + back_asy.image_asymm[1]**2)
        print "ASYMMETRY, ERROR and flag_out ", \
               ASY, ASY_ERROR, asy.image_asymm[5]
        ########################
        #   CLUMPNESS          #
        ########################
        sigma = int(sigma)
        if(sigma / 2.0 == int(sigma / 2.0)):
            sigma = sigma + 1.0
        clump = clumpness(z, asy.image_asymm[2], asy.image_asymm[3], 0, 0, \
                          extraction_radius, sigma, sky, 1)
        S1 = 10.0 * clump.image_clumpness[0] / clump.image_clumpness[2]
        error_S1 = sqrt((clump.image_clumpness[1] + clump.image_clumpness[3] / \
                         clump.image_clumpness[4]) * S1**2.0)
        if(sigma > back_extraction_radius):
            back_extraction_radius = sigma + 2.0
        back_clump = clumpness(z, back_ini_xcntr, back_ini_ycntr, 0, 0,\
                               back_extraction_radius, sigma, sky, 0)
        S2 = 10.0 * back_clump.image_clumpness[0] / clump.image_clumpness[2]
        error_S2 = sqrt((back_clump.image_clumpness[1] \
                        + clump.image_clumpness[3] \
                        / clump.image_clumpness[4]) * S2**2.0)
        try:
            back_clump1 = clumpness(z, back_ini_xcntr1, back_ini_ycntr1, 0, 0,\
                                     back_extraction_radius, sigma, sky, 0)
            S3 = 10.0 * back_clump1.image_clumpness[0] / \
                 clump.image_clumpness[2]
            error_S3 = sqrt((back_clump1.image_clumpness[1] + \
               clump.image_clumpness[3]  / clump.image_clumpness[4]) * S3**2.0)
            S = S1 - (S2 +S3) / 2.0
            ERROR_SMOO = sqrt(error_S1**2.0 + error_S2**2.0 + error_S3**2.0)
        except:
            S = S1 - S2
            ERROR_SMOO = sqrt(error_S1**2.0 + error_S2**2.0)
        print "SMOTHNESS AND ERROR ", S, ERROR_SMOO

        ########################
        #   GINI COEFFICIENT   #
        ########################

        extraction_radius = con.total_rad
        gin = gini(z, xcntr, ycntr, 0, 0, extraction_radius, sky, skysig)
        gini_coef = gin.gini_coef
        print "GINI COEFFI ",gini_coef
        for myfile in ['segmentation.fits']:
            if os.access(myfile,os.F_OK):
                os.remove(myfile)
        # Write Model galaxy image
        hdu = pyfits.PrimaryHDU(gin.segmentation.astype(Float32))
        hdu.writeto('segmentation.fits')

        ########################
        #   MOMENT CALCULATION #
        ########################
        mo=moment(gin.segmentation, xcntr, ycntr)
        M=mo.moment_of_light[0]
        print "MOMENT ",M
        return con.concen, con.error_con, ASY, ASY_ERROR, S, ERROR_SMOO, gini_coef, M
