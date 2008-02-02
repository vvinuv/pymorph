import pyfits,os,time,sys
from os.path import exists
import numpy as n
import numpy.core.ma as ma
import ndimage as im
import config as c
from concfunc import *
from asymfunc import *
from clumfunc import *
from ginifunc_modi import *
#from momentfunc import *
from pyraf import iraf

class CasGm:
    """The class which will find CASGM parameters. The algorithm for each 
       parameters can be found in the corresponding class files. This will 
       also write a file agm_result_with_radius.csv which gives the Asymmetry
       Gini Coefficient and M20 parameters at different radius from the center
       of the galaxy"""
    def __init__(self, cutimage, maskimage, xcntr, ycntr, bxcntr, bycntr, eg, pa, sky, skysig):
        self.cutimage   = cutimage
        self.maskimage  = maskimage
        self.xcntr = xcntr
        self.ycntr = ycntr
        self.bxcntr = bxcntr
        self.bycntr = bycntr
        self.eg = eg
        self.pa = pa
        self.sky = sky
        self.skysig = skysig
        self.casgm = casgm(cutimage, maskimage, xcntr, ycntr, bxcntr, bycntr, eg, pa, sky, skysig)
        return

def casgm(cutimage, maskimage, xcntr, ycntr, back_ini_xcntr, back_ini_ycntr, eg, pa, sky, skysig):
    xcntr = xcntr-1 #this is because python index statrs from 0
    ycntr = ycntr-1
    angle = c.angle
    back_extraction_radius = c.back_extraction_radius
    tstart = time.clock() # Start clock to measure time
    f = pyfits.open(cutimage)
    z = f[0].data
    header = f[0].header
    if (header.has_key('sky')):
        sky = header['sky']
    f.close()
#    print cutimage
    nxpts = z.shape[0]
    nypts = z.shape[1]
    f_err = open('error.log', 'a')
    try:
        print "back_ini_xcntr, back_ini_ycntr", back_ini_xcntr, back_ini_ycntr
        casgmrun = 1
    except:
        casgmrun = 0
        f_err.writelines(['Failed to find the background region\n'])
    f_err.close()
    z = z - sky
    f=pyfits.open(maskimage)
    mask = f[0].data
    f.close()
    maskedgalaxy = ma.masked_array(z, mask)
    z = ma.filled(maskedgalaxy, value=0.0)
    ########################
    #   CONCENTRATION      #
    ########################
    if(casgmrun):
        con=concentration(z, xcntr, ycntr, nxpts, nypts, 0, 0, sky)
        extraction_radius=con.total_rad
        r20 = con.r20
        r50 = con.r50
        r80 = con.r80
    else:
        extraction_radius == 9999
    if(extraction_radius == 9999):
        return 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999
    else:
        sigma=0.25*extraction_radius/1.5

        print "r20, r50, r80, EXTRACTION RADIUS ",r20, r50, r80, con.total_rad
#        print "CONCENTRATIN AND ERROR ", con.concen,con.error_con
        
        ########################
        #   ASYMMETRY          #
        ########################
        try:
            asy = asymmetry(cutimage, maskimage, xcntr, ycntr, 0, 0, r50, \
                  extraction_radius, sky, angle, 1, 0)
            extraction_radius = asy.image_asymm[8]
            ABS_ZSUM = asy.image_asymm[6] * (back_extraction_radius * \
                       back_extraction_radius) / (extraction_radius * \
                       extraction_radius * 1.0)
            asy_r20 = asymmetry(cutimage, maskimage, xcntr, ycntr, 0, 0, r50,\
                                r20, sky, angle, 1, 0)
            ABS_ZSUM_r20 = asy.image_asymm[6] * (r20 * r20) / \
                       (extraction_radius * extraction_radius * 1.0)
            asy_r20_zsum = asymmetry(cutimage, maskimage, xcntr, ycntr, 0, 0,\
                                r50, r20, sky, angle, 0, ABS_ZSUM_r20)

            back_asy = asymmetry(cutimage, maskimage, back_ini_xcntr, \
                                back_ini_ycntr, 0, 0, r50, \
                                back_extraction_radius, \
                                sky, angle, 0, ABS_ZSUM)
            try:
                back_asy1 = asymmetry(cutimage, maskimage, back_ini_xcntr1, \
                            back_ini_ycntr1, 0, 0, r50,\
                            back_extraction_radius, sky, angle, 0, ABS_ZSUM)
                ASY = asy.image_asymm[0] - (back_asy.image_asymm[0] +\
                                      back_asy1.image_asymm[0]) / 2.0
                ASY_ERROR = 2 * n.sqrt(asy.image_asymm[1]**2 + \
                       back_asy.image_asymm[1]**2 + back_asy1.image_asymm[1]**2)
            except:
                ASY = asy.image_asymm[0] - back_asy.image_asymm[0]
                ASY_ERROR = 2 * n.sqrt(asy.image_asymm[1]**2 \
                            + back_asy.image_asymm[1]**2)
            try:
                ASY_ERROR = round(ASY_ERROR, 4)    
            except:
                ASY_ERROR = 9999
#            print "ASYMMETRY, ERROR and flag_out A_wo_back BA A20 A20_ZSUM", \
#                   ASY, ASY_ERROR, asy.image_asymm[5], asy.image_asymm[0], \
#                   back_asy.image_asymm[0], asy_r20.image_asymm[0],\
#                   asy_r20_zsum.image_asymm[0]
        except:
            ASY = ASY_ERROR = 9999 
        ########################
        #   CLUMPNESS          #
        ########################
        try:
            sigma = int(sigma)
            if(sigma / 2.0 == int(sigma / 2.0)):
                sigma = sigma + 1.0
            clump = clumpness(z, asy.image_asymm[2], asy.image_asymm[3], 0, 0, \
                              extraction_radius, sigma, sky, 1)
            S1 = 10.0 * clump.image_clumpness[0] / clump.image_clumpness[2]
            error_S1 = n.sqrt((clump.image_clumpness[1] + \
                             clump.image_clumpness[3] / \
                             clump.image_clumpness[4]) * S1**2.0)
            if(sigma > back_extraction_radius):
                back_extraction_radius = sigma + 2.0
            back_clump = clumpness(z, back_ini_xcntr, back_ini_ycntr, 0, 0,\
                                   back_extraction_radius, sigma, sky, 0)
            S2 = 10.0 * back_clump.image_clumpness[0] / clump.image_clumpness[2]
            error_S2 = n.sqrt((back_clump.image_clumpness[1] \
                            + clump.image_clumpness[3] \
                            / clump.image_clumpness[4]) * S2**2.0)
            try:
                back_clump1 = clumpness(z, back_ini_xcntr1, back_ini_ycntr1, \
                              0, 0, back_extraction_radius, sigma, sky, 0)
                S3 = 10.0 * back_clump1.image_clumpness[0] / \
                     clump.image_clumpness[2]
                error_S3 = n.sqrt((back_clump1.image_clumpness[1] + \
                           clump.image_clumpness[3]  / \
                           clump.image_clumpness[4]) * S3**2.0)
                S = S1 - (S2 +S3) / 2.0
                ERROR_SMOO = n.sqrt(error_S1**2.0 + error_S2**2.0 + error_S3**2.0)
            except:
                S = S1 - S2
                ERROR_SMOO = n.sqrt(error_S1**2.0 + error_S2**2.0)
            try:
                ERROR_SMOO = round(ERROR_SMOO, 4)
            except:
                ERROR_SMOO = 9999
#            print "SMOTHNESS AND ERROR ", S, ERROR_SMOO
        except:
             S = ERROR_SMOO = 9999

        ########################
        #   GINI COEFFICIENT   #
        ########################

        extraction_radius = con.total_rad
        print 'sky sigma ', skysig
        gin = gini(z, xcntr, ycntr, 0, 0, r20, r50, r80, extraction_radius, sky, skysig)
#        gini_coef = gin.gini_coef
	gini_coef = gin.segmentation
#        print "GINI COEFFI ",gini_coef
#        for myfile in ['segmentation.fits']:
#            if os.access(myfile,os.F_OK):
#                os.remove(myfile)
        # Write Model galaxy image
#        hdu = pyfits.PrimaryHDU(gin.segmentation.astype(Float32))
#        hdu.writeto('segmentation.fits')

        ########################
        #   MOMENT CALCULATION #
        ########################
#        mo=moment(gin.segmentation, xcntr, ycntr)
#        M=mo.moment_of_light[0]
#        print "MOMENT ",M
        to_remove = len(c.rootname) + 2
        if exists("agm_result_with_radius.csv"):
            pass
        else:
            f_tmp = open("agm_result_with_radius.csv", "ab")
            f_tmp.writelines(['gal_id', '\t', 'C_err', '\t', 'A', '\t', 'A_err', '\t', 'A_flag', '\t', 'image_A', '\t', 'back_A', '\t', 'A_20', '\t', 'A_20_with_zsum', '\t', 'S', '\t', 'S_err', '\t', 'r20', '\t', 'r50', '\t', 'r80', '\t', 'extraction_radius', '\t', 'G', '\t', 'G_res', '\t', 'G80', '\t', 'G50', '\t', 'G20', '\t', 'M', '\t', 'M_res', '\t', 'M80', '\t', 'M50', '\t', 'M20\n'])
            f_tmp.close()
        f_tmp = open("agm_result_with_radius.csv", "ab")
        f_tmp.writelines([str(cutimage)[to_remove:-5], '\t', str(con.concen), '\t', str(con.error_con), '\t', str(ASY), '\t', str(ASY_ERROR), '\t',str(asy.image_asymm[5]), '\t',str(asy.image_asymm[0]), '\t',str(back_asy.image_asymm[0]), '\t',str(asy_r20.image_asymm[0]), '\t',str(asy_r20_zsum.image_asymm[0]), '\t', str(S), '\t', str(ERROR_SMOO), '\t', str(con.r20), '\t', str(con.r50), '\t', str(con.r80), '\t', str(extraction_radius), '\t', str(gini_coef[0]), '\t', str(gini_coef[1]), '\t', str(gini_coef[2]), '\t', str(gini_coef[3]), '\t', str(gini_coef[4]), '\t', str(gini_coef[5]), '\t', str(gini_coef[6]), '\t', str(gini_coef[7]), '\t', str(gini_coef[8]), '\t', str(gini_coef[9]), '\n'])
        f_tmp.close()

        return con.concen, con.error_con, ASY, ASY_ERROR, S, ERROR_SMOO, gini_coef[0], gini_coef[5]

#CasGm('n5585_lR.fits', 'BMask.fits', 192.03, 157.42, 40.0, 40.0, 0.0, 0.0, 0.0)
