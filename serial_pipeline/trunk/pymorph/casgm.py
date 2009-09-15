import pyfits,os
import csv
from os.path import exists
import numpy as n
import numpy.ma as ma
import ndimage as im
import config as c
from concfunc import *
from asymfunc import *
from clumfunc import *
from ginifunc_modi import *
#from momentfunc import *
from pyraf import iraf
from runsexfunc import *

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
    xcntr = xcntr #this is because python index statrs from 0
    ycntr = ycntr
#    cmd = 'ds9 ' + str(maskimage)
#    os.system(cmd)
    FoundNewCntr = 0
    if xcntr > 35.0 or ycntr > 35.0:
        dectThre = 18.0
    else:
	dectThre = 12.0
    while FoundNewCntr == 0:
	RunSex(cutimage, 'None', 'CaSsEx.cat', dectThre, dectThre, 1)
        for line in open('CaSsEx.cat', 'r'):
  	    try:
	        values = line.split()
	        if abs(float(values[1]) - xcntr) < 4.001 and abs(float(values[2]) - ycntr) < 4.001:
                    xcntr = float(values[1]) - 1.0
                    ycntr = float(values[2]) - 1.0
		    FoundNewCntr = 1
            except:
	        pass
	os.system('rm -f CaSsEx.cat CaSsEx.cat.sex')
        dectThre -= 2.0
    angle = c.angle
    back_extraction_radius = c.back_extraction_radius
    f = pyfits.open(cutimage)
    z = f[0].data
    header = f[0].header
    if (header.has_key('sky')):
        sky = header['sky']
    f.close()
    z = n.swapaxes(z, 0, 1)
#    print cutimage
    nxpts = z.shape[0]
    nypts = z.shape[1]
    f_err = open('error.log', 'a')
    try:
        print "Initial background Center >>> (", back_ini_xcntr, \
               back_ini_ycntr, ")"
        casgmrun = 1
    except:
        casgmrun = 0
        f_err.writelines(['Failed to find the background region!!!\n'])
    f_err.close()
    z = z - sky
    f=pyfits.open(maskimage)
    mask = f[0].data
    f.close()
    mask = n.swapaxes(mask, 0, 1)
    maskedgalaxy = ma.masked_array(z, mask)
    z = ma.filled(maskedgalaxy, 0.0)
    ########################
    #   CONCENTRATION      #
    ########################
    if(casgmrun):
	try:
	    ApErTuRe = c.aperture
	except:
	    print 'aperture keyword is not found in config.py. Setting '\
	          'circular'
	    ApErTuRe = 1
	if ApErTuRe:
            con=concentration(z, mask, xcntr, ycntr, nxpts, nypts, 0.0, 0.0, sky)
	else:
	    con=concentration(z, mask, xcntr, ycntr, nxpts, nypts, pa - 90.0, eg, sky)
        extraction_radius=con.total_rad
        r20 = con.r20
        r50 = con.r50
        r80 = con.r80
        r90 = con.r90
    else:
        extraction_radius == 9999
    if(extraction_radius == 9999):
        return 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999
    else:
        sigma=0.25*extraction_radius/1.5

        print "R20 R50 R80 R90 Extraction Radius >>> ", str(r20)[:5], \
              str(r50)[:5], str(r80)[:5], str(r90)[:5], str(con.total_rad)[:5]
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
#            asy_r20_zsum = asymmetry(cutimage, maskimage, xcntr, ycntr, 0, 0,\
#                                r50, r20, sky, angle, 0, ABS_ZSUM_r20) This was commented on sep13 as i forgot what this is
            asy_r20_zsum = 0 # This line is added to compensate the above commenting of line. I have replaced the corresponding value to 0 at the line 239
            back_asy = asymmetry(cutimage, maskimage, back_ini_xcntr, \
                                back_ini_ycntr, 0, 0, r50, \
                                back_extraction_radius, \
                                sky, angle, 0, ABS_ZSUM)
            if asy.image_asymm[4] > 20 or back_asy.image_asymm[4] > 20:
                c.Flag += 262144
            if asy.image_asymm[5] == 1:
                c.Flag += 524288 
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
#		print asy.image_asymm[0] ,  back_asy.image_asymm[0]
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
#        print 'sky sigma ', skysig
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
            tmp_writer = csv.writer(f_tmp)
            tmp_ParamToWrite = ['gal_id','C','C_err','A','A_err','A_flag','image_A','back_A','A_20','A_20_with_zsum','S','S_err','r20','r50','r80','r90','extraction_radius','G','G_res','G80','G50','G20','M','M_res','M80','M50','M20']
            tmp_writer.writerow(tmp_ParamToWrite)
#            f_tmp.writelines(['gal_id', '\t', 'C', '\t','C_err', '\t', 'A', '\t', 'A_err', '\t', 'A_flag', '\t', 'image_A', '\t', 'back_A', '\t', 'A_20', '\t', 'A_20_with_zsum', '\t', 'S', '\t', 'S_err', '\t', 'r20', '\t', 'r50', '\t', 'r80', '\t', 'r90', '\t','extraction_radius', '\t', 'G', '\t', 'G_res', '\t', 'G80', '\t', 'G50', '\t', 'G20', '\t', 'M', '\t', 'M_res', '\t', 'M80', '\t', 'M50', '\t', 'M20\n'])
            f_tmp.close()
        f_tmp = open("agm_result_with_radius.csv", "ab")
        tmp_writer = csv.writer(f_tmp)
        tmp_ParamToWrite = [str(cutimage)[to_remove:-5], str(con.concen), str(con.error_con), str(ASY), str(ASY_ERROR), str(asy.image_asymm[5]), str(asy.image_asymm[0]), str(back_asy.image_asymm[0]), str(asy_r20.image_asymm[0]), str(0.0), str(S), str(ERROR_SMOO), str(con.r20), str(con.r50), str(con.r80), str(con.r90), str(extraction_radius), str(gini_coef[0]), str(gini_coef[1]), str(gini_coef[2]), str(gini_coef[3]), str(gini_coef[4]), str(gini_coef[5]), str(gini_coef[6]), str(gini_coef[7]), str(gini_coef[8]), str(gini_coef[9])]
        tmp_writer.writerow(tmp_ParamToWrite)
#        f_tmp.writelines([str(cutimage)[to_remove:-5], '\t', str(con.concen), '\t', str(con.error_con), '\t', str(ASY), '\t', str(ASY_ERROR), '\t',str(asy.image_asymm[5]), '\t',str(asy.image_asymm[0]), '\t',str(back_asy.image_asymm[0]), '\t',str(asy_r20.image_asymm[0]), '\t',str(asy_r20_zsum.image_asymm[0]), '\t', str(S), '\t', str(ERROR_SMOO), '\t', str(con.r20), '\t', str(con.r50), '\t', str(con.r80), '\t', str(con.r90), '\t', str(extraction_radius), '\t', str(gini_coef[0]), '\t', str(gini_coef[1]), '\t', str(gini_coef[2]), '\t', str(gini_coef[3]), '\t', str(gini_coef[4]), '\t', str(gini_coef[5]), '\t', str(gini_coef[6]), '\t', str(gini_coef[7]), '\t', str(gini_coef[8]), '\t', str(gini_coef[9]), '\n'])
        f_tmp.close()
        if str(con.concen) in ('nan', 'inf', '-inf', '-nan'):
            con.concen = 9999
        else:
            pass
        if str(con.error_con) in ('nan', 'inf', '-inf', '-nan'):
            con.error_con = 9999
        else:
            pass
        if str(ASY) in ('nan', 'inf', '-inf', '-nan'):
            ASY = 9999
        else:
            pass
        if str(ASY_ERROR) in ('nan', 'inf', '-inf', '-nan'):
            ASY_ERROR = 9999
        else:
            pass
        if str(S) in ('nan', 'inf', '-inf', '-nan'):
            S = 9999
        else:
            pass
        if str(ERROR_SMOO)  in ('nan', 'inf', '-inf', '-nan'):
            ERROR_SMOO = 9999
        else:
            pass
        if str(gini_coef[0])  in ('nan', 'inf', '-inf', '-nan'):
            Gini_Coef = 9999
        else:
            Gini_Coef = gini_coef[0]
        if str(gini_coef[5])  in ('nan', 'inf', '-inf', '-nan'):
            M20_coef = 9999
        else:
            M20_coef = float(gini_coef[5])
        return con.concen, con.error_con, ASY, ASY_ERROR, S, ERROR_SMOO, Gini_Coef, M20_coef

#CasGm('n5585_lR.fits', 'BMask.fits', 192.03, 157.42, 40.0, 40.0, 0.0, 0.0, 0.0)
