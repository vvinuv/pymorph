import os
import csv
from os.path import exists
import pyfits
import numpy as np
import numpy.ma as ma
from concfunc import concentration
from asymfunc import asymmetry
from clumfunc import clumpness 
from ginifunc_modi import gini
from runsexfunc import RunSex
from flagfunc import GetFlag, SetFlag
import pymorphutils as ut
import config as c

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

    def casgm(self):
        # following till END will find better center of image
        FoundNewCntr = 0
        if self.xcntr > 35.0 or self.ycntr > 35.0:
            dectThre = 18.0
        else:
            dectThre = 12.0
        while FoundNewCntr == 0:
            RunSex(sex_config, SEX_PATH, 
                   os.path.join(c.datadir, self.cutimage), 'None', 
                   'CaSsEx.cat', dectThre, dectThre, 1)

            for line in open('CaSsEx.cat', 'r'):
                try:
                    values = line.split()
                    if abs(float(values[1]) - xcntr) < 4.001 and \
                           abs(float(values[2]) - ycntr) < 4.001:
                        xcntr = float(values[1]) - 1.0
                        ycntr = float(values[2]) - 1.0
                    FoundNewCntr = 1
                except:
                    pass

            for myfile in ['CaSsEx.cat', 'CaSsEx.cat.sex']:
                if os.access(myfile, os.F_OK):
                    os.remove(myfile)

            if dectThre < 2.0:
                dectThre -= 0.5
            else:
                dectThre -= 2.0

            if dectThre == 0:
                xcntr = xcntr
                ycntr = ycntr
                FoundNewCntr = 1
        # END

        angle = c.angle
        back_extraction_radius = c.back_extraction_radius
        # open cutimage
        f = pyfits.open(os.path.join(c.datadir, self.cutimage))
        z = f[0].data
        header = f[0].header
        if ('sky' in header):
            #sky = header['sky']
            print('Header has a sky value key. If you want to use that') 
            print('uncomment line number 75 in casgm.py')
        f.close()

        try:
            print("Initial background Center >>> ({}, {})".format( self.back_ini_xcntr, self.back_ini_ycntr))
            casgmrun = 1
        except:
            casgmrun = 0
            ut.WriteError('Failed to find the background region!!!\n')

        z = z - sky
        f = pyfits.open(maskimage)
        mask = f[0].data
        f.close()

        maskedgalaxy = ma.masked_array(z, mask)
        z = ma.filled(maskedgalaxy, 0.0) # filling 0 in mask regions

        ########################
        #   CONCENTRATION      #
        ########################
        if(casgmrun):
            try:
                ApErTuRe = c.aperture
            except:
                print('aperture keyword is not found in config.py. Setting '\
                      'circular')
                ApErTuRe = 1
            if ApErTuRe:
                con = concentration(z, mask, xcntr, ycntr, 0.0, 0.0, sky)
            else:
                con = concentration(z, mask, xcntr, ycntr, pa - 90.0, eg, sky)
            extraction_radius = con.TotRad
            r20 = con.r20
            r50 = con.r50
            r80 = con.r80
            r90 = con.r90
        else:
            extraction_radius == 9999
        if(extraction_radius == 9999):
            return 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999
        else:
            sigma = 0.25 * extraction_radius / 1.5 # The kernal size for 
                                                   # clumpiness

            print('R20: {}'.format(round(r20, 2)))
            print('R50: {}'.format(round(r50, 2))) 
            print('R80: {}'.format(round(r80, 2)))
            print('R90: {}'.format(round(r90, 2)))
            print('Total Radius: {}'.format(round(con.TotRad, 2)))

            
            ########################
            #   ASYMMETRY          #
            ########################
            try:
                asy = asymmetry(cutimage, maskimage, xcntr, ycntr, 0, 0, r50,
                                extraction_radius, sky, angle, 1, 0)
                extraction_radius = asy.image_asymm[8]
                ABS_ZSUM = asy.image_asymm[6] * (back_extraction_radius * \
                           back_extraction_radius) / (extraction_radius * \
                           extraction_radius * 1.0)
                asy_r20 = asymmetry(cutimage, maskimage, xcntr, ycntr, 0, 0, 
                                    r50, r20, sky, angle, 1, 0)
                ABS_ZSUM_r20 = asy.image_asymm[6] * (r20 * r20) / \
                           (extraction_radius * extraction_radius * 1.0)
                # asy_r20_zsum = asymmetry(cutimage, maskimage, xcntr, ycntr, 0, \
                # 0, r50, r20, sky, angle, 0, ABS_ZSUM_r20) This was commented \
                # on sep13 as i forgot what this is
                asy_r20_zsum = 0 # This line is added to compensate the above
                                 # commenting of line. I have replaced the 
                                 # corresponding value to 0 at the line 239
                back_asy = asymmetry(cutimage, maskimage, back_ini_xcntr, 
                                     back_ini_ycntr, 0, 0, r50, 
                                     back_extraction_radius, 
                                     sky, angle, 0, ABS_ZSUM)
                # asymmetry is not converging w. r. t. the center
                if asy.image_asymm[4] > 20 or back_asy.image_asymm[4] > 20:
                    c.Flag = SetFlag(c.Flag, GetFlag('ASYM_NOT_CONV'))
                # the extraction radius is larger than the image size
                if asy.image_asymm[5] == 1:
                    c.Flag = SetFlag(c.Flag, GetFlag('ASYM_OUT_FRAME'))
                try:
                    back_asy1 = asymmetry(cutimage, maskimage, 
                                          back_ini_xcntr1, back_ini_ycntr1, 
                                          0, 0, r50, back_extraction_radius, 
                                          sky, angle, 0, ABS_ZSUM)
                    ASY = asy.image_asymm[0] - (back_asy.image_asymm[0] +\
                                          back_asy1.image_asymm[0]) / 2.0
                    ASY_ERROR = 2 * np.sqrt(asy.image_asymm[1]**2 + \
                         back_asy.image_asymm[1]**2 + back_asy1.image_asymm[1]**2)
                except:
                    ASY = asy.image_asymm[0] - back_asy.image_asymm[0]
                    ASY_ERROR = 2 * np.sqrt(asy.image_asymm[1]**2 \
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
                ASY, ASY_ERROR = 9999, 9999


            ########################
            #   CLUMPNESS          #
            ########################
            try:
                sigma = int(sigma)
                if sigma / 2.0 == int(sigma / 2.0):
                    sigma = sigma + 1.0
                clump = clumpness(z, asy.image_asymm[2], asy.image_asymm[3], 
                                  0, 0, extraction_radius, sigma, sky, 1)
                S1 = 10.0 * clump.clumpness[0] / clump.clumpness[2]
                error_S1 = np.sqrt((clump.clumpness[1] + \
                                 clump.clumpness[3] / \
                                 clump.clumpness[4]) * S1**2.0)
                if sigma > back_extraction_radius:
                    back_extraction_radius = sigma + 2.0
                back_clump = clumpness(z, back_ini_xcntr, back_ini_ycntr, 0, 0,
                                       back_extraction_radius, sigma, sky, 0)
                S2 = 10.0 * back_clump.clumpness[0] / clump.clumpness[2]
                error_S2 = np.sqrt((back_clump.clumpness[1] \
                                + clump.clumpness[3] \
                                / clump.clumpness[4]) * S2**2.0)
                try:
                    back_clump1 = clumpness(z, back_ini_xcntr1, 
                                            back_ini_ycntr1, 0, 0, 
                                            back_extraction_radius, sigma, 
                                            sky, 0)
                    S3 = 10.0 * back_clump1.clumpness[0] / \
                         clump.clumpness[2]
                    error_S3 = np.sqrt((back_clump1.clumpness[1] + \
                               clump.clumpness[3]  / \
                               clump.clumpness[4]) * S3**2.0)
                    S = S1 - (S2 +S3) / 2.0
                    ERROR_SMOO = np.sqrt(error_S1**2.0 + error_S2**2.0 + \
                                 error_S3**2.0)
                except:
                    S = S1 - S2
                    ERROR_SMOO = np.sqrt(error_S1**2.0 + error_S2**2.0)
                try:
                    ERROR_SMOO = round(ERROR_SMOO, 4)
                except:
                    ERROR_SMOO = 9999
    #            print "SMOTHNESS AND ERROR ", S, ERROR_SMOO
            except:
                 S, ERROR_SMOO = 9999, 9999


            ###########################
            #   GINI COEFFICIENT  M20 #
            ###########################

            extraction_radius = con.TotRad # ext. rad was over riden by asym.
            gin = gini(z, xcntr, ycntr, 0, 0, r20, r50, r80, 
                       extraction_radius, sky, skysig)
            gini_coef = gin.segmentation
            # for myfile in ['segmentation.fits']:
            #     if os.access(myfile,os.F_OK):
            #         os.remove(myfile)
            # Write Model galaxy image
            # hdu = pyfits.PrimaryHDU(gin.segmentation.astype(Float32))
            # hdu.writeto('segmentation.fits')
            
            # Writing all the casgm parameters to agm_result_with_radius.csv
            if exists("agm_result_with_radius.csv"):
                pass
            else:
                f_tmp = open("agm_result_with_radius.csv", "ab")
                tmp_writer = csv.writer(f_tmp)
                tmp_ParamToWrite = ['gal_id', 'C', 'C_err', 'A', 'A_err', 
                                    'A_flag', 'image_A', 'back_A', 'A_20', 
                                    'A_20_with_zsum', 'S', 'S_err', 'r20', 
                                    'r20e', 'r50', 'r50e', 'r80', 'r80e', 
                                    'r90', 'r90e', 'extraction_radius', 'G', 
                                    'G_res', 'G80', 'G50', 'G20', 'M', 
                                    'M_res', 'M80', 'M50', 'M20']
                tmp_writer.writerow(tmp_ParamToWrite)
                f_tmp.close()

            f_tmp = open("agm_result_with_radius.csv", "ab")
            tmp_writer = csv.writer(f_tmp)
            tmp_ParamToWrite = [c.fstring, con.concen, con.error_con, 
                                ASY, ASY_ERROR, asy.image_asymm[5], 
                                asy.image_asymm[0], back_asy.image_asymm[0], 
                                asy_r20.image_asymm[0], 0.0, S, ERROR_SMOO, 
                                con.r20, con.r20e, con.r50, con.r50e, con.r80, 
                                con.r80e, con.r90, con.r90e,
                                extraction_radius, gini_coef[0], gini_coef[1],
                                gini_coef[2], gini_coef[3], gini_coef[4], 
                                gini_coef[5], gini_coef[6], gini_coef[7], 
                                gini_coef[8], gini_coef[9]]
            tmp_writer.writerow(tmp_ParamToWrite)
            f_tmp.close()

            # check for nan and inf values and overmasking
            if np.isnan(con.concen) or np.isinf(con.concen):
                con.concen = 9999
            if np.isnan(con.error_con) or np.isinf(con.error_con) or  isinstance(con.error_con, np.ma.core.MaskedConstant):
            #the last test is to see if the error failed because it is overmasked 
                con.error_con = 9999
            if np.isnan(ASY) or np.isinf(ASY):
                ASY = 9999
            if np.isnan(ASY_ERROR) or np.isinf(ASY_ERROR):
                ASY_ERROR = 9999
            if np.isnan(S) or np.isinf(S):
                S = 9999
            if np.isnan(ERROR_SMOO) or np.isinf(ERROR_SMOO):
                ERROR_SMOO = 9999

            if np.isnan(gini_coef[0]) or np.isinf(gini_coef[0]):
               Gini_Coef = 9999
            else:
                Gini_Coef = gini_coef[0]

            if np.isnan(gini_coef[5]) or np.isinf(gini_coef[5]):
                M20_coef = 9999
            else:
                M20_coef = float(gini_coef[5])

            return con.concen, con.error_con, ASY, ASY_ERROR, S, ERROR_SMOO, \
                   Gini_Coef, M20_coef

    #CasGm('n5585_lR.fits', 'BMask.fits', 192.03, 157.42, 40.0, 40.0, 0.0, 0.0, 0.0)
