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
    """The class which will find CASGM parameters."""
    def __init__(self, cutimage, maskimage, xcntr, ycntr, bxcntr, bycntr, eg, pa, sky, skysig):
        self.cutimage = cutimage
        self.maskimage = maskimage
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
    # Find better center of image
    FoundNewCntr = 0
    if xcntr > 35.0 or ycntr > 35.0:
        dectThre = 18.0
    else:
        dectThre = 12.0

    while FoundNewCntr == 0:
        RunSex(c.datadir + cutimage, 'None', 'CaSsEx.cat', dectThre, dectThre, 1)
        if exists('CaSsEx.cat'):
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
            
        if dectThre <= 0:
            FoundNewCntr = 1

    angle = c.angle
    back_extraction_radius = c.back_extraction_radius
    
    # Open cutimage
    f = pyfits.open(c.datadir + cutimage)
    z = f[0].data
    header = f[0].header
    f.close()

    try:
        print("Initial background Center >>> (", back_ini_xcntr, back_ini_ycntr, ")")
        casgmrun = 1
    except:
        casgmrun = 0
        ut.WriteError('Failed to find the background region!!!\n')

    z = z - sky
    f = pyfits.open(maskimage)
    mask = f[0].data
    f.close()
    maskedgalaxy = ma.masked_array(z, mask)
    z = ma.filled(maskedgalaxy, 0.0) 

    extraction_radius = 9999
    
    ########################
    #   CONCENTRATION      #
    ########################
    if casgmrun:
        try:
            ApErTuRe = c.aperture
        except:
            print('aperture keyword not found in config.py. Setting circular')
            ApErTuRe = 1
            
        if ApErTuRe:
            con = concentration(z, mask, xcntr, ycntr, 0.0, 0.0, sky)
        else:
            con = concentration(z, mask, xcntr, ycntr, pa - 90.0, eg, sky)
        
        extraction_radius = con.TotRad
        r20, r50, r80, r90 = con.r20, con.r50, con.r80, con.r90
    else:
        return 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999

    if extraction_radius == 9999:
        return 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999

    sigma = 0.25 * extraction_radius / 1.5 
    print("R20 R50 R80 R90 Extraction Radius >>> ", str(r20)[:5], \
          str(r50)[:5], str(r80)[:5], str(r90)[:5], str(con.TotRad)[:5])
    
    ########################
    #   ASYMMETRY          #
    ########################
    try:
        asy = asymmetry(cutimage, maskimage, xcntr, ycntr, 0, 0, r50, \
              extraction_radius, sky, angle, 1, 0)
        extraction_radius = asy.image_asymm[8]
        ABS_ZSUM = asy.image_asymm[6] * (back_extraction_radius**2) / (extraction_radius**2)
        
        asy_r20 = asymmetry(cutimage, maskimage, xcntr, ycntr, 0, 0, r50, r20, sky, angle, 1, 0)
        
        back_asy = asymmetry(cutimage, maskimage, back_ini_xcntr, \
                             back_ini_ycntr, 0, 0, r50, \
                             back_extraction_radius, \
                             sky, angle, 0, ABS_ZSUM)

        if asy.image_asymm[4] > 20 or back_asy.image_asymm[4] > 20:
            c.Flag = SetFlag(c.Flag, GetFlag('ASYM_NOT_CONV'))
        if asy.image_asymm[5] == 1:
            c.Flag = SetFlag(c.Flag, GetFlag('ASYM_OUT_FRAME'))

        ASY = asy.image_asymm[0] - back_asy.image_asymm[0]
        ASY_ERROR = 2 * np.sqrt(asy.image_asymm[1]**2 + back_asy.image_asymm[1]**2)
        ASY_ERROR = round(ASY_ERROR, 4)
    except:
        ASY, ASY_ERROR = 9999, 9999

    ########################
    #   CLUMPNESS          #
    ########################
    try:
        sigma = int(sigma)
        if sigma % 2 == 0:
            sigma = sigma + 1.0
        clump = clumpness(z, asy.image_asymm[2], asy.image_asymm[3], 0, 0, extraction_radius, sigma, sky, 1)
        S1 = 10.0 * clump.clumpness[0] / clump.clumpness[2]
        error_S1 = np.sqrt((clump.clumpness[1] + clump.clumpness[3] / clump.clumpness[4]) * S1**2.0)
        
        back_clump = clumpness(z, back_ini_xcntr, back_ini_ycntr, 0, 0, back_extraction_radius, sigma, sky, 0)
        S2 = 10.0 * back_clump.clumpness[0] / clump.clumpness[2]
        error_S2 = np.sqrt((back_clump.clumpness[1] + clump.clumpness[3] / clump.clumpness[4]) * S2**2.0)
        
        S = S1 - S2
        ERROR_SMOO = np.sqrt(error_S1**2.0 + error_S2**2.0)
        ERROR_SMOO = round(ERROR_SMOO, 4)
    except:
        S, ERROR_SMOO = 9999, 9999

    ###########################
    #   GINI COEFFICIENT  M20 #
    ###########################
    extraction_radius = con.TotRad 
    gin = gini(z, xcntr, ycntr, 0, 0, r20, r50, r80, extraction_radius, sky, skysig)
    gini_coef = gin.segmentation
    
    # Write result CSV
    if not exists("agm_result_with_radius.csv"):
        with open("agm_result_with_radius.csv", "w", newline='') as f_tmp:
            tmp_writer = csv.writer(f_tmp)
            tmp_writer.writerow(['gal_id','C','C_err','A','A_err','A_flag', \
                                'image_A','back_A','A_20','A_20_with_zsum', \
                                'S','S_err','r20', 'r20e', 'r50', 'r50e', \
                                'r80', 'r80e', 'r90', 'r90e', \
                                'extraction_radius','G','G_res','G80', \
                                'G50','G20','M','M_res','M80','M50','M20'])

    with open("agm_result_with_radius.csv", "a", newline='') as f_tmp:
        tmp_writer = csv.writer(f_tmp)
        tmp_writer.writerow([c.fstring, con.concen, con.error_con, \
                            ASY, ASY_ERROR, asy.image_asymm[5], \
                            asy.image_asymm[0], back_asy.image_asymm[0], \
                            asy_r20.image_asymm[0], 0.0, S, ERROR_SMOO, \
                            con.r20, con.r20e, con.r50, con.r50e, con.r80, \
                            con.r80e, con.r90, con.r90e,\
                            extraction_radius, gini_coef[0], gini_coef[1],\
                            gini_coef[2], gini_coef[3], gini_coef[4], \
                            gini_coef[5], gini_coef[6], gini_coef[7], \
                            gini_coef[8], gini_coef[9]])

    # Data Validation
    results = [con.concen, con.error_con, ASY, ASY_ERROR, S, ERROR_SMOO]
    final_res = []
    for val in results:
        if np.isnan(val) or np.isinf(val) or isinstance(val, np.ma.core.MaskedConstant):
            final_res.append(9999)
        else:
            final_res.append(val)

    Gini_Coef = 9999 if np.isnan(gini_coef[0]) or np.isinf(gini_coef[0]) else gini_coef[0]
    M20_coef = 9999 if np.isnan(gini_coef[5]) or np.isinf(gini_coef[5]) else float(gini_coef[5])

    return (*final_res, Gini_Coef, M20_coef)