import os
import config as c
import re

class RunSex:
    """The class for running SExtractor, if the pipeline doesn't find any
       SExtractor catalogue. It uses the default.* files for doing that. """
    def __init__(self, cutimage, whtimage, sex_cata, detect_thr, ana_thr):
        self.cutimage = cutimage
        self.whtimage = whtimage
        self.sex    = sex(cutimage, whtimage, sex_cata, detect_thr, ana_thr)

def sex(cutimage, whtimage, sex_cata, detect_thr, ana_thr):
    if sex_cata =='None':
        sex_cata = c.sex_cata
    mag_zero = c.mag_zero #magnitude zero point
    if detect_thr == 9999:
        SEx_DETECT_MINAREA = c.SEx_DETECT_MINAREA 
    if ana_thr == 9999:
        SEx_DETECT_THRESH = c.SEx_DETECT_THRESH
    SEx_ANALYSIS_THRESH = c.SEx_ANALYSIS_THRESH 
    SEx_FILTER = c.SEx_FILTER 
    SEx_FILTER_NAME = c.SEx_FILTER_NAME 
    SEx_DEBLEND_NTHRESH = c.SEx_DEBLEND_NTHRESH
    SEx_DEBLEND_MINCONT = c.SEx_DEBLEND_MINCONT
    SEx_PHOT_FLUXFRAC = c.SEx_PHOT_FLUXFRAC 
    SEx_GAIN = c.SEx_GAIN
    SEx_PIXEL_SCALE = c.SEx_PIXEL_SCALE
    SEx_SEEING_FWHM = c.SEx_SEEING_FWHM 
    SEx_BACK_SIZE = c.SEx_BACK_SIZE 
    SEx_BACK_FILTERSIZE = c.SEx_BACK_FILTERSIZE 
    SEx_BACKPHOTO_TYPE = c.SEx_BACKPHOTO_TYPE 
    SEx_BACKPHOTO_THICK = c.SEx_BACKPHOTO_THICK 
    SEx_WEIGHT_TYPE = c.SEx_WEIGHT_TYPE
    pymorph_path = c.PYMORPH_PATH
    if(whtimage == 'None'):
        f_tpl = open(str(c.PYMORPH_PATH) + '/SEx/default_wow.sex','r')
    else:
        if SEx_WEIGHT_TYPE == 'DECIDE':
            if re.search("rms", whtimage.lower()):
                SEx_WEIGHT_TYPE = 'MAP_RMS'
            elif re.search("weight", whtimage.lower()):
                SEx_WEIGHT_TYPE = 'MAP_WEIGHT'
            else:
                SEx_WEIGHT_TYPE = 'MAP_RMS'
        f_tpl = open(str(c.PYMORPH_PATH) + '/SEx/default.sex','r')
    template = f_tpl.read()
    f_tpl.close()
    sex_conf = str(sex_cata) + '.sex'
    f_sex = open(sex_conf, 'w')
    f_sex.write(template %vars())
    f_sex.close()
    cmd = str(c.SEX_PATH) + ' ' + str(cutimage) + ' -c ' + str(sex_conf)
    os.system(cmd)    
