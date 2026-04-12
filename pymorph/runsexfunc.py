import os
import config as c
import re

class RunSex:
    """
    The class for running SExtractor if the pipeline doesn't find any
    SExtractor catalogue. It uses the default.* files for doing that.
    """
    def __init__(self, cutimage, whtimage, sex_cata, detect_thr, ana_thr, cas):
        self.cutimage = cutimage
        self.whtimage = whtimage
        self.sex_cata = sex_cata
        self.detect_thr = detect_thr
        self.ana_thr = ana_thr
        # Trigger the sex function
        self.sex = sex(cutimage, whtimage, sex_cata, detect_thr, ana_thr, cas)

def sex(cutimage, whtimage, sex_cata, detect_thr, ana_thr, cas):
    if sex_cata == 'None':
        sex_cata = c.sex_cata
        
    # Standard SExtractor parameters from config
    mag_zero = c.mag_zero 
    SEx_DETECT_MINAREA = c.SEx_DETECT_MINAREA
    
    if detect_thr == 9999:
        SEx_DETECT_THRESH = c.SEx_DETECT_THRESH
    else:
        SEx_DETECT_THRESH = detect_thr
        
    if ana_thr == 9999:
        SEx_ANALYSIS_THRESH = c.SEx_ANALYSIS_THRESH 
    else:
        SEx_ANALYSIS_THRESH = ana_thr
        
    SEx_FILTER = c.SEx_FILTER
    SEx_FILTER_NAME = c.SEx_FILTER_NAME 
    SEx_DEBLEND_NTHRESH = c.SEx_DEBLEND_NTHRESH
    SEx_DEBLEND_MINCONT = c.SEx_DEBLEND_MINCONT
    SEx_PHOT_FLUXFRAC = c.SEx_PHOT_FLUXFRAC 
    
    # Python 3 logic: GAIN must be handled carefully
    if cas:
        SEx_GAIN = 1
    else:
        SEx_GAIN = c.SEx_GAIN
        
    SEx_PIXEL_SCALE = c.SEx_PIXEL_SCALE
    SEx_SEEING_FWHM = c.SEx_SEEING_FWHM 
    SEx_BACK_SIZE = c.SEx_BACK_SIZE 
    SEx_BACK_FILTERSIZE = c.SEx_BACK_FILTERSIZE 
    SEx_BACKPHOTO_TYPE = c.SEx_BACKPHOTO_TYPE 
    SEx_BACKPHOTO_THICK = c.SEx_BACKPHOTO_THICK 
    SEx_WEIGHT_TYPE = c.SEx_WEIGHT_TYPE
    pymorph_path = c.PYMORPH_PATH

    # Decide which template to open
    base_wht_name = whtimage.split('/')[-1]
    if 'None' in base_wht_name:
        template_path = os.path.join(str(c.PYMORPH_PATH), 'SEx/default_wow.sex')
    else:
        if SEx_WEIGHT_TYPE == 'DECIDE':
            if re.search("rms", whtimage.lower()):
                SEx_WEIGHT_TYPE = 'MAP_RMS'
            elif re.search("weight", whtimage.lower()):
                SEx_WEIGHT_TYPE = 'MAP_WEIGHT'
            else:
                SEx_WEIGHT_TYPE = 'MAP_RMS'
        template_path = os.path.join(str(c.PYMORPH_PATH), 'SEx/default.sex')

    # Read the template
    if os.path.exists(template_path):
        with open(template_path, 'r') as f_tpl:
            template = f_tpl.read()
    else:
        print(f"Error: SExtractor template not found at {template_path}")
        return

    # Write the specific .sex configuration file
    sex_conf = str(sex_cata) + '.sex'
    with open(sex_conf, 'w') as f_sex:
        # %vars() maps local variables into the template placeholders
        f_sex.write(template % vars())

    print('SExtractor Detecting Objects (Deep)....')
    cmd = f"{str(c.SEX_PATH)} {str(cutimage)} -c {str(sex_conf)} > /dev/null 2>&1"
    os.system(cmd)    

def SexShallow(cutimage, whtimage, sex_cata, detect_thr, ana_thr, cas):
    if sex_cata == 'None':
        sex_cata = f"{c.sex_cata}.Shallow"
        
    DeepCata = c.sex_cata
    mag_zero = c.mag_zero 
    SEx_WEIGHT_TYPE = c.SEx_WEIGHT_TYPE
    pymorph_path = c.PYMORPH_PATH

    base_wht_name = whtimage.split('/')[-1]
    if 'None' in base_wht_name:
        template_path = os.path.join(str(c.PYMORPH_PATH), 'SEx/default_wow_shallow.sex')
    else:
        if SEx_WEIGHT_TYPE == 'DECIDE':
            if re.search("rms", whtimage.lower()):
                SEx_WEIGHT_TYPE = 'MAP_RMS'
            elif re.search("weight", whtimage.lower()):
                SEx_WEIGHT_TYPE = 'MAP_WEIGHT'
            else:
                SEx_WEIGHT_TYPE = 'MAP_RMS'
        template_path = os.path.join(str(c.PYMORPH_PATH), 'SEx/default_shallow.sex')

    if os.path.exists(template_path):
        with open(template_path, 'r') as f_tpl:
            template = f_tpl.read()
    else:
        print(f"Error: Shallow template not found at {template_path}")
        return

    sex_conf = f"{str(sex_cata)}.sex"
    with open(sex_conf, 'w') as f_sex:
        f_sex.write(template % vars())

    print('SExtractor Detecting Objects (Shallow)....')
    cmd = f"{str(c.SEX_PATH)} {str(cutimage)} -c {str(sex_conf)} > /dev/null 2>&1"
    os.system(cmd)