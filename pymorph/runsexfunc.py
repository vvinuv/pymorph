import os
import time

#class SexParams_back:

    #    """
    #
    #    The class for running SExtractor. It runs 
    #    
    #    1. if the pipeline doesn't find any SExtractor catalogue. It uses 
    #    the default.* files for doing that. 
    #
    #    2. Needs to find segmentation map
    #
    #    3. Need to find a shallow run
    #
    #    """
    #
    #    def __init__(self, sex_params, mag_zero):
    #
    #        self.SEx_DETECT_MINAREA = sex_params[0]
    #        self.SEx_DETECT_THRESH = sex_params[1]
    #        self.SEx_ANALYSIS_THRESH = sex_params[2]
    #        self.SEx_FILTER = sex_params[3]
    #        self.SEx_FILTER_NAME = sex_params[4]
    #        self.SEx_DEBLEND_NTHRESH = sex_params[5]
    #        self.SEx_DEBLEND_MINCONT = sex_params[6]
    #        self.SEx_PHOT_FLUXFRAC = sex_params[7]
    #        self.SEx_BACK_SIZE = sex_params[8]
    #        self.SEx_BACK_FILTERSIZE = sex_params[9]
    #        self.SEx_BACKPHOTO_TYPE = sex_params[10]
    #        self.SEx_BACKPHOTO_THICK = sex_params[11]
    #        self.SEx_WEIGHT_TYPE = sex_params[12]
    #        self.SEx_PIXEL_SCALE = sex_params[13]
    #        self.SEx_SEEING_FWHM = sex_params[14]
    #        self.mag_zero = mag_zero
    #    
    #
    #    def RunSex(self, SP, gimg, wimg, scat, SEx_GAIN, sconfig='default',
    #               SEx_DETECT_MINAREA = SP.SEx_DETECT_MINAREA, 
    #               SEx_DETECT_THRESH = SP.SEx_DETECT_THRESH,
    #               SEx_ANALYSIS_THRESH = SP.SEx_ANALYSIS_THRESH, 
    #               SEx_FILTER = SP.SEx_FILTER,
    #               SEx_FILTER_NAME = SP.SEx_FILTER_NAME,
    #               SEx_DEBLEND_NTHRESH = SP.SEx_DEBLEND_NTHRESH,
    #               SEx_DEBLEND_MINCONT = SP.SEx_DEBLEND_MINCONT,
    #               SEx_PHOT_FLUXFRAC = SP.SEx_PHOT_FLUXFRAC,
    #               SEx_BACK_SIZE = SP.SEx_BACK_SIZE,
    #               SEx_BACK_FILTERSIZE = SP.SEx_BACK_FILTERSIZE,
    #               SEx_BACKPHOTO_TYPE = SP.SEx_BACKPHOTO_TYPE,
    #               SEx_BACKPHOTO_THICK = SP.SEx_BACKPHOTO_THICK,
    #               SEx_WEIGHT_TYPE = SP.SEx_WEIGHT_TYPE,
    #               SEx_PIXEL_SCALE = SP.SEx_PIXEL_SCALE,
    #               SEx_SEEING_FWHM = SP.SEx_SEEING_FWHM,
    #               mag_zero = SP.mag_zero
    #              ):
    #
    #        PYMORPH_PATH = os.path.dirname(__file__)
    #
    #        if sconfig == 'default':
    #            if wimg is None:
    #                fsex = 'default_wow.sex'
    #            else:
    #                fsex = 'default.sex'
    #            print('SExtractor Detecting Objects (Deep)')
    #        elif sconfig == 'seg':
    #            if wimg is None:
    #                fsex = 'default_seg_wow.sex'
    #            else:
    #                fsex = 'default_seg.sex'
    #            print('SExtractor Detecting segmentation')
    #        elif sconfig == 'shallow':
    #            if wimg is None:
    #                fsex = 'default_wow_shallow.sex'
    #            else:
    #                fsex = 'default_shallow.sex'
    #            print('SExtractor Detecting Objects (Shallow)')
    #
    #        sconfig = os.path.join(PYMORPH_PATH, 'SEx', fsex)
    #
    #        print(vars())
    #        f_tpl = open(sconfig, 'r')
    #        template = f_tpl.read()
    #        f_tpl.close()
    #
    #        f_sex = open(fsex, 'w')
    #        f_sex.write(template %vars())
    #        f_sex.close()
    #
    #        cmd = '{} {} -c {} > /dev/null'.format(SEX_PATH, gimg, fsex)
    #        os.system(cmd)
    #


class PySex(object):
    def __init__(self, SEX_PATH):
        self.SEX_PATH = SEX_PATH

    def RunSex(self, sex_params, SEX_PATH, gimg, wimg, sex_cat, 
               SEx_GAIN, check_fits='check.fits', sconfig='default'):
        

        SEx_DETECT_MINAREA = sex_params[0]
        SEx_DETECT_THRESH = sex_params[1]
        SEx_ANALYSIS_THRESH = sex_params[2]
        SEx_FILTER = sex_params[3]
        SEx_FILTER_NAME = sex_params[4]
        SEx_DEBLEND_NTHRESH = sex_params[5]
        SEx_DEBLEND_MINCONT = sex_params[6]
        SEx_PHOT_FLUXFRAC = sex_params[7]
        SEx_BACK_SIZE = sex_params[8]
        SEx_BACK_FILTERSIZE = sex_params[9]
        SEx_BACKPHOTO_TYPE = sex_params[10]
        SEx_BACKPHOTO_THICK = sex_params[11]
        SEx_WEIGHT_TYPE = sex_params[12]
        SEx_PIXEL_SCALE = sex_params[13]
        SEx_SEEING_FWHM = sex_params[14]
        mag_zero = sex_params[15]

        PYMORPH_PATH = os.path.dirname(__file__)


        if sconfig == 'default':
            if wimg is None:
                fsex = 'default_wow.sex'
            else:
                fsex = 'default.sex'
            print('SExtractor Detecting Objects (Deep)')
        elif sconfig == 'seg':
            if wimg is None:
                fsex = 'default_seg_wow.sex'
            else:
                fsex = 'default_seg.sex'
            print('SExtractor Detecting segmentation')
        elif sconfig == 'shallow':
            if wimg is None:
                fsex = 'default_wow_shallow.sex'
            else:
                fsex = 'default_shallow.sex'
            print('SExtractor Detecting Objects (Shallow)')

        sconfig = os.path.join(PYMORPH_PATH, 'SEx', fsex)

        #print(vars())
        f_tpl = open(sconfig, 'r')
        template = f_tpl.read()
        f_tpl.close()

        f_sex = open(fsex, 'w')
        f_sex.write(template %vars())
        f_sex.close()

        cmd = '{} {} -c {} > /dev/null'.format(SEX_PATH, gimg, fsex)
        print(cmd)
        os.system(cmd)

        check_fits_not_exists = True
        sleep = 0
        ti = time.time()
        while (check_fits_not_exists) & (sleep < 10):
            if os.path.exists(check_fits):
                check_fits_not_exists = False
            time.sleep(1)
            sleep = time.time() - ti 

    def get_sexobj(self):
        with open(filename) as f:
            sex_params = f.readlines()
        sex_params = [x.strip().split() for x in sex_params] 
        sex_params = np.array(sex_params).astype(float)
        xpix = sex_params[:, 1]
        ypix = sex_params[:, 2]
        xarc = sex_params[:, 3]
        yarc = sex_params[:, 4]
            
        if alpha_j == 9999:
            distance = np.sqrt((ximg - xpix)**2 + (yimg - ypix)**2) 
            SexTargets = distance[distance <= dmin].shape[0]
        else:
            distance = np.sqrt((alpha_j - xarc)**2 + (delta_j - yarc)**2)
            SexTargets = distance[distance <= dmin].shape[0]

        self.target = mf.GetSExObj(values=sex_params[np.argmin(distance)])







