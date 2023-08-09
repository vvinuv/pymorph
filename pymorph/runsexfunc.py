import os
import time
import fitsio
import sys

class PySex(object):
    def __init__(self, SEX_PATH):
        self.SEX_PATH = SEX_PATH



    #def RunSex(self, sex_config, gimg, wimg, sex_cata, 
    #           SEx_GAIN, check_fits=None, sconfig='default'):
    def RunSex(self, sex_config, cutimage, wimg, sex_cata, 
               SEx_GAIN, check_fits=None, sconfig='default'):
        print('cutimage', cutimage)   
        print('sex_cata', sex_cata)
        SEx_DETECT_MINAREA= sex_config['SEx_DETECT_MINAREA']
        SEx_DETECT_THRESH = sex_config['SEx_DETECT_THRESH']
        SEx_ANALYSIS_THRESH = sex_config['SEx_ANALYSIS_THRESH']
        SEx_FILTER = sex_config['SEx_FILTER']
        SEx_FILTER_NAME = sex_config['SEx_FILTER_NAME']
        SEx_DEBLEND_NTHRESH = sex_config['SEx_DEBLEND_NTHRESH']
        SEx_DEBLEND_MINCONT = sex_config['SEx_DEBLEND_MINCONT']
        SEx_PHOT_FLUXFRAC = sex_config['SEx_PHOT_FLUXFRAC']
        SEx_BACK_SIZE = sex_config['SEx_BACK_SIZE']
        SEx_BACK_FILTERSIZE = sex_config['SEx_BACK_FILTERSIZE']
        SEx_BACKPHOTO_TYPE = sex_config['SEx_BACKPHOTO_TYPE']
        SEx_BACKPHOTO_THICK = sex_config['SEx_BACKPHOTO_THICK']
        SEx_WEIGHT_TYPE = sex_config['SEx_WEIGHT_TYPE']
        SEx_PIXEL_SCALE = sex_config['SEx_PIXEL_SCALE']
        SEx_SEEING_FWHM = sex_config['SEx_SEEING_FWHM']
        mag_zero = sex_config['SEx_MAG_ZEROPOINT']

        PYMORPH_PATH = os.path.dirname(__file__)

        #print('Sex 1')

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
            shallow_cata = f'{sex_cata}.Shallow'
            if wimg is None:
                fsex = 'default_wow_shallow.sex'
            else:
                fsex = 'default_shallow.sex'
            print('SExtractor Detecting Objects (Shallow)')
            print(shallow_cata, sex_cata)
        print('fsex', fsex)
        sconfig = os.path.join(PYMORPH_PATH, 'SEx', fsex)

        #print(vars())
        f_tpl = open(sconfig, 'r')
        template = f_tpl.read()
        f_tpl.close()

        f_sex = open(fsex, 'w')
        f_sex.write(template %vars())
        f_sex.close()
     
        #print('Sex 2')

        d, h = fitsio.read(cutimage, header=True)
        #print('Header cutimage', h)
        cmd = '{} {} -c {} > /dev/null'.format(self.SEX_PATH, cutimage, fsex)
        #print(cmd)
        os.system(cmd)
        #print('Sex 2')
        check_fits_not_exists = False
        sleep = 0
        ti = time.time()
        #print('check_fits', check_fits)
        while (check_fits_not_exists) & (sleep < 10):
            if os.path.exists(check_fits):
                check_fits_not_exists = False
            time.sleep(1)
            sleep = time.time() - ti



    def get_sexvalue(self, target_id, sex_cata, param='sky'):

        with open(sex_cata) as f:
            sex_params = f.readlines()
        sex_params = [x.strip().split() for x in sex_params]
        sex_params = np.array(sex_params).astype(float)

        sex_params = sex_params[sex_params[:, 19] == target_id]
        if param == 'sky':
            sky = sex_params[10]
        return sky




