import os
import config as c

class RunSex:
    """The class for running SExtractor, if the pipeline doesn't find any
       SExtractor catalogue. It uses the default.* files for doing that. """
    def __init__(self, cutimage, whtimage):
        self.cutimage = cutimage
        self.whtimage = whtimage
        self.sex    = sex(cutimage, whtimage)

def sex(cutimage, whtimage):
    sex_cata = c.sex_cata
    mag_zero = c.mag_zero #magnitude zero point
    pymorph_path = c.PYMORPH_PATH
    if(whtimage == 'None'):
        f_tpl = open(str(c.PYMORPH_PATH) + '/default_wow.sex','r')
    else:
        f_tpl = open(str(c.PYMORPH_PATH) + '/default.sex','r')
    template = f_tpl.read()
    f_tpl.close()
    sex_conf = str(sex_cata) + '.sex'
    f_sex = open(sex_conf, 'w')
    f_sex.write(template %vars())
    f_sex.close()
    cmd = str(c.SEX_PATH) + ' ' + str(cutimage) + ' -c ' + str(sex_conf)
    os.system(cmd)    
