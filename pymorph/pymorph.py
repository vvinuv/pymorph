
"""PyMorph [Py MOrphological Parameters' Hunter], is a pipeline to find the Morphological parameters of galaxy. Authors: Vinu Vikram , Yogesh Wadadekar, Ajit K. Kembhavi. 2008 Feb, Alan Meert 2010"""

import os
import sys
import time
import csv
import traceback
import argparse
import fitsio
import numpy as np
import re
import configparser
from multiprocessing import Pool

import pymorphutils as ut
from flagfunc import GetFlag, isset, SetFlag

from ellimaskfunc_easy import ElliMaskFunc

from maskfunc_easy import MaskFunc

from configfunc import GalfitConfigFunc
#from configtwostep import ConfigIter
from yetbackfunc import FindYetSky
from plotfunc import PlotFunc
from runsexfunc import RunSex
from writehtmlfunc import WriteHtmlFunc
import psffunc as pf  
import mask_or_fit as mf
from pipeline import Pipeline
        


def run_test(option, opt, value, parser):
    print("Using directory {} for output\n".format(value))
    center_deviated = 0
    starthandle = 0
    FindAndFit()
    main()
    if crashhandler:
        starthandle = 1
        os.system('mv restart.cat CRASH.CAT')
        clus_cat = 'CRASH.CAT' 
        main()
    
    sys.exit(0)
    return


def SExtractorConfEdit():
    SEx_DETECT_MINAREA = input("DETECT_MINAREA (6) >>> ")
    try:
        SEx_DETECT_MINAREA = float(SEx_DETECT_MINAREA)
        SEx_DETECT_MINAREA = int(SEx_DETECT_MINAREA)
    except:
        SEx_DETECT_MINAREA = 6
    SEx_DETECT_THRESH = input('DETECT_THRESH (1.5) >>> ')
    try:
        SEx_DETECT_THRESH = float(SEx_DETECT_THRESH)
    except:
        SEx_DETECT_THRESH = 1.5
    SEx_ANALYSIS_THRESH = input('ANALYSIS_THRESH (1.5) >>> ')
    try:
        SEx_ANALYSIS_THRESH = float(SEx_ANALYSIS_THRESH)
    except:
        SEx_ANALYSIS_THRESH = 1.5
    SEx_FILTER = input('FILTER (Y/N) >>> ')
    while SEx_FILTER != 'Y' and SEx_FILTER != 'N' and SEx_FILTER != '':
        SEx_FILTER = input('FILTER (Y/N) >>> ')
    if len(SEx_FILTER) == 0:
        SEx_FILTER = 'Y'
    else:
        SEx_FILTER = SEx_FILTER
    print('Available options for convolve filter are gauss_1.5_3x3.conv(1) '\
          'gauss_2.0_3x3.conv(2) gauss_2.0_5x5.conv(3) gauss_2.5_5x5.conv(4) '\
          'gauss_3.0_5x5.conv(5) gauss_3.0_7x7.conv(6) gauss_4.0_7x7.conv(7) '\
          'gauss_5.0_9x9.conv(8) default(0)')
    SEx_FILTER_NAME  = input('FILTER_NAME (default.conv) >>> ')
    if len(SEx_FILTER_NAME) == 0 or SEx_FILTER_NAME == '0':
        SEx_FILTER_NAME = 'default.conv'
    elif SEx_FILTER_NAME == '1':
        SEx_FILTER_NAME = 'gauss_1.5_3x3.conv'
    elif SEx_FILTER_NAME == '2':
        SEx_FILTER_NAME = 'gauss_2.0_3x3.conv'
    elif SEx_FILTER_NAME == '3':
        SEx_FILTER_NAME = 'gauss_2.0_5x5.conv'
    elif SEx_FILTER_NAME == '4':
        SEx_FILTER_NAME = 'gauss_2.5_5x5.conv'
    elif SEx_FILTER_NAME == '5':
        SEx_FILTER_NAME = 'gauss_3.0_5x5.conv'
    elif SEx_FILTER_NAME == '6':
        SEx_FILTER_NAME = 'gauss_3.0_7x7.conv'
    elif SEx_FILTER_NAME == '7':
        SEx_FILTER_NAME = 'gauss_4.0_7x7.conv'
    elif SEx_FILTER_NAME == '8':
        SEx_FILTER_NAME = 'gauss_5.0_9x9.conv'
    SEx_DEBLEND_NTHRESH = input('DEBLEND_NTHRESH (32) >>> ')
    try:
        SEx_DEBLEND_NTHRESH = float(SEx_DEBLEND_NTHRESH)
        SEx_DEBLEND_NTHRESH = int(SEx_DEBLEND_NTHRESH)
    except:
        SEx_DEBLEND_NTHRESH = 32
    SEx_DEBLEND_MINCONT = input('DEBLEND_MINCONT (0.005) >>> ')
    try:
        SEx_DEBLEND_MINCONT = float(SEx_DEBLEND_MINCONT)
    except:
        SEx_DEBLEND_MINCONT = 0.005
    SEx_PHOT_FLUXFRAC = input('PHOT_FLUXFRAC (0.5) >>> ')
    try:
        SEx_PHOT_FLUXFRAC = float(SEx_PHOT_FLUXFRAC)
    except:
        SEx_PHOT_FLUXFRAC = 0.5
    SEx_pix_scale_disp = 'PIXEL_SCALE (' + str(SEx_PIXEL_SCALE) + ') >>> '
    SEx_PIXEL_SCALE = input(SEx_pix_scale_disp)
    try:
        SEx_PIXEL_SCALE = float(SEx_PIXEL_SCALE)
    except:
        SEx_PIXEL_SCALE = pixelscale
    SEx_SEEING_FWHM = input('SEEING_FWHM (0.11) >>> ')
    try:
        SEx_SEEING_FWHM = float(SEx_SEEING_FWHM )
    except:
        SEx_SEEING_FWHM = pixelscale * 3.37
    SEx_BACK_SIZE = input('BACK_SIZE (64) >>> ')
    try:
        SEx_BACK_SIZE = float(SEx_BACK_SIZE)
        SEx_BACK_SIZE = int(SEx_BACK_SIZE)
    except:
        SEx_BACK_SIZE = 64
    SEx_BACK_FILTERSIZE = input('BACK_FILTERSIZE (3) >>> ')
    try:
        SEx_BACK_FILTERSIZE = float(SEx_BACK_FILTERSIZE)
        SEx_BACK_FILTERSIZE = int(SEx_BACK_FILTERSIZE)
    except:
        SEx_BACK_FILTERSIZE = 3
    SEx_BACKPHOTO_TYPE = input('BACKPHOTO_TYPE (G)LOBAL/(L)OCAL) >>> ')
    while SEx_BACKPHOTO_TYPE != 'G' and SEx_BACKPHOTO_TYPE != 'L' \
          and SEx_BACKPHOTO_TYPE != '':
        SEx_BACKPHOTO_TYPE = input('BACKPHOTO_TYPE (G)LOBAL/(L)OCAL) >>> ')
    if len(SEx_BACKPHOTO_TYPE) == 0:
        SEx_BACKPHOTO_TYPE = 'GLOBAL'
    elif SEx_BACKPHOTO_TYPE == 'G':
        SEx_BACKPHOTO_TYPE = 'GLOBAL'
    elif SEx_BACKPHOTO_TYPE == 'L':
        SEx_BACKPHOTO_TYPE = 'LOCAL'
    if SEx_BACKPHOTO_TYPE == 'LOCAL':
        SEx_BACKPHOTO_THICK = input('BACKPHOTO_THICK (24) >>> ')
    try:
        SEx_BACKPHOTO_THICK = float(SEx_BACKPHOTO_THICK)
        SEx_BACKPHOTO_THICK = int(SEx_BACKPHOTO_THICK)
    except:
        SEx_BACKPHOTO_THICK = 24
    SEx_WEIGHT_TYPE = input('WEIGHT_TYPE (MAP_RMS) >>> ')
    SEx_WEIGHT_TYPE = SEx_WEIGHT_TYPE

    return [SEx_DETECT_MINAREA, SEx_DETECT_THRESH, SEx_ANALYSIS_THRESH, SEx_FILTER, SEx_FILTER_NAME, SEx_DEBLEND_NTHRESH, SEx_DEBLEND_MINCONT, SEx_PHOT_FLUXFRAC, SEx_BACK_SIZE, SEx_BACK_FILTERSIZE, SEx_BACKPHOTO_TYPE, SEx_BACKPHOTO_THICK, SEx_WEIGHT_TYPE, SEx_PIXEL_SCALE, SEx_SEEING_FWHM]



def SExtractorConfDefault():
    SEx_DETECT_MINAREA = 6
    SEx_DETECT_THRESH = 1.5
    SEx_ANALYSIS_THRESH = 1.5
    SEx_FILTER = 'Y'
    SEx_FILTER_NAME = 'default.conv'
    SEx_DEBLEND_NTHRESH = 32
    SEx_DEBLEND_MINCONT = 0.005
    SEx_PHOT_FLUXFRAC = 0.5
    SEx_BACK_SIZE = 64
    SEx_BACK_FILTERSIZE = 3
    SEx_BACKPHOTO_TYPE = 'GLOBAL'
    SEx_BACKPHOTO_THICK = 24
    SEx_WEIGHT_TYPE = 'MAP_RMS'

    return [SEx_DETECT_MINAREA, SEx_DETECT_THRESH, SEx_ANALYSIS_THRESH, SEx_FILTER, SEx_FILTER_NAME, SEx_DEBLEND_NTHRESH, SEx_DEBLEND_MINCONT, SEx_PHOT_FLUXFRAC, SEx_BACK_SIZE, SEx_BACK_FILTERSIZE, SEx_BACKPHOTO_TYPE, SEx_BACKPHOTO_THICK, SEx_WEIGHT_TYPE]



def run_SExtractorConf():
    try:
        SExtractorConf()
    except:
        raise OptionValueError("failure in SextractorConf()")
    return


def rm_sex_cat(sex_cat):
    if os.path.exists(sex_cat):
        os.remove(sex_cat)
    else:
        pass
    return





class PyMorph:
    
    __version__ = 3.0

    def __init__(self, config_file='config.ini'):

        self.InitilizeParams(config_file) 



    def InitilizeParams(self, config_file):

        '''
        Reading configuration file
        '''

        c = configparser.ConfigParser()
        c.read(config_file)
        
        self.DATADIR = c.get('imagecata', 'datadir')
        self.OUTDIR = c.get('imagecata', 'outdir')

        self.imagefile = c.get('imagecata', 'imagefile')
        self.imagefile = os.path.join(self.DATADIR, self.imagefile)

        self.whtfile = c.get('imagecata', 'whtfile')
        self.whtfile = os.path.join(self.DATADIR, self.whtfile)

        self.sex_cat = c.get('imagecata', 'sex_cat')
        self.sex_cat = os.path.join(self.DATADIR, self.sex_cat)

        self.clus_cat = c.get('imagecata', 'clus_cat')
        self.clus_cat = os.path.join(self.DATADIR, self.clus_cat)

        self.out_cat = c.get('imagecata', 'out_cat')
        self.out_cat = os.path.join(self.OUTDIR, self.out_cat)


        self.rootname = c.get('imagecata', 'rootname')

        self.GALFIT_PATH = c.get('external', 'GALFIT_PATH')
        self.SEX_PATH = c.get('external', 'SEX_PATH')
        self.PYMORPH_PATH = c.get('external', 'PYMORPH_PATH')

        self.repeat = c.getboolean('modes', 'repeat')                       
        self.galcut = c.getboolean('modes', 'galcut')                        
        self.decompose = c.getboolean('modes', 'decompose')
        self.detail = c.getboolean('modes', 'detail')
        self.galfit = c.getboolean('modes', 'galfit')
        self.cas = c.getboolean('modes', 'cas')
        self.crashhandler = c.getboolean('modes', 'crashhandler')
        

        try:
            self.components = c.get('galfit', 'components').split(',')
        except:
            print("components undefined. Asuming bulge+disk model")
            self.components = ['bulge', 'disk']

        self.devauc = c.get('galfit', 'devauc')
        self.fitting = c.get('galfit', 'fitting')
        self.fitting = [int(tf) for tf in self.fitting.split(',')]
        
        self.chi2sq = c.getfloat('diagnosis', 'chi2sq')
        self.goodness = c.getfloat('diagnosis', 'goodness')
        self.center_deviation = c.getfloat('diagnosis', 'center_deviation')
        self.center_constrain = c.getfloat('diagnosis', 'center_constrain')
     
       
        try:
            size = c.get('size', 'size_list')
            size = [int(s) for s in size.split(',')]
            self.ReSize = size[0]
            self.VarSize = size[1]
            self.FracRad = size[2]
            self.Square = size[3]
            self.FixSize = size[4] 
        except:
            self.ReSize = c.getint('size', 'size_list')
      
            if self.ReSize:
                self.VarSize = 1
            else:
                self.VarSize = 0
            self.Square = 1
            self.FracRad = 20
            self.FixSize = 120
        try:
            self.searchrad = c.get('size', 'searchrad')
        except:
            self.searchrad = None


        self.pixelscale = c.getfloat('cosmology', 'pixelscale')
        self.H0 = c.getfloat('cosmology', 'H0')
        self.WM = c.getfloat('cosmology', 'WM')
        self.WV = c.getfloat('cosmology', 'WV')
        self.redshift = c.getfloat('cosmology', 'redshift')
 
        self.psfselect = c.getint('psf', 'psfselect')
        self.stargal_prob = c.getfloat('psf', 'stargal_prob')
        self.star_size = c.getint('psf', 'star_size')
        self.psflist = c.get('psf', 'psflist')
        self.which_psf = c.getint('psf', 'which_psf')
        self.area_obj = c.getint('psf', 'which_psf')

        self.manual_mask = c.getint('mask', 'manual_mask')
        self.mask_reg = c.getfloat('mask', 'mask_reg')
        self.thresh_area = c.getfloat('mask', 'thresh_area')
        self.threshold = c.getfloat('mask', 'threshold')

        self.mag_zero = c.getfloat('mag', 'mag_zero')
        self.maglim = c.get('mag', 'maglim').split(',')
        self.maglim = [float(mlim) for mlim in self.maglim]


        self.host = c.get('db', 'host')
        self.database = c.get('db', 'database')
        self.table = c.get('db', 'table')
        self.usr = c.get('db', 'usr')
        self.pword = c.get('db', 'pword')
        self.dbparams = c.get('db', 'dbparams')

        self.galfitv = c.get('version', 'galfitv')

        self.FirstCreateDB = 1 #Won't create table when FirstCreateDB=0
        self.FILTER = 'UNKNOWN'

        self.SEx_PIXEL_SCALE = self.pixelscale
        self.SEx_SEEING_FWHM = self.pixelscale * 3.37
        self.SEx_MAG_ZEROPOINT = self.mag_zero


    def set_sexparams(self, 
                      SEx_DETECT_MINAREA = 6, 
                      SEx_DETECT_THRESH = 1.5,
                      SEx_ANALYSIS_THRESH =1.5,
                      SEx_FILTER = 'Y',        
                      SEx_FILTER_NAME = 'default.conv',
                      SEx_DEBLEND_NTHRESH = 32,
                      SEx_DEBLEND_MINCONT = 0.005,
                      SEx_PHOT_FLUXFRAC = 0.5, 
                      SEx_BACK_SIZE = 64,      
                      SEx_BACK_FILTERSIZE = 3,
                      SEx_BACKPHOTO_TYPE = 'GLOBAL',
                      SEx_BACKPHOTO_THICK = 24,
                      SEx_WEIGHT_TYPE = 'MAP_RMS',     
                      ):
        '''

        Setting SExtractor default.sex

        FHWH is 3.37 * SEx_PIXEL_SCALE

        '''

        self.SEx_DETECT_MINAREA =  SEx_DETECT_MINAREA
        self.SEx_DETECT_THRESH =   SEx_DETECT_THRESH 
        self.SEx_ANALYSIS_THRESH = SEx_ANALYSIS_THRESH
        self.SEx_FILTER =          SEx_FILTER  
        self.SEx_FILTER_NAME =     SEx_FILTER_NAME 
        self.SEx_DEBLEND_NTHRESH = SEx_DEBLEND_NTHRESH
        self.SEx_DEBLEND_MINCONT = SEx_DEBLEND_MINCONT
        self.SEx_PHOT_FLUXFRAC =   SEx_PHOT_FLUXFRAC 
        self.SEx_BACK_SIZE =       SEx_BACK_SIZE 
        self.SEx_BACK_FILTERSIZE = SEx_BACK_FILTERSIZE
        self.SEx_BACKPHOTO_TYPE =  SEx_BACKPHOTO_TYPE
        self.SEx_BACKPHOTO_THICK = SEx_BACKPHOTO_THICK
        self.SEx_WEIGHT_TYPE =     SEx_WEIGHT_TYPE

        self.sex_params = [SEx_DETECT_MINAREA, SEx_DETECT_THRESH, 
                           SEx_ANALYSIS_THRESH, SEx_FILTER, 
                           SEx_FILTER_NAME, SEx_DEBLEND_NTHRESH, 
                           SEx_DEBLEND_MINCONT, SEx_PHOT_FLUXFRAC, 
                           SEx_BACK_SIZE, SEx_BACK_FILTERSIZE, 
                           SEx_BACKPHOTO_TYPE, SEx_BACKPHOTO_THICK, 
                           SEx_WEIGHT_TYPE]

    def set_limits(self,
                   LMag=500, UMag=-500, LN=0.1, UN=20., 
                   LRe=0., URe=500., LRd=0., URd=500., 
                   bdbox=False, bbox=False, dbox=False, devauc=False, 
                   avoidme=50., FILTER='UNKNOWN', 
                   which_psf=0, stargal_prob=0.9, area_obj=40., 
                   no_mask=False, norm_mask=False):
        '''

        Setting pymorph parameters limits

        '''

        self.LMag  = LMag        
        self.UMag  = UMag  
        self.LN    = LN    
        self.UN    = UN    
        self.LRe   = LRe   
        self.URe   = URe   
        self.LRd   = LRd   
        self.URd   = URd   
        self.bdbox = bdbox 
        self.bbox  = bbox  
        self.dbox  = dbox  
        self.devauc = devauc
                     
        self.avoidme = avoidme
        self.FILTER = FILTER
        self.which_psf = which_psf
        self.stargal_prob = stargal_prob
        self.area_obj = area_obj 
        self.no_mask = no_mask
        self.norm_mask = norm_mask



        
        def argparse_input():

            '''

            The below is the idea to get arguments through parser. However, it
            should not be used. Use either set_sexparams() and set_limits()

            Note I set defaults here and call them in creating the options. 
            This is so I can use them in determining whether to retain 
            the default or not.

            '''
            
            
            usage = "Usage: pymorph [--edit-conf[-e]] [--with-psf [-p]]" \
                    "[--force[-f]] [--help[-h]] [--lmag] [--umag]"\
                    "[--ln] [--un]"\
                    "[--lre] [--ure] [--lrd] [--urd] [--with-in] "\
                    "[--with-filter] [--with-db] [--with-area]  [--no-mask]"\
                    "[--norm-mask] [--with-sg] [--bdbox] [--bbox] [--dbox]"\
                    "[--test [-t]] [--devauc] [--outdir] [--datadir]"


            parser = argparse.ArgumentParser(prog='PROG', usage=usage)
            parser.add_argument("-e", "--sex_conf", action="store_const",
                              const=SExtractorConfEdit, 
                              default=SExtractorConfDefault,
                              help="runs SExtractor configuration")
            parser.add_argument("-f", "--force", action="store_const",
                              const=rm_sex_cat,
                              help="removes SExtractor catalog")
            parser.add_argument("-t", "--test", action="store_const", 
                              const=run_test,  
                              help="runs the test instance-OVERRIDES ALL                                    OTHER INPUT-User must supply a directory                                     for output")
            parser.add_argument("--lmag", action="store", 
                              dest="LMag",default = 500.0,
                              help="lower magnitude cutoff")
            parser.add_argument("--umag", action="store", 
                              dest="UMag", default = -500.0,
                              help="upper magnitude cutoff")
            parser.add_argument("--ln", action="store", 
                              dest="LN", default = 0.1, help="Lower Sersic")
            parser.add_argument("--un", action="store", 
                              dest="UN", default = 20.0, help="Upper Sersic")
            parser.add_argument("--lre", action="store", 
                              dest="LRe", default = 0.0,
                              help="Lower Bulge Radius")
            parser.add_argument("--ure", action="store", 
                              dest="URe", default = 500.0,
                              help="Upper Bulge Radius")
            parser.add_argument("--lrd", action="store", 
                              dest="LRd", default = 0.0, help="Lower Disk \
                                                               Radius")
            parser.add_argument("--urd", action="store", 
                              dest="URd", default = 500.0, help="Upper Disk \
                                                                 Radius")

            parser.add_argument("--bdbox", action="store_true", dest="bdbox",
                              default = False, help="turns on bdbox")
            parser.add_argument("--bbox", action="store_true", dest="bbox",
                              default = False, help="turns on bbox")
            parser.add_argument("--dbox", action="store_true", dest="dbox",
                              default = False, help="turns on dbox")
            parser.add_argument("--devauc", action="store_true", dest="devauc",
                              default = False,
                              help="turns on DeVacouleur's bulge fitting\
                                      (sersic index of bulge frozen at 4.0\
                                      for fit)")


            parser.add_argument("--with-in", action="store", 
                              dest="avoidme", default = 50.0, help="avoid me!")
            parser.add_argument("--with-filter", action="store", 
                              dest="FILTER", default = 'UNKNOWN', 
                              help="Filter used")
            
            parser.add_argument("-p", "--with-psf", action="store",
                              dest="which_psf", default = 0,
                              help="Nearest/farthest PSF")
            parser.add_argument("--with-sg", action="store", 
                              dest="stargal_prob", default = 0.9,
                              help="for psf identification")
            parser.add_argument("--with-area", action="store", 
                              dest="area_obj", default = 40.0,
                              help="min Area of psf for selection")
            parser.add_argument("--no_mask", action="store_true", 
                              dest="no_mask",
                              default=False, help="turns off masking")
            parser.add_argument("--norm_mask", action="store_true", 
                              dest="norm_mask",
                              default=False, help="turns on Normal masking")
            


            parser.add_argument("--with-host", action="store", 
                               dest="host", default = 'localhost',
                              help="mysql host used")
            parser.add_argument("--with-db", action="store", 
                               dest="database", default = 'UNKNOWN',
                              help="database used")
            # parses command line aguments for pymorph
            #print('B3')
            args = parser.parse_args()

            #print('B4')
            self.sex_params = args.sex_conf()
            if len(self.sex_params) == 15:
                self.sex_params.append(self.mag_zero)
            else:
                self.sex_params.append(self.pixelscale)
                self.sex_params.append(self.pixelscale * 3.37)        
                self.sex_params.append(self.mag_zero)

            #print('B5')


            self.SEx_DETECT_MINAREA = self.sex_params[0]
            self.SEx_DETECT_THRESH = self.sex_params[1]
            self.SEx_ANALYSIS_THRESH = self.sex_params[2]
            self.SEx_FILTER = self.sex_params[3]
            self.SEx_FILTER_NAME = self.sex_params[4]
            self.SEx_DEBLEND_NTHRESH = self.sex_params[5]
            self.SEx_DEBLEND_MINCONT = self.sex_params[6]
            self.SEx_PHOT_FLUXFRAC = self.sex_params[7]
            self.SEx_BACK_SIZE = self.sex_params[8]
            self.SEx_BACK_FILTERSIZE = self.sex_params[9]
            self.SEx_BACKPHOTO_TYPE = self.sex_params[10]
            self.SEx_BACKPHOTO_THICK = self.sex_params[11]
            self.SEx_WEIGHT_TYPE = self.sex_params[12]
            self.SEx_PIXEL_SCALE = self.sex_params[13]
            self.SEx_SEEING_FWHM = self.sex_params[14]

            args = parser.parse_args()

            self.LMag  = args.LMag
            self.UMag  = args.UMag
            self.LN    = args.LN
            self.UN    = args.UN
            self.LRe   = args.LRe
            self.URe   = args.URe
            self.LRd   = args.LRd
            self.URd   = args.URd
            self.bdbox = args.bdbox
            self.bbox  = args.bbox
            self.dbox  = args.dbox
            self.devauc= args.devauc

            self.avoidme  = args.avoidme
            self.FILTER   = args.FILTER
            self.which_psf = args.which_psf
            self.stargal_prob = args.stargal_prob
            self.area_obj = args.area_obj
            self.no_mask  = args.no_mask
            self.norm_mask= args.norm_mask

            #self.DATADIR  = args.DATADIR
            #self.OUTDIR   = args.OUTDIR

            self.host     = args.host
            self.database = args.database

            #self.SP = SexParams(self.sex_params, self.mag_zero)






    def main_thread(self):

        #print(self.pixelscale, self.ReSize, self.FixSize)

        os.system('rm -f TmpElliMask.fits TmpElliMask1.fits')

        

        #Initialize index.html
        if os.path.exists('index.html'):
            pass
        else:
            indexfile = open('index.html', 'w')
            indexfile.writelines(['<HTML>\n<BODY>\n'])
            indexfile.writelines(['</BODY></HTML>'])
            indexfile.close()


        #Reading image and weigh files
        if self.repeat == False & self.galcut == False:
            print("Using large image. imagefile >>> {}".format(self.imagefile))

            if os.path.exists(self.whtfile):
                wht = fitsio.FITS(self.whtfile)

                if re.search("rms", self.whtfile.lower()):
                    self.P.weightdata = wht[0].read()
                    print("whtfile >>> {}".format(self.whtfile))
                elif re.search("weight", self.whtfile.lower()):
                    self.P.weightdata = 1 / np.sqrt(wht[0].read())
                    print("whtfile >>> {}".format(self.whtfile))
                else:
                    print('Weight file is not understood. Please include ' + \
                          'the word weight/rms to the weight file name. ' + \
                          'If it is weight the rms will be found by 1/sqrt(w)') 
                wht.close()
                self.P.weightexists = True
            else:
                self.P.weightexists = False
                print('No weight image found\n')
        
        #print('B7')
        #Initializing psf array. ie. creating psflist from file 
        self.psflist = pf.PSFArr(self.DATADIR, self.psflist)
        #print(self.psflist)
        #Updating the header RA and DEC of PSF if the filename is in a
        #specific format
        if self.decompose:
            for p in self.psflist:
                #print('psf', self.DATADIR, p)
                pf.UpdatePsfRaDec(self.DATADIR, p)

        self.P.psflist = self.psflist

        #print('B8')
        # Writing csv header and finding the output parameters
        params_to_write = ut.PyMorphOutputParams(self.dbparams, 
                                                 self.decompose)    
        #print(params_to_write)
        self.P.params_to_write = params_to_write
        if os.path.exists('result.csv'):
            pass
        else:
            f_res = open("result.csv", "a")
            csvhead = ['{}_{:d}'.format(params_to_write[par_key][0], par_key)
                       for par_key in params_to_write.keys()]
            writer = csv.writer(f_res)
            writer.writerow(csvhead)
            f_res.close()

        self.f_cat = open(self.out_cat, 'w')
        self.f_failed = open('restart.cat', 'w')
        #obj_file = open(self.clus_cat, 'r')  # The file contains the 
                                                     # objects of interest
        #pnames = obj_file.readline().split() #The names of the parameters given 
                                             #in the first line in the clus_cat

        with open(self.clus_cat, 'r') as f:
            obj_file = f.read().splitlines()

        # writing a input catalogue (restart.cat) for failed objects
        #for pname in pnames:
        pnames = obj_file.pop(0).split()
        for pname in pnames:
            self.f_failed.writelines(['{} '.format(pname)])
        self.f_failed.writelines(['flag \n'])

        ##print(pnames)
        ##print(obj_file)
        self.P.pnames = pnames 
        self.P.psfcounter = 0
        self.psfcounter = 0                  #For getting psf in the case of unknown ra
        #print('Before multi')
        pool = Pool(processes=1)
        results = pool.map(self.P.main, obj_file)
        pool.close()
        pool.join()

        print(results)



    def FindAndFit(self):
        fc = open(self.clus_cat, 'w')
        if self.redshift == 9999:
            fc.writelines(['gal_id ra dec mag\n'])
        else:
            fc.writelines(['gal_id ra dec mag z\n'])

        for line in open(self.sex_cat, 'r'):

            values = line.split()
            if 1:
                if float(values[17]) > self.UMag and \
                   float(values[17]) < self.LMag and \
                   float(values[16]) < self.stargal_prob:
                    if self.redshift == 9999:
                        cline = '{} {} {} {}\n'.format(values[0], 
                                                    float(values[3]) / 15.0,
                                                    values[4], 
                                                    values[17])
                    else:
                        cline = '{} {} {} {} {}\n'.format(values[0], 
                                                      float(values[3]) / 15.0,
                                                      values[4], 
                                                      values[17],
                                                      self.redshift)
                    fc.writelines([cline])
            #except:
            #    pass
        fc.close()

        searchrad = '0.05arc'


    def pymorph(self):

        self.P = Pipeline([self.DATADIR, self.OUTDIR, self.imagefile, 
                           self.whtfile, self.sex_cat, self.clus_cat, 
                           self.out_cat, self.rootname, self.GALFIT_PATH, 
                           self.SEX_PATH, self.repeat, self.galcut, 
                           self.decompose, self.detail, self.galfit,
                           self.cas, self.crashhandler, self.components,
                           self.devauc, self.fitting, self.chi2sq, 
                           self.goodness, self.center_deviation, 
                           self.center_constrain, self.ReSize,
                           self.VarSize, self.FracRad, self.Square, 
                           self.FixSize, 
                           self.searchrad, self.pixelscale, self.H0, self.WM, 
                           self.WV, self.redshift, self.psfselect, 
                           self.stargal_prob,
                           self.star_size, self.psflist, self.which_psf, 
                           self.area_obj, self.manual_mask, self.mask_reg, 
                           self.thresh_area, self.threshold, self.mag_zero,
                           self.maglim, self.host, self.database, self.table,
                           self.usr, self.pword, self.dbparams, self.galfitv, 
                           self.FirstCreateDB, self.FILTER, 
                           self.LMag, self.UMag, 
                           self.LN, self.UN, self.LRe, self.URe, 
                           self.LRd, self.URd, 
                           self.bdbox, self.bbox, self.dbox, self.devauc, 
                           self.avoidme, self.no_mask, self.norm_mask])
                           
        self.P.sex_params = self.sex_params 

        # now change dir to the 
        thisdir = os.getcwd()
        print("Current directory is ", thisdir)
        print("Output directory is ", self.OUTDIR)

        os.chdir(self.OUTDIR)

        try:
            if self.repeat == False & self.galcut == False:
                fimg = fitsio.FITS(self.imagefile)
                self.imagedata = fimg[0].read()
                self.header0 = fimg[0].read_header()
                fimg.close()
                
                gheader = ut.CheckHeader(self.header0)
                self.P.EXPTIME = gheader[0]
                self.P.RDNOISE = gheader[1],
                self.P.GAIN = gheader[2]
                self.SEx_GAIN = gheader[3]
                self.P.SEx_GAIN = gheader[3]
                self.P.NCOMBINE = gheader[4]
            else:
                self.SEx_GAIN = 1.
                self.P.SEx_GAIN = 1.

        except IOError as e:
            print(self.imagefile, "I/O error ({}): {}".format(e.args[0], e.args[1]))
            os._exit(0)
        
        ##print(self.sex_cat)
        self.P.imagedata = self.imagedata
        self.P.center_deviated = 0

        if os.path.exists(self.sex_cat):
            pass
        elif self.galcut == False:
            print('The SExtractor catalogue for your frame is NOT found. ' \
                  'One is being made using the default values. It is always '\
                  'recommended to make SExtractor catalogue by YOURSELF as '\
                  'the pipeline keeps the sky value at the SExtractor value '\
                  'during the decomposition.')
            RunSex(self.sex_params, self.SEX_PATH, 
                   self.imagefile, self.whtfile,
                   self.sex_cat, self.SEx_GAIN,
                   check_fits='check.fits', sconfig='default')

    #The old function for psfselect = 2
    #    if psfselect == 2:
    #        center_deviated = 0
    #        starthandle = 0
    #        os.system('ds9 &')
    #        time.sleep(2)
    #        run_psfselect(imagefile, DATADIR, clus_cat, self.galcut)
    #        os.system('xpaset -p ds9 quit')
    #        psflist = '@psflist.list'
    #        FindAndFit()
    #        main()
    #        if crashhandler:
    #            starthandle = 1
    #            os.system('mv restart.cat CRASH.CAT')
    #            clus_cat = 'CRASH.CAT' 
    #            main()
    #New function for psfselect=2 non-interactive for webservice
        ##print(self.psfselect)
        ##print(self.clus_cat)
        if self.psfselect == 2:
            self.Interactive = 0
            self.starthandle = 0
            run_psfselect(self.imagefile, 
                          self.DATADIR, 
                          self.clus_cat, 
                          self.galcut)
            psflist = os.path.join(self.DATADIR, '@psflist.list')
            self.FindAndFit()
            self.main_thread()
            if self.crashhandler:
                starthandle = 1
                os.system('mv restart.cat CRASH.CAT')
                self.clus_cat = 'CRASH.CAT' 
                self.main_thread()

    #The old function for psfselect=1
        elif self.psfselect == 1:
            self.Interactive = 1
            os.system('ds9 &')
            time.sleep(2)
            run_psfselect(self.imagefile, 
                          self.DATADIR, 
                          self.clus_cat, 
                          self.galcut)
            os.system('xpaset -p ds9 quit')
    #new function for psfselect=1 non-interactive for webservice (abhishek rawat)
    #    elif psfselect == 1:
    #        Interactive = 0
    #        run_psfselect(imagefile, DATADIR, clus_cat, self.galcut)
            



        elif self.psfselect == 0:
            self.starthandle = 0
            #self.FindAndFit()
            #print(self.sex_cat)
            print('B6')
            self.main_thread()
            print('B7')
            if self.crashhandler:
                self.starthandle = 1
                os.system('mv restart.cat CRASH.CAT')
                self.clus_cat = 'CRASH.CAT' 
                self.main_thread()

        os.chdir(thisdir)





if __name__ == '__main__':
    print('Beginning1')
    p = PyMorph()
    p.set_sexparams()
    p.set_limits()
    print('Beginning2')
    #p.pymorph()
    #p.main()
