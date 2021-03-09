
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

#import config as c
import pymorphutils as ut
from flagfunc import GetFlag, isset, SetFlag

from ellimaskfunc_easy import ElliMaskFunc
#from ellimaskfunc import ElliMaskFunc

from maskfunc_easy import MaskFunc
#from maskfunc import MaskFunc

from configfunc import GalfitConfigFunc
#from configtwostep import ConfigIter
from yetbackfunc import FindYetSky
#from plotfunc import PlotFunc
from runsexfunc import RunSex
#from writehtmlfunc import WriteParams
import psffunc as pf  
import mask_or_fit as mf
       
        


def run_test(option, opt, value, parser):
    print("Using directory {} for output\n".format(value))
    center_deviated = 0
    starthandle = 0
    FindAndFit()
    main()
    if crashhandler:
        starthandle = 1
        os.system('mv restart.cat CRASH.CAT')
        clus_cata = 'CRASH.CAT' 
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


def rm_sex_cata(sex_cata):
    if os.path.exists(sex_cata):
        os.remove(sex_cata)
    else:
        pass
    return

class PyMorph:
    
    __version__ = 3.0

    def __init__(self, config_file='config.ini'):

        self.InitilizeParams(config_file) 


    def InitilizeParams(self, config_file):

        c = configparser.ConfigParser()
        c.read(config_file)
        
        self.DATADIR = c.get('imagecata', 'datadir')
        self.OUTDIR = c.get('imagecata', 'outdir')

        self.imagefile = c.get('imagecata', 'imagefile')
        self.imagefile = os.path.join(self.DATADIR, self.imagefile)

        self.whtfile = c.get('imagecata', 'whtfile')
        self.whtfile = os.path.join(self.DATADIR, self.whtfile)

        self.sex_cata = c.get('imagecata', 'sex_cata')
        self.sex_cata = os.path.join(self.DATADIR, self.sex_cata)

        self.clus_cata = c.get('imagecata', 'clus_cata')
        self.clus_cata = os.path.join(self.DATADIR, self.clus_cata)

        self.out_cata = c.get('imagecata', 'out_cata')
        self.out_cata = os.path.join(self.OUTDIR, self.out_cata)


        self.rootname = c.get('imagecata', 'rootname')

        self.GALFIT_PATH = c.get('external', 'GALFIT_PATH')
        self.SEX_PATH = c.get('external', 'SEX_PATH')

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
        self.Goodness = c.getfloat('diagnosis', 'Goodness')
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

        # Note I set defaults here and call them in creating the options. 
        #This is so I can use them in determining whether to retain 
        #the default or not.

        usage = "Usage: pymorph [--edit-conf[-e]] [--with-psf [-p]]" \
                "[--force[-f]] [--help[-h]] [--lmag] [--umag] [--ln] [--un]"\
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
                          const=rm_sex_cata,
                          help="removes SExtractor catalog")
        parser.add_argument("-t", "--test", action="store_const", 
                          const=run_test,  
                          help="runs the test instance-OVERRIDES ALL OTHER INPUT-User must supply a directory for output")
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
                          dest="LRd", default = 0.0, help="Lower Disk Radius")
        parser.add_argument("--urd", action="store", 
                          dest="URd", default = 500.0, help="Upper Disk Radius")

        parser.add_argument("--bdbox", action="store_true", dest="bdbox",
                          default = False, help="turns on bdbox")
        parser.add_argument("--bbox", action="store_true", dest="bbox",
                          default = False, help="turns on bbox")
        parser.add_argument("--dbox", action="store_true", dest="dbox",
                          default = False, help="turns on dbox")
        parser.add_argument("--devauc", action="store_true", dest="devauc",
                          default = False,
                          help="turns on DeVacouleur's bulge fitting (sersic index of bulge frozen at 4.0 for fit)")


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
        parser.add_argument("--no_mask", action="store_true", dest="no_mask",
                          default=False, help="turns off masking")
        parser.add_argument("--norm_mask", action="store_true", dest="norm_mask",
                          default=False, help="turns on Normal masking")
        


        #parser.add_argument("-d","--datadir", action="store", dest="DATADIR",
        #                   default = os.getcwd()+'/',
        #                  help="path to directory containing all input. MUST end in '/'")
        #parser.add_argument("-o","--outdir", action="store", 
        #                  dest="OUTDIR", default = os.getcwd()+'/',
        #                  help="path to directory that will contain all output. MUST end in '/'")
        

        parser.add_argument("--with-host", action="store", 
                           dest="host", default = 'localhost',
                          help="mysql host used")
        parser.add_argument("--with-db", action="store", 
                           dest="database", default = 'UNKNOWN',
                          help="database used")
        # parses command line aguments for pymorph
        args = parser.parse_args()

        self.sex_params = args.sex_conf()
        if len(self.sex_params) == 15:
            self.sex_params.append(self.mag_zero)
        else:
            self.sex_params.append(self.pixelscale)
            self.sex_params.append(self.pixelscale * 3.37)        
            self.sex_params.append(self.mag_zero)



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



    def FindCutSize(self, SizeXX, SizeYY, SexHalfRad, pa, axis_rat):
        """Return the size of the cutout. SizeXX, SizeYY are the size of 
          cutout if c.galcut is True"""

        pa = pa * np.pi / 180. #pa in radian
        SizeX = SexHalfRad * self.FracRad * np.abs(np.cos(pa)) + \
                axis_rat * SexHalfRad * self.FracRad * np.abs(np.sin(pa))
        SizeY = SexHalfRad * self.FracRad * np.abs(np.sin(pa)) + \
                axis_rat * SexHalfRad * self.FracRad * np.abs(np.cos(pa))
        SizeX = int(SizeX)
        SizeY = int(SizeY)

        if self.Square:
            SizeX = max(SizeX, SizeY)
            SizeY = max(SizeX, SizeY)
        if self.galcut:
            if self.ReSize:
                if self.VarSize:
                    pass
                else:
                    SizeX = self.FixSize
                    SizeY = self.FixSize
            else:
                # FIX
                SizeX = SizeXX
                SizeY = SizeYY
                # END
        else:
            if self.VarSize:
                pass
            else:
                SizeX = self.FixSize
                SizeY = self.FixSize
        return SizeX, SizeY

    def distance(self, psffile, ra, dec):
        """
        Find the distance between psf and object in arcsec. 
        Ra and dec is the position of the object. 
        Psf coordinates will be read from the header

        """
        
        try:
            p = fitsio.FITS(psffile, 'r')
            header = p[0].read_header()
            p.close()

            if ('RA_TARG' in header):
                ra_p = header['RA_TARG']
            elif 'RA' in header:
                ra_p = header['RA']

            if ('DEC_TARG' in header):
                dec_p = header['DEC_TARG']
            elif 'DEC' in header:
                dec_p = header['DEC']

            p.close()
            distance = 3600.0 * np.sqrt((dec - dec_p)**2.0 + \
                       ((ra - ra_p) * np.cos(dec * Get_R()))**2.0)
        except:
            distance = 9999
            print('Distance between psf and image is 9999')

        return distance

    def HandlePsf(self, UserGivenPsf, ra, dec):

        """Determine the psf used for fitting"""

        if self.repeat:
            if np.abs(ra) == 9999 or np.abs(dec) == 9999:
                distance = 9999
                psffile = c.pfile
            else:
                psffile = c.pfile
                pf.UpdatePsfRaDec(c.pfile)
                distance = self.distance(c.pfile, ra, dec)

        else:
            if UserGivenPsf != 'None':
                psffile = UserGivenPsf
                pf.UpdatePsfRaDec(self.DATADIR, psffile)
                distance = self.distance(psffile, ra, dec)
            elif np.abs(ra) == 9999 or np.abs(dec) == 9999:
                psffile = self.psflist[self.psfcounter]
                pf.UpdatePsfRaDec(self.DATADIR, psffile)
                distance = 9999
                self.psfcounter += 1
            else:
                psffile, distance = pf.getpsf(DATADIR,
                                              psflist, which_psf,
                                              ra, dec)
                distance = distance * 60. * 60.

        return psffile, distance


    def MakeCutOut(self, xcntr, ycntr, 
                   alpha_j, delta_j, 
                   SizeX, SizeY, 
                   TX, TY,
                   gimg, wimg,
                   weightexists):
        """

        Make cutout image. The xcntr, ycntr are like iraf.
        SizeX, SizeY are half size

        """

        ExceedSize = 0
        #All are floor to make the size even number
        xmin = np.floor(xcntr - SizeX)
        ymin = np.floor(ycntr - SizeY)
        xmax = np.floor(xcntr + SizeX)
        ymax = np.floor(ycntr + SizeY)
        if xmin < 0:
            xmin = 0
            cut_xcntr = xcntr
            ExceedSize = 1
        else:
            cut_xcntr = SizeX + np.modf(xcntr)[0]
        if ymin < 0:
            cut_ycntr = ycntr
            ymin = 0
            ExceedSize = 1
        else:
            cut_ycntr = SizeY + np.modf(ycntr)[0]
        if xmax > TX - 1:
            xmax = TX
            ExceedSize = 1
        if ymax > TY - 1:
            ymax = TY
            ExceedSize = 1

        ymin = int(ymin)
        ymax = int(ymax)
        xmin = int(xmin)
        xmax = int(xmax)
        # FIX check c.imagedata or c.ggimage
        data = self.imagedata[ymin:ymax, xmin:xmax]
        # END
        SizeY, SizeX = data.shape

        hdict = {}
        try:
            hdict['RA_TARG'] = alpha_j
            hdict['DEC_TARG'] = delta_j
        except:
            print('Problem updating the Ra and Dec in cutout image')
        if self.EXPTIME != -9999:
            hdict['EXPTIME'] = self.EXPTIME
        else:
            print('EXPTIME have value -9999. Something wrong?')
        if self.RDNOISE != -9999:
            hdict['RDNOISE'] = self.RDNOISE
        else:
            print('RDNOISE have value -9999. Something wrong?')
        if self.GAIN != -9999:
            hdict['GAIN'] = self.GAIN
        else:
            print('GAIN have value -9999. Something wrong?')
        if self.NCOMBINE != -9999:
            hdict['NCOMBINE'] = self.NCOMBINE
        else:
            print('NCOMBINE have value -9999. Something wrong?')

        print(13, gimg)
        print(14, data)
        fits = fitsio.FITS(os.path.join(self.DATADIR, gimg), 'rw')
        fits.write(data, header=hdict)
        fits.close()


        # FIX
        #Making weight image cut
        if weightexists:
            data_wht = self.weightdata[ymin:ymax,xmin:xmax].copy()

            fits = fitsio.FITS(os.path.join(self.DATADIR, wimg), 'rw')
            fits.write(data_wht)
            fits.close()
        else:
            print('Cannot creat weight image. If you supply weight image \
                   please check whether it exists or report a bug')
        #END
        return cut_xcntr, cut_ycntr, SizeX, SizeY, ExceedSize


    def main(self):

        print(self.pixelscale, self.ReSize, self.FixSize)

        #try:
        os.system('rm -f TmpElliMask.fits TmpElliMask1.fits')
        #except:
        #    pass

        

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
                    self.weightdata = wht[0].read()
                    print("whtfile >>> {}".format(self.whtfile))
                elif re.search("weight", self.whtfile.lower()):
                    self.weightdata = 1 / np.sqrt(wht[0].read())
                    print("whtfile >>> {}".format(self.whtfile))
                else:
                    print('Weight file is not understood. Please include ' + \
                          'the word weight/rms to the weight file name. ' + \
                          'If it is weight the rms will be found by 1/sqrt(w)') 
                wht.close()
                weightexists = True
            else:
                weightexists = False
                print('No weight image found\n')
        
        #Initializing psf array. ie. creating psflist from file 
        self.psflist = pf.PSFArr(self.DATADIR, self.psflist)
        if self.decompose:
            for p in self.psflist:
                print('psf', self.DATADIR, p)
                pf.UpdatePsfRaDec(self.DATADIR, p)

        # Writing csv header and finding the output parameters
        params_write = ut.PyMorphOutputParams(self.dbparams, self.decompose)    
        if os.path.exists('result.csv'):
            pass
        else:
            f_res = open("result.csv", "a")
            csvhead = ['{}_{:d}'.format(params_write[par_key][0], par_key)
                       for par_key in params_write.keys()]
            writer = csv.writer(f_res)
            writer.writerow(csvhead)
            f_res.close()

        f_cat = open(self.out_cata, 'w')
        f_failed = open('restart.cat', 'w')
        obj_file = open(self.clus_cata, 'r')  # The file contains the 
                                                     # objects of interest
        pnames = obj_file.readline().split() #The names of the parameters given 
                                             #in the first line in the clus_cata

        # writing a input catalogue (restart.cat) for failed objects
        for pname in pnames:
            f_failed.writelines(['{} '.format(pname)])
        f_failed.writelines(['flag \n'])


        self.psfcounter = 0                  #For getting psf in the case of unknown ra

        pdb = {}                        #The parameter dictionary
        for line_j in obj_file:
            # declare the flag
            self.flag = 0
            
            try:
                values = line_j.split()
            except Exception as e:
                print('No lines in {}'.format(self.clus_cata)) 
                print('Exiting PyMorph')
                sys.exit()

            k = 0
            for pname in pnames:
                pdb[pname] = values[k]
                k += 1

            try:
                gal_id = pdb["gal_id"]
            except:
                print("No gal_id. Try to use gimg")

                try:
                    gal_id = pdb["gimg"].split('.')[0] #gal_id will be 
                                                       #filename without .fits
                except:
                    print("No image or gal_id found in the object catalogue." \
                          "Exiting")
                    os._exit(0)

            self.fstring = '{}_{}'.format(self.rootname, gal_id)
            print('fstring {}'.format(self.fstring))


            try:
                gimg = pdb["gimg"]    #Galaxy cutout
            except:
                print("No gimg given.")
                if os.path.exists(os.path.join(self.DATADIR, 
                                 'I{}.fits'.format(self.fstring))):
                    gimg = 'I{}.fits'.format(self.fstring)
                elif os.path.exists(os.path.join(self.DATADIR, 
                                   '{}.fits'.format(gal_id))): 
                    gimg = '{}.fits'.format(gal_id)
                else:
                    print("No possible gimg found")
                    gimg = 'None'

            try:
                wimg = pdb["wimg"]   #Weight cut
            except:
                wimg = 'None'
                print('Search for weight image (wimg) in galfit config')


            try:
                alpha_j = pdb["ra"]
                h, m, s = alpha.split(':')
                alpha_j = ut.HMSToDeg(int(h), int(m), float(s)) 
            except:
                print("No ra is given")
                alpha_j = -9999

            try:
                delta_j = pdb["dec"]
                d, m, s = delta.split(':')
                delta_j = ut.DMSToDeg(int(d), int(m), float(s))
            except:
                print("No dec is given")
                delta_j = -9999

            if alpha_j == -9999 or delta_j == -9999:
                position = 0
            else:
                position = 1 #Understood position is given

            try:
                z = float(pdb["z"])
            except:
                print("No z is given")
                z = 9999
            try:
                galfit_conf = pdb["cfile"]  #GALFIT configuration file
            except:
                print("No cfile given")
                if self.repeat == True:
                    galfit_conf = 'G_{}.in'.format(self.fstring)
                else:
                    print("Repeat is false. No possible cfile (galfit " + \
                          "config file) is using")
                    galfit_conf = 'None'

            #Reading galfit config file, if it exits, to know the gimg etc.
            if os.path.exists(galfit_conf):
                galfit_params = ut.ReadGalfitConfig(galfit_conf)
                gimg = galfit_params[0].split('/')[-1]
                oimg = galfit_params[1].split('/')[-1]
                wimg = galfit_params[2].split('/')[-1]
                pfile = galfit_params[3].split('/')[-1]
                mimg = galfit_params[4].split('/')[-1]
                cofile = galfit_params[5].split('/')[-1]
            else:
                galfit_conf = 'None'

            #Handling image cutout names
            if self.galcut == True:
                print('Image is >>> {}'.format(gimg))

                gfits = fitsio.FITS(os.path.join(self.DATADIR, gimg), 'r')
                self.imagedata = gfits[0].read()
                self.header0 = gfits[0].read_header()
                gfits.close()

                gheader = ut.CheckHeader(self.header0) #Will set up global header parameters
                self.EXPTIME = gheader[0]
                self.RDNOISE = gheader[1], 
                self.GAIN = gheader[2]
                self.SEx_GAIN = gheader[3]
                self.NCOMBINE = gheader[4]

                if os.path.exists(os.path.join(self.DATADIR, wimg)):
                    wfits = fitsio.FITS(os.path.join(self.DATADIR, wimg))
                    self.weightdata = wfits[0].read()
                    wfits.close()
                    weightexists = True
                else:
                    wimg = None
                    weightexists = False

                print('Using cutouts')

            print(1, gimg)
            TX = self.imagedata.shape[0]
            TY = self.imagedata.shape[1]

            if self.ReSize:
                gimg = 'I{}.fits'.format(self.fstring)
                wimg = 'W{}.fits'.format(self.fstring)

            #The sextractor runs on the cutout before resizing to estimate 
            #shallow sky
            if self.galcut == True:   
                #Given galaxy cutouts if the user provides sextractor catalogue
                #then it will not run SExtractor else do!
                if not os.path.exists(self.sex_cata):
                    try:

                        rs.RunSex(self.sex_params, 
                                  self.PYMORPH_PATH, 
                                  os.path.join(self.DATADIR, gimg),
                                  os.path.join(self.DATADIR, wimg),
                                  os.path.join(self.DATADIR, self.sex_cata),
                                  self.SEx_GAIN, 
                                  sconfig='default')

                        rs.RunSex(self.sex_params,
                                  self.PYMORPH_PATH,
                                  os.path.join(self.DATADIR, gimg),
                                  os.path.join(self.DATADIR, wimg),
                                  os.path.join(self.DATADIR, 
                                              '.Shallow'.format(self.sex_cata)),
                                  self.SEx_GAIN,
                                  sconfig='shallow')

                    except Exception as e:
                        print(type(e))     # the exception instance
                        print(e.args)      # arguments stored in\
                                             # .args
                        print(e)           # __str__ allows args\
                                             # to printed directly
                        print("Something bad happened (Sextractor)!!!!\n\n")
                        print(traceback.print_exc())
                sys.exit() 
            #Lets see whether decorator can be used with ximg and yimg
            try:
                ximg = float(pdb["ximg"])
                if self.ReSize and self.galcut and self.repeat:
                    ximg = TX / 2.0
            except:
                print('No ximg is given or either ReSize or galcut or repeat '\
                      'keywords are False. Trying to find from the cutout if '\
                      'no ra dec in the image header')
                if self.galcut == True & position == 0:
                    ximg = TX / 2.0
                else:
                    ximg = -9999

            try:
                yimg = float(pdb["yimg"])
                if self.ReSize and self.galcut and self.repeat:
                    yimg = TY / 2.0
            except:
                print('No yimg is given or either ReSize or galcut or repeat '\
                      'keywords are False. Trying to find from the cutout if '\
                      'no ra dec in the image header')
                if self.galcut == True & position == 0:
                    yimg = TY / 2.0
                else:
                    yimg = -9999

            try:
                bxcntr = float(pdb["bxcntr"])
            except:
                bxcntr = -9999

            try:
                bycntr = float(pdb["bycntr"])
            except:
                bycntr = -9999

            try:
                mag_zero = float(pdb["mzero"])
            except:
                mag_zero = self.mag_zero

            try:
                UserGivenPsf = pdb["star"]
                print('PSF is assigned individually to galaxies')
            except:
                UserGivenPsf = 'None'

            try:
                UserGivenSky = float(pdb['sky'])
                print('Sky is assigned individually to galaxies')
            except:
                UserGivenSky = np.nan

            # Crashhandling starts
            if self.crashhandler:# & starthandle:
                print("CrashHandler is invoked")

                ut.CrashHandlerToRemove(gal_id, self.fstring, self.OUTDIR)

                try:
                    CrashFlag = float(pdb["flag"])
                    CrashFlag = int(CrashFlag)
                except:
                    CrashFlag = 0

                try:
                    CrashFitFlag = float(pdb["FitFlag"])
                    CrashFitFlag = int(CrashFitFlag)
                except:
                    CrashFitFlag = 0

                if isset(CrashFlag, GetFlag("GALFIT_FAIL")) or \
                   isset(CrashFitFlag, Get_FitFlag("IE_AT_LIMIT")) or \
                   isset(CrashFitFlag, Get_FitFlag("RE_AT_LIMIT")) or \
                   isset(CrashFitFlag, Get_FitFlag("N_AT_LIMIT")) or \
                   isset(CrashFitFlag, Get_FitFlag("EB_AT_LIMIT")) or \
                   isset(CrashFitFlag, Get_FitFlag("ID_AT_LIMIT")) or \
                   isset(CrashFitFlag, Get_FitFlag("RD_AT_LIMIT")) or \
                   isset(CrashFitFlag, Get_FitFlag("ED_AT_LIMIT")):
                    if isset(CrashFlag, GetFlag("FIT_SKY")):
                        fitting[2] = 0
                    else:
                        fitting[2] = 1
                    if isset(CrashFlag, GetFlag("FIT_BULGE_CNTR")) and\
                       isset(CrashFlag, GetFlag("FIT_DISK_CNTR")):
                        pass
                    else:
                        fitting[0] = 1
                        fitting[1] = 1
                if isset(CrashFitFlag, Get_FitFlag("LARGE_CHISQ")):
                    if isset(CrashFlag, GetFlag("FIT_BULGE_CNTR")) and\
                       isset(CrashFlag, GetFlag("FIT_DISK_CNTR")):
                        pass
                    else:
                        fitting[0] = 1
                        fitting[1] = 1
                if isset(CrashFitFlag, Get_FitFlag("FAKE_CNTR")):
                    self.center_deviated = 1


            # Determine Search Radius 
            if self.searchrad is None:
                if position:
                    self.searchrad = '1arc'
                    print('No search radius found. Setting to 1 arc sec')
                else:
                    self.searchrad = '10pix'
                    print('No search radius found. Setting to 10 pix')

            if self.searchrad.endswith('arc'):
                SeaDeg = float(self.searchrad[:-3]) / (60.0 * 60.0)
                SeaPix = 10.0
                print('Arcsec')
            elif self.searchrad.endswith('pix'):
                SeaPix = float(self.searchrad[:-3])
                SeaDeg = pixelscale * SeaPix  / (60.0 * 60.0)

            # first count the number of "potential" targets in the search radius
            SexTargets = 0
            good_objects = []
            bad_objects = []
            good_distance = []
            bad_distance = []
            good_object = ''

            os.system('cp {} sex_{}.txt'.format(self.sex_cata, self.fstring))
            center_distance = 9999 #the distance from the center to the best target

            print(2, gimg)
            for line_s in open(self.sex_cata, 'r'):
                if 1:
                    values = line_s.split()
                    alpha_s = float(values[3])
                    delta_s = float(values[4])
                    if position == 0:
                        alpha_s = 9999
                        delta_s = 9999
                    sex_id = values[0]
                    xcntr  = float(values[1])
                    ycntr  = float(values[2])
                    
                    if (abs(alpha_j - alpha_s) < SeaDeg and \
                       abs(delta_s - delta_j) < SeaDeg) or \
                       (abs(xcntr - ximg) < SeaPix and \
                       abs(ycntr - yimg) < SeaPix):
                        print(xcntr, ximg, ycntr, yimg)   
                        curr_distance = np.sqrt(np.min([(xcntr - ximg)**2 + \
                                (ycntr - yimg)**2, (alpha_j - alpha_s)**2 + \
                                (delta_s - delta_j)**2]))
                        print("Candidate distance: {:3f}".format(curr_distance)) 

                        SexTargets +=1

                        if curr_distance < center_distance:
                            center_distance = curr_distance
                            print("New Preferred target!!")
                            good_object = line_s

                        print("Target distance: {:.3f}".format(center_distance))
                #except:
                #    if values[0].strip().isdigit():
                #        print('Something happend in the pipeline. ' + \
                        #              'Check error.log')
                    #else:
                    #    pass
            if len(good_object) == 0:
                # No suitable target found
                print("No Target Found!!!!")
                good_object = ' 9999  9999 9999  9999 9999  9999 9999  9999 9999  9999 9999  0 9999  9999 9999  9999 9999 9999 9999\n'
                self.flag = SetFlag(self.flag, GetFlag('NO_TARGET'))  

            # now fit best object            
            #XXX
            if 1:#try:
                target = mf.GetSExObj(line_s=good_object)

                sex_id = target.sex_num
                print("SExtractor ID >>> {}".format(sex_id))

                if position == 0:
                    alpha_s = 9999
                    delta_s = 9999
                
                SexHalfRad = target.radius #Sex halfrad 
                SexPosAng = target.pos_ang
                pos_ang = target.pos_ang
                axis_rat = target.bbya #axis ration b/a
                eg = target.eg
                if eg <= 0.05:
                    eg = 0.07
                major_axis = target.maj_axis

                SexMagAuto = target.mag
                SexMagAutoErr = target.mag_err
                UMag = -1e5
                LMag = -1e5
                URe = 1e0
                URd = 1e0
                if UMag > -9999.0:
                    UMag = SexMagAuto - 7.0
                if LMag < 9999.0:
                    LMag = SexMagAuto + 7.0

                # FIX
                if URe < 9999.0:
                    URe = SexHalfRad * 50.0
                if URd < 9999.0:
                    URd = SexHalfRad * 50.0
                # END
                print(1)
                # The shallow sky from the shallow run. If no shallow
                # sky it uses the deep sky.
                if os.path.exists('{}.Shallow'.format(self.sex_cata)):
                    f_sex_shallow = open('{}.Shallow'.format(self.sex_cata), 'r')
                    for line_shallow in f_sex_shallow:
                        v_shallow = line_shallow.split()
                        try:
                            if v_shallow[19] == values[0]:
                                SexSky = float(v_shallow[10])
                        except:
                            pass
                    f_sex_shallow.close()
                else:
                    SexSky = float(values[10])

                if not np.isnan(UserGivenSky):
                    skyval = UserGivenSky #using the sky supplied by the user
                else:
                    skyval = SexSky #use the sky value from sextractor
                    
                if (position & alpha_j == -9999) | (position & delta_j == -9999):
                    alpha_j = alpha_s
                    delta_j = delta_s

                ut.WriteError('\n\n###########   {} ###########\n'.format(gal_id))
                print(0, gimg)

                run = 1 #run =1 if pipeline runs sucessfuly
                # Adding initial setup to the flag
                if self.repeat:
                    self.flag = SetFlag(self.flag, GetFlag('REPEAT'))
                if self.fitting[0] and 'bulge' in self.components:
                    self.flag = SetFlag(self.flag, GetFlag('FIT_BULGE_CNTR'))
                if self.fitting[1] and 'disk' in self.components:
                    self.flag = SetFlag(self.flag, GetFlag('FIT_DISK_CNTR'))
                if self.fitting[2]:
                    self.flag = SetFlag(self.flag, GetFlag('FIT_SKY'))
                
                # Calculating the cutout size (half size).
                # SizeX, SizeY return from MakeCutOut are full sizes
                SizeX, SizeY = self.FindCutSize(TX / 2, TY / 2,
                                                SexHalfRad, 
                                                pos_ang, axis_rat) 
                print('Calculated half sizes {}, {}'.format(SizeX, SizeY))

                # Finding psf and the distance between psf and image
                if self.decompose:
                    psffile, distance = self.HandlePsf(UserGivenPsf, 
                                                       alpha_j, delta_j)
                else:
                    psffile, distance = 'None', 9999

                print('psffile: {}, distance: {}'.format(psffile, distance))

                # For the new run
                if self.repeat:
                    cut_xcntr = TX / 2.0
                    cut_ycntr = TY / 2.0
                    if os.path.exists(mimg):
                        pass
                    else:
                        MF = MaskFunc(mimg, 
                                      cut_xcntr, cut_ycntr, 
                                      SizeX, SizeY, 
                                      good_object)
                        MF.gmask(self.threshold, self.thresh_area,
                                 0, NoMask=False)

                else:
                    if self.galcut:
                        if self.ReSize: 
                            #Can take a function
                            #Can check whether a decorator can be used
                            if os.path.exists(gimg):
                                ut.WriteError('The file {} exists\n'.format(gimg))
                                run = 0
                                break #Breaking the sextractor loop
                            try:
                                cut_params = self.MakeCutOut(xcntr, ycntr,
                                                             alpha_j, delta_j,
                                                             SizeX, SizeY,
                                                             TX, TY,
                                                             gimg, wimg,
                                                             weightexists)
                                cut_xcntr = cut_params[0]
                                cut_ycntr = cut_params[1]
                                SizeX = cut_params[2]
                                SizeY = cut_params[3]
                                ExceedSize = cut_params[4]
                            except  Exception as e:
                                print(type(e))     # the exception instance
                                print(e.args)      # arguments stored in\
                                                     # .args
                                print(e)           # __str__ allows args\
                                                     # to printed directly
                                ut.WriteError('Cutout exists!')
                                print(traceback.print_exc())
                                break
                        else:
                            cut_xcntr = xcntr
                            cut_ycntr = ycntr 
                            SizeX = TX
                            SizeY = TY
                            ExceedSize = 0
                    else:
                        print(11, gimg)
                        gimg = 'I{}.fits'.format(self.fstring)
                        wimg = 'W{}.fits'.format(self.fstring)
                        #Can take a function
                        if os.path.exists(gimg):
                            ut.WriteError('The file {} exists\n'.format(gimg))
                            run = 0
                            break #Breaking the sextractor loop
                            print(3)

                        # Sizes are total size
                        print(12, gimg)
                        #XXX
                        if 1:#try: 
                            cut_params = self.MakeCutOut(xcntr, ycntr, 
                                                       alpha_j, delta_j, 
                                                       SizeX, SizeY, 
                                                       TX, TY, 
                                                       gimg, wimg,
                                                       weightexists)
                            cut_xcntr = cut_params[0]
                            cut_ycntr = cut_params[1]
                            SizeX = cut_params[2]
                            SizeY = cut_params[3]
                            ExceedSize = cut_params[4]
                            
                        #except  Exception as e:
                        else:
                            print(type(e))     # the exception instance
                            print(e.args)      # arguments stored in\
                                                 # .args
                            print(e)           # __str__ allows args\
                                                 # to printed directly
                            ut.WriteError('Cutout exists!')
                            print(traceback.print_exc())
                            break

                    print(2, gimg)

                    print('Center of cutimage and exceed size ', \
                          cut_xcntr, cut_ycntr, ExceedSize)
                    print('Full Sizes ', SizeX, SizeY)

                    if self.galcut & self.ReSize == 0:
                        pass
                    elif ExceedSize:
                        self.flag = SetFlag(self.flag, GetFlag('EXCEED_SIZE'))

                    # Runs sextractor to find the segmentation map
                    ##RunSegSex(os.path.join(self.DATADIR, gimg))

                    # This creates ellipse mask and used for 
                    # ellipse fitting and casgm
                    seg_cata = os.path.join(self.DATADIR,
                                            '{}.seg'.format(self.sex_cata))

                    
                    gimg = os.path.join(self.DATADIR, gimg)
                    wimg = os.path.join(self.DATADIR, wimg)

                    seg_file = os.path.join(self.OUTDIR, 'seg.fits')

                    #if os.path.exists(seg_file):
                    #    os.remove(seg_file)
                    if os.path.exists(seg_cata):
                        os.remove(seg_cata)

                    RunSex(self.sex_params, self.SEX_PATH, gimg, wimg, 
                           seg_cata, self.SEx_GAIN, sconfig='seg')
                    #sys.exit()
                    EM = ElliMaskFunc(cut_xcntr, cut_ycntr, 1,
                                      center_limit=5., seg_limit=1e-5)
                    EM.emask(seg_file, seg_cata, self.fstring)

                    sys.exit()
                    # Fitting ellipse task or the manual 1d finder
                    mimg = 'M_{}.fits'.format(self.fstring)
                    print(3, gimg)
                    if self.decompose:
                        print(4, gimg)
                        ut.FindEllipse(gimg, xcntr, ycntr, SexHalfRad,
                                      SexPosAng, axis_rat, skyval, 
                                      self.fstring, output=False)

                        MaskFunc(gimg, 
                                 cut_xcntr, cut_ycntr, 
                                 SizeX, SizeY, 
                                 good_object)

                        galfit_conf = 'G_{}.in'.format(self.fstring)
                        oimg = 'O_{}.fits'.format(self.fstring)

                        print(4, gimg)

                        CF = GalfitConfigFunc(self.DATADIR,
                                              gimg, wimg,  
                                              cut_xcntr, cut_ycntr, 
                                              SizeX, SizeY,
                                              self.components, self.fitting,
                                              psffile,
                                              good_object, 
                                              self.sex_cata,
                                              skyval, self.mag_zero, self.flag)
                        CF.write_config(self.fstring, 
                                        self.threshold, self.thresh_area, 
                                        self.center_deviated, self.galfitv, 
                                        center_constrain=2.0,
                                        avoidme=self.avoidme,
                                        LMag = self.LMag,
                                        UMag = self.UMag,
                                        LN = self.LN,
                                        UN = self.UN,
                                        LRe = self.LRe,
                                        URe = self.URe,
                                        LRd = self.LRd,
                                        URd = self.URd,
                                        bdbox = self.bdbox,
                                        bbox = self.bbox,
                                        dbox = self.dbox,
                                        devauc = self.devauc)
                        #continue
                

                # Estimates sky parameters
                #XXX
                if 1:
                    sky_values = FindYetSky(self.sex_params, self.SEX_PATH, 
                                            os.path.join(self.DATADIR, gimg),
                                            os.path.join(self.DATADIR, wimg),
                                            seg_file,
                                            cut_xcntr, cut_ycntr,
                                            os.path.join(self.DATADIR,
                                                '{}.seg'.format(self.sex_cata)),
                                            self.SEx_GAIN,
                                            center_err=5., median_std=1.3,
                                            sconfig='seg')
                    #sky_values = [9999, 9999, 9999, 9999, 9999, 9999]
                    SexySky = sky_values[0]
                    SkyYet = sky_values[1]
                    SkyMed = sky_values[2]
                    SkyMin = sky_values[3]
                    SkyQua = sky_values[4]
                    SkySig = sky_values[5]

                    if SkyMin != 9999:
                        SkyMin = SkyYet * 1.0
                        SkySig = SkySig * 1.0
                    else:
                        SkyMin = SexSky * 1.0 
                        SkySig = np.sqrt(np.abs(SexSky))
                    print('Sky Sigma {} '.format(SkySig))
                    print('SkyMin {} SexSky {}'.format(SkyMin, SexSky))
                else:
                #except Exception as e:
                    print(type(e))     # the exception instance
                    print(e.args)      # arguments stored in\
                                         # .args
                    print(e)           # __str__ allows args\
                                         # to printed directly
                    ut.WriteError('YetSky estimation failed\n')
                    print(traceback.print_exc())

                # Estimate CASGM  
                if self.cas:
                    cas_values = ut.HandleCasgm(gimg, 
                                                cut_xcntr, cut_ycntr, 
                                                alpha_j, delta_j, 
                                                z, 
                                                SizeX, SizeY, 
                                                good_object, 
                                                bxcntr, bycntr)

                    C = cas_values[0]
                    C_err = cas_values[1]
                    A = cas_values[2]
                    A_err = cas_values[3]
                    S = cas_values[4]
                    S_err = cas_values[5]
                    G = cas_values[6]
                    M = cas_values[7]
                else:
                    C, C_err, A, A_err = 9999, 9999, 9999, 9999
                    S, S_err, G, M = 9999, 9999, 9999, 9999

                # Removing CASGM temp files
                for f in ['BMask.fits', 'MRotated.fits', 
                          'MaskedGalaxy.fits', 'Rotated.fits']:
                    if os.access(f, os.F_OK):
                        os.remove(f)

                # Decomposition
                if os.access('fit.log', os.F_OK):
                    os.remove('fit.log')

                if self.decompose:
                    try:
                        if self.galfit & self.detail:
                            ConfigIter(gimg, wimg, 
                                       cut_xcntr, cut_ycntr, 
                                       SizeX, SizeY, 
                                       good_object, 
                                       psffile, z)
                        elif self.galfit:
                            cmd = '{} {}'.format(self.GALFIT_PATH, galfit_conf)
                            os.system(cmd)
                            f_fit = open('fit2.log','a')
                            if os.path.exists('fit.log'):
                                for line in open('fit.log','r'):
                                    f_fit.writelines([str(line)])
                            f_fit.close()
                            # FIX
                            # Mainly the FindEllipse problem
                            ut.OImgFindEllipse(gimg, oimg, 
                                            cut_xcntr, cut_ycntr, 
                                            self.fstring, self.repeat,
                                            self.flag)
                            #                good_object, self.fstring,
                            #                self.components, 
                            #                SexSky, self.repeat, 
                            #                self.flag, run, 1)
                            #END
                    except Exception as e:
                        print(type(e))     # the exception instance
                        print(e.args)      # arguments stored in\
                                             # .args
                        print(e)           # __str__ allows args\
                                             # to printed directly
                        print("Something bad happened (GALFIT)\n\n")
                        print(traceback.print_exc())

                
                            
                    if run == 1:
                        ut.WriteError('((((( Decomposition Successful )))))\n')

                sys.exit()
                try:
                    Goodness = -9999
                    if os.access('P_' + self.fstring + '.png', os.F_OK):	
                        os.remove('P_' + self.fstring + '.png')
                    GoodNess = PlotFunc(oimg, mimg, 
                                        cut_xcntr, cut_ycntr, 
                                        SexSky, SkySig)
                    Goodness = GoodNess.plot_profile
                except:
                    ut.WriteError('Error in plotting \n')
                    if mimg == 'None':
                        ut.WriteError('Could not find Mask image for plottong \n')
                    run = 0	
                    self.flag = SetFlag(self.flag, GetFlag('PLOT_FAIL'))
                
                if isset(self.flag, GetFlag("GALFIT_FAIL")): #or \
                   #isset(flag, GetFlag("LARGE_CHISQ")) or \
                   #isset(flag, GetFlag("FAKE_CNTR")) or \
                   #isset(flag, GetFlag("BULGE_AT_LIMIT")) or \
                   #isset(flag, GetFlag("DISK_AT_LIMIT")):
                    failedvalues = line_j.split()
                    for fv in failedvalues:
                        f_failed.writelines([' '.format(fv)])
                    f_failed.writelines(['{}\n'.format(self.flag)])

                f_cat.writelines(['{} '.format(gal_id)])
                f_cat.write(good_object)


                sys.exit()
                try:
                    WriteParams(ParamToWrite, 
                                gimg, 
                                cut_xcntr, cut_ycntr, 
                                distance, 
                                alpha_j, delta_j, 
                                z, Goodness, 
                                C, C_err, A, A_err, S, S_err, G, M, 
                                EXPTIME)
                except Exception as e:
                    print(type(e))     # the exception instance
                    print(e.args)      # arguments stored in\
                    # .args
                    print(e)        # __str__ allows args\
                    # to printed directly
                    print("Something bad happened (Writing)!!!!\n\n")
                    print(traceback.print_exc())

                    #except:
                    #    ut.WriteError('Error in writing html\n')
                    #    run = 0

                #The following removes all the temporary files 
                # after every fit
                ToClean = 0

                if ToClean:
                    ClaCli = ['E*fits', 'E*txt', 'G_*', 'I*fits', '*.pl',
                              '*.con', 'M_*', 'O*fits', 'O*txt', 'P*png', 
                              'R*html', 'error.log', 'galfit.*', 'Tmp*', 'SO*', 
                              'agm_r*', 'BMask.fits', 'MaskedGalaxy.fits', 
                              'MRotated.fits', 'B.fits', 'GalEllFit.fits', 
                              'AResidual.fits', 'ellip', 'err', 'BackMask.fits']
                    for f in ClaCli:
                        if os.access(f, os.F_OK):
                            os.remove(f)

                for myfile in ['ellip','err','test.tab']:
                    if os.access(myfile,os.F_OK):
                        os.remove(myfile)

            #except Exception as e:
            else:
                print(type(e))     # the exception instance
                print(e.args)      # arguments stored in\
                                            # .args
                print(e)           # __str__ allows args\
                                             # to printed directly
                print("Something bad happened (general sextractor)!!!!\n\n")
                print(traceback.print_exc())

                #NOTE CHANGE THIS ADD ADDITIONAL IF len > 0 HERE

            if values[0].strip().isdigit():
                print('Something happend in the pipeline. ' + \
                      'Check error.log 2')
            else:
                pass

            if self.galcut:
                if os.access(self.sex_cata, os.F_OK):
                    os.remove(self.sex_cata)

        #except Exception as e:
        #    print(type(e))     # the exception instance
        #    print(e.args)      # arguments stored in\
            # .args
        #    print(e)           # __str__ allows args\
            # to printed directly
        #    print("something bad happened (general object search)!!!!\n\n")
        #    print(traceback.print_exc())
        f_cat.close()
        f_failed.close()






    def FindAndFit(self):
        fc = open(self.clus_cata, 'w')
        if self.redshift == 9999:
            fc.writelines(['gal_id ra1 dec1 mag\n'])
        else:
            fc.writelines(['gal_id ra1 dec1 mag z\n'])

        for line in open(self.sex_cata, 'r'):

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
                self.EXPTIME = gheader[0]
                self.RDNOISE = gheader[1],
                self.GAIN = gheader[2]
                self.SEx_GAIN = gheader[3]
                self.NCOMBINE = gheader[4]
            else:
                self.SEx_GAIN = 1.

        except IOError as e:
            print(self.imagefile, "I/O error ({}): {}".format(e.args[0], e.args[1]))
            os._exit(0)
        
        print(self.sex_cata)
        if os.path.exists(self.sex_cata):
            pass
        elif self.galcut == False:
            print('The SExtractor catalogue for your frame is NOT found. ' \
                  'One is being made using the default values. It is always '\
                  'recommended to make SExtractor catalogue by YOURSELF as '\
                  'the pipeline keeps the sky value at the SExtractor value '\
                  'during the decomposition.')
            RunSex(self.sex_params, self.SEX_PATH, 
                   self.imagefile, self.whtfile,
                   self.sex_cata, self.SEx_GAIN,
                   sconfig='default')

    #The old function for psfselect = 2
    #    if psfselect == 2:
    #        center_deviated = 0
    #        starthandle = 0
    #        os.system('ds9 &')
    #        time.sleep(2)
    #        run_psfselect(imagefile, DATADIR, clus_cata, self.galcut)
    #        os.system('xpaset -p ds9 quit')
    #        psflist = '@psflist.list'
    #        FindAndFit()
    #        main()
    #        if crashhandler:
    #            starthandle = 1
    #            os.system('mv restart.cat CRASH.CAT')
    #            clus_cata = 'CRASH.CAT' 
    #            main()
    #New function for psfselect=2 non-interactive for webservice
        print(self.psfselect)
        print(self.clus_cata)
        if self.psfselect == 2:
            self.Interactive = 0
            self.center_deviated = 0
            self.starthandle = 0
            run_psfselect(self.imagefile, 
                          self.DATADIR, 
                          self.clus_cata, 
                          self.galcut)
            psflist = os.path.join(self.DATADIR, '@psflist.list')
            self.FindAndFit()
            self.main()
            if self.crashhandler:
                starthandle = 1
                os.system('mv restart.cat CRASH.CAT')
                self.clus_cata = 'CRASH.CAT' 
                self.main()

    #The old function for psfselect=1
        elif self.psfselect == 1:
            self.Interactive = 1
            os.system('ds9 &')
            time.sleep(2)
            run_psfselect(self.imagefile, 
                          self.DATADIR, 
                          self.clus_cata, 
                          self.galcut)
            os.system('xpaset -p ds9 quit')
    #new function for psfselect=1 non-interactive for webservice (abhishek rawat)
    #    elif psfselect == 1:
    #        Interactive = 0
    #        run_psfselect(imagefile, DATADIR, clus_cata, self.galcut)
            



        elif self.psfselect == 0:
            self.center_deviated = 0
            self.starthandle = 0
            self.FindAndFit()
            print(self.sex_cata)
            self.main()
            if self.crashhandler:
                self.starthandle = 1
                os.system('mv restart.cat CRASH.CAT')
                self.clus_cata = 'CRASH.CAT' 
                self.main()

        os.chdir(thisdir)







p = PyMorph()
p.pymorph()
#p.main()
