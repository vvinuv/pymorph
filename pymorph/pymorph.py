
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
import subprocess
from multiprocessing import Pool
#import pymorphutils as ut
from .pymorphutils import get_header, output_params
from .flagfunc import GetFlag, isset, SetFlag

from .ellimaskfunc_easy import ElliMaskFunc

from .maskfunc_easy import MaskFunc

from .configfunc import GalfitConfigFunc
#from configtwostep import ConfigIter
from .yetbackfunc import FindYetSky
#from plotfunc import PlotFunc
from .runsexfunc import PySex
#from writehtmlfunc import WriteParams
from .psffunc import psfarr, update_psf_ra_dec  
from .weightimage import return_sigma
#import mask_or_fit as mf
from .pipeline import Pipeline
       
        


class InitializeParams(object):

    def _initilize_params(self, config_file='config.ini'):

        c = configparser.ConfigParser()
        c.read(config_file)
        
        try:
            self.DATADIR = c.get('imagecata', 'datadir')
        except:
            self.DATADIR = os.getcwd()

        try:
            self.OUTDIR = c.get('imagecata', 'outdir')
        except:
            self.OUTDIR = os.getcwd()

        self.imagefile = c.get('imagecata', 'imagefile')
        self.imagefile = os.path.join(self.DATADIR, self.imagefile)

        self.whtfile = c.get('imagecata', 'whtfile')
        self.whtfile = os.path.join(self.DATADIR, self.whtfile)

        self.sex_cata = c.get('imagecata', 'sex_cata')
        self.sex_cata = os.path.join(self.DATADIR, self.sex_cata)

        self.obj_cata = c.get('imagecata', 'obj_cata')
        self.obj_cata = os.path.join(self.DATADIR, self.obj_cata)

        self.out_cata = c.get('imagecata', 'out_cata')
        self.out_cata = os.path.join(self.OUTDIR, self.out_cata)


        self.rootname = c.get('imagecata', 'rootname')

        self.GALFIT_PATH = c.get('external', 'GALFIT_PATH')
        self.SEX_PATH = c.get('external', 'SEX_PATH')

        self.findandfit = c.getboolean('modes', 'findandfit')
        self.repeat = c.getboolean('modes', 'repeat')
        self.galcut = c.getboolean('modes', 'galcut')
        self.decompose = c.getboolean('modes', 'decompose')
        self.detail = c.getboolean('modes', 'detail')
        self.galfit = c.getboolean('modes', 'galfit')
        self.cas = c.getboolean('modes', 'cas')
        self.crashhandler = c.getboolean('modes', 'crashhandler')


        try:
            components = c.get('galfit', 'components').split(',')
            self.components = [cm.strip() for cm in components]
        except:
            print("components undefined. Asuming bulge+disk model")
            self.components = ['bulge', 'disk']

        self.devauc = c.getboolean('galfit', 'devauc')
        self.fitting = c.get('galfit', 'fitting')
        self.fitting = [int(tf) for tf in self.fitting.split(',')]

        self.chi2nu_limit = c.getfloat('diagnosis', 'chi2nu_limit')
        self.goodness_limit = c.getfloat('diagnosis', 'goodness_limit')
        self.center_deviation_limit = c.getfloat('diagnosis', 'center_deviation_limit')
        self.remove_obj_boundary = c.getfloat('diagnosis', 'remove_obj_boundary')
        self.do_plot = c.get('diagnosis', 'do_plot')

        try:
            size = c.get('size', 'size_list')
            size = [int(s) for s in size.split(',')]
            self.ReSize = size[0]
            self.ImSize = size[1]
            self.FracRad = size[2]
            self.Square = size[3]
            self.FixSize = size[4]
        except:
            self.ReSize = c.getint('size', 'size_list')
            self.ImSize = 0
            self.Square = 1
            self.FracRad = 10
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
        #Read PSF 
        psflist = c.get('psf', 'psflist')
        if psflist.startswith('@'):
            with open(psflist[1: ], 'r') as fpsf:
                dpsf = fpsf.read()
            self.psflist = [tdpsf.strip() for tdpsf in dpsf.split('\n')[:-1]]
        else:
            self.psflist = [fpsf.strip() for fpsf in psflist.split(',')]
        
        self.which_psf = c.getint('psf', 'which_psf')
        self.area_obj = c.getint('psf', 'which_psf')

        self.manual_mask = c.getint('mask', 'manual_mask')
        self.mask_reg = c.getfloat('mask', 'mask_reg')
        self.thresh_area = c.getfloat('mask', 'thresh_area')
        self.threshold = c.getfloat('mask', 'threshold')
        self.avoidme = c.getfloat('mask', 'avoidme')
        self.mag_zero = c.getfloat('mag', 'mag_zero')
        self.maglim = c.get('mag', 'maglim').split(',')
        self.maglim = [float(mlim) for mlim in self.maglim]

        self.stargal_prob_lim = c.getfloat('findfit', 'stargal_prob_lim')

        self.host = c.get('db', 'host')
        self.database = c.get('db', 'database')
        self.table = c.get('db', 'table')
        self.user = c.get('db', 'user')
        self.password = c.get('db', 'password')
        self.dbparams = c.get('db', 'dbparams')

        self.galfitv = c.get('version', 'galfitv')
        self.FirstCreateDB = 1 #Won't create table when FirstCreateDB=0
        self.FILTER = 'UNKNOWN'


        #if self.set_lim_mag_rad_sex == False:
        self.NMag = c.getfloat('params_limit', 'NMag')
        self.NRadius = c.getfloat('params_limit', 'NRadius')

        self.LN = c.getfloat('params_limit', 'LN')
        self.UN = c.getfloat('params_limit', 'UN')
        
        self.center_constrain = c.getfloat('params_limit', 'center_constrain')

        self.bdbox = c.get('params_limit', 'bdbox')
        self.bbox = c.get('params_limit', 'bbox')
        self.dbox = c.get('params_limit', 'dbox')
        self.no_mask = c.get('mask', 'no_mask')
        self.norm_mask = c.get('mask', 'norm_mask')

        self.SEx_PIXEL_SCALE = self.pixelscale
        self.SEx_SEEING_FWHM = self.pixelscale * 3.37
        self.SEx_MAG_ZEROPOINT = self.mag_zero

        if c.has_section('sextractor'):
            self.SEx_DETECT_MINAREA = c.get('sextractor', 'DETECT_MINAREA') 
            self.SEx_DETECT_THRESH =  c.get('sextractor', 'DETECT_THRESH')
            self.SEx_ANALYSIS_THRESH = c.get('sextractor', 'ANALYSIS_THRESH')
            self.SEx_FILTER = c.get('sextractor', 'FILTER')
            self.SEx_FILTER_NAME = c.get('sextractor', 'FILTER_NAME') 
            self.SEx_DEBLEND_NTHRESH = c.get('sextractor', 'DEBLEND_NTHRESH') 
            self.SEx_DEBLEND_MINCONT = c.get('sextractor', 'DEBLEND_MINCONT')
            self.SEx_PHOT_FLUXFRAC = c.get('sextractor', 'PHOT_FLUXFRAC') 
            self.SEx_BACK_SIZE = c.get('sextractor', 'BACK_SIZE')
            self.SEx_BACK_FILTERSIZE = c.get('sextractor', 'BACK_FILTERSIZE') 
            self.SEx_BACKPHOTO_TYPE = c.get('sextractor', 'BACKPHOTO_TYPE')
            self.SEx_BACKPHOTO_THICK = c.get('sextractor', 'BACKPHOTO_THICK')
            self.SEx_WEIGHT_TYPE = c.get('sextractor', 'WEIGHT_TYPE')
        else:
            #Sextractor configuration
            self.SEx_DETECT_MINAREA = 6
            self.SEx_DETECT_THRESH = 1.5
            self.SEx_ANALYSIS_THRESH =1.5
            self.SEx_FILTER = 'Y'
            self.SEx_FILTER_NAME = 'default.conv'
            self.SEx_DEBLEND_NTHRESH = 32
            self.SEx_DEBLEND_MINCONT = 0.005
            self.SEx_PHOT_FLUXFRAC = 0.5
            self.SEx_BACK_SIZE = 64
            self.SEx_BACK_FILTERSIZE = 3
            self.SEx_BACKPHOTO_TYPE = 'GLOBAL'
            self.SEx_BACKPHOTO_THICK = 24
            self.SEx_WEIGHT_TYPE = 'MAP_RMS'


        self.pymorph_config = self.get_pymorph_config_dict()
        
        self.sex_config = self.get_sex_config_dict()


    def get_pymorph_config_dict(self):

        pymorph_config = dict()

        pymorph_config['DATADIR'] = self.DATADIR
        pymorph_config['OUTDIR'] = self.OUTDIR
        pymorph_config['imagefile'] = self.imagefile
        pymorph_config['whtfile'] = self.whtfile
        pymorph_config['sex_cata'] = self.sex_cata
        pymorph_config['obj_cata'] = self.obj_cata
        pymorph_config['out_cata'] = self.out_cata
        pymorph_config['rootname'] = self.rootname
        pymorph_config['GALFIT_PATH'] = self.GALFIT_PATH
        pymorph_config['SEX_PATH'] = self.SEX_PATH
        pymorph_config['findandfit'] = self.findandfit
        pymorph_config['repeat'] = self.repeat
        pymorph_config['galcut'] = self.galcut
        pymorph_config['decompose'] = self.decompose
        pymorph_config['detail'] = self.detail
        pymorph_config['galfit'] = self.galfit
        pymorph_config['cas'] = self.cas
        pymorph_config['crashhandler'] = self.crashhandler
        pymorph_config['components'] = self.components
        pymorph_config['devauc'] = self.devauc
        pymorph_config['fitting'] = self.fitting
        pymorph_config['chi2nu_limit'] = self.chi2nu_limit
        pymorph_config['goodness_limit'] = self.goodness_limit
        pymorph_config['center_deviation_limit'] = self.center_deviation_limit
        pymorph_config['center_constrain'] = self.center_constrain
        pymorph_config['remove_obj_boundary'] = self.remove_obj_boundary
        pymorph_config['do_plot'] = self.do_plot
        pymorph_config['ReSize'] = self.ReSize
        pymorph_config['ImSize'] = self.ImSize
        pymorph_config['FracRad'] = self.FracRad
        pymorph_config['Square'] = self.Square
        pymorph_config['FixSize'] = self.FixSize
        pymorph_config['searchrad'] = self.searchrad
        pymorph_config['pixelscale'] = self.pixelscale
        pymorph_config['H0'] = self.H0
        pymorph_config['WM'] = self.WM
        pymorph_config['WV'] = self.WV
        pymorph_config['redshift'] = self.redshift
        pymorph_config['psfselect'] = self.psfselect
        pymorph_config['stargal_prob'] = self.stargal_prob
        pymorph_config['star_size'] = self.star_size
        pymorph_config['psflist'] = self.psflist
        pymorph_config['which_psf'] = self.which_psf
        pymorph_config['area_obj'] = self.area_obj
        pymorph_config['manual_mask'] = self.manual_mask
        pymorph_config['mask_reg'] = self.mask_reg
        pymorph_config['thresh_area'] = self.thresh_area
        pymorph_config['threshold'] = self.threshold
        pymorph_config['mag_zero'] = self.mag_zero
        pymorph_config['maglim'] = self.maglim
        pymorph_config['host'] = self.host
        pymorph_config['database'] = self.database
        pymorph_config['table'] = self.table
        pymorph_config['user'] = self.user
        pymorph_config['password'] = self.password
        pymorph_config['dbparams'] = self.dbparams
        pymorph_config['galfitv'] = self.galfitv
        pymorph_config['FirstCreateDB'] = self.FirstCreateDB
        pymorph_config['FILTER']  = self.FILTER
        pymorph_config['NMag'] = self.NMag
        pymorph_config['NRadius']  = self.NRadius
        pymorph_config['LN'] = self.LN
        pymorph_config['UN'] = self.UN
        pymorph_config['bdbox'] = self.bdbox
        pymorph_config['bbox'] = self.bbox
        pymorph_config['dbox'] = self.dbox
        pymorph_config['avoidme']  = self.avoidme
        pymorph_config['no_mask'] = self.no_mask
        pymorph_config['norm_mask'] = self.norm_mask

        return pymorph_config


    def get_sex_config_dict(self):

        sex_config = dict()

        sex_config['SEx_PIXEL_SCALE'] = self.pixelscale
        sex_config['SEx_SEEING_FWHM'] = self.pixelscale * 3.37
        sex_config['SEx_MAG_ZEROPOINT'] = self.mag_zero
        sex_config['SEx_DETECT_MINAREA'] = self.SEx_DETECT_MINAREA
        sex_config['SEx_DETECT_THRESH'] = self.SEx_DETECT_THRESH
        sex_config['SEx_ANALYSIS_THRESH'] = self.SEx_ANALYSIS_THRESH
        sex_config['SEx_FILTER'] = self.SEx_FILTER
        sex_config['SEx_FILTER_NAME'] = self.SEx_FILTER_NAME
        sex_config['SEx_DEBLEND_NTHRESH'] = self.SEx_DEBLEND_NTHRESH
        sex_config['SEx_DEBLEND_MINCONT'] = self.SEx_DEBLEND_MINCONT
        sex_config['SEx_PHOT_FLUXFRAC'] = self.SEx_PHOT_FLUXFRAC
        sex_config['SEx_BACK_SIZE'] = self.SEx_BACK_SIZE
        sex_config['SEx_BACK_FILTERSIZE'] = self.SEx_BACK_FILTERSIZE
        sex_config['SEx_BACKPHOTO_TYPE'] = self.SEx_BACKPHOTO_TYPE
        sex_config['SEx_BACKPHOTO_THICK'] = self.SEx_BACKPHOTO_THICK
        sex_config['SEx_WEIGHT_TYPE'] = self.SEx_WEIGHT_TYPE

        return sex_config



class PyMorph(InitializeParams):
    
    #__version__ = 3.0

    def __init__(self, config_file='config.ini'):

        super()._initilize_params(config_file=config_file) 
        galfitv = self.pymorph_config['galfitv'].split('.')
        galfitv = galfitv[0] + '.' + galfitv[1]
        self.galfitv = float(galfitv)
        self.pymorph_config['galfitv'] = float(galfitv)
        self.final_result_file = os.path.join(self.OUTDIR, 'result.csv')
        self.restart_file = os.path.join(self.OUTDIR, 'restart.cat') 

    def _indexfile(self):
        #Initialize index.html
        if os.path.exists('index.html'):
            pass
        else:
            indexfile = open('index.html', 'w')
            indexfile.writelines(['<HTML>\n<BODY>\n'])
            indexfile.writelines(['</BODY></HTML>'])
            indexfile.close()



    def main_thread(self):

        #print(self.pixelscale, self.ReSize, self.FixSize)

        #print('PSF 2', self.psflist)
        self._indexfile()
        

        print(self.imagefile)
        #print(self.imagedata)
        print('Image done')
        print(self.whtfile)
        print('Weight done')
        #print('PSF 1', self.psflist)
        #Initializing psf array. ie. creating psflist from file 
        self.psflist = psfarr(self.DATADIR, self.psflist)
        #if self.decompose:
        #    for p in self.psflist:
        #        print('psf', self.DATADIR, p)
        #        update_psf_ra_dec(self.DATADIR, p)

        #XXX
        #self.P.psflist = self.psflist

        # Writing csv header and finding the output parameters
        params_to_write = output_params(self.dbparams, self.decompose)    

        #print('P4')
        #XXX
        #self.P.params_to_write = params_to_write

        #print(self.final_result_file) 
        if os.path.exists(self.final_result_file):
            pass
        else:
            f_res = open(self.final_result_file, "w")
            csvhead = ['{}_{:d}'.format(params_to_write[par_key][0], par_key)
                       for par_key in params_to_write.keys()]
            writer = csv.writer(f_res)
            writer.writerow(csvhead)
            f_res.close()

        #print(self.out_cata)

        fout = open(self.out_cata, 'w')
        #obj_cata = open(self.obj_cata, 'r')  # The file contains the 
                                                     # objects of interest
        #pnames = obj_cata.readline().split() #The names of the parameters given 
                                             #in the first line in the obj_cata


        #print('Before multi')
        #pool = Pool(processes=1)
        #results = pool.map(self.P.main, obj_cata)
        #pool.close()
        #pool.join()


        #self.obj_cata = obj_cata
        #pnames = obj_file.readline().split() #The names of the parameters given 
                                             #in the first line in the obj_cata
        #self.obj_cata = obj_cata

        with open(self.obj_cata, 'r') as f_objs:
            obj_lines = f_objs.readlines()

        pnames = obj_lines[0].strip().split()
        obj_values = []
        for oval in obj_lines[1:]:
            if len(oval) == 1:
                continue
            oval2 = []
            for i in oval.strip().split():
                try:
                    oval2.append(float(i))
                except:
                    oval2.append(i)
        
            obj_values.append(oval2)

        print('pnames', pnames)
        #print(obj_values)
        #obj_values = np.genfromtxt(self.obj_cata, names=True)

        #pnames = obj_values.dtype.names
        #obj_values = list(obj_values[i] for i in pnames)
        #obj_values = np.column_stack(obj_values)                              
        #print(obj_values)

        #print('pnames', pnames)
        #print(obj_values)
        #print(obj_values.size)

        if len(obj_values) == 0:
            print('No lines in {}'.format(self.obj_cata))
            print('Exiting PyMorph')
            sys.exit()


        # writing a input catalogue (restart.cat) for failed objects
        f_failed = open(self.restart_file, 'w')
        for pname in pnames:
            f_failed.writelines(['{} '.format(pname)])
        f_failed.writelines(['flag \n'])
        f_failed.close()

        self.psfcounter = 0       #For getting psf in the case of unknown ra

        #print(1, imgsize, whtsize)

        #print('P5')
        P = Pipeline()
        #print('P6')
        P.pymorph_config = self.pymorph_config
        P.sex_config = self.sex_config

        P.center_deviated = self.center_deviated

        P.final_result_file = self.final_result_file
        P.restart_file = self.restart_file

        #It is to give all the parameter values to the instance P
        for k, v in self.pymorph_config.items():
            setattr(P, k, v)

        P.pnames = pnames
        P.params_to_write = params_to_write

        
        #If galcut then it passes to pipeline otherwise it reads the image
        if self.galcut:
            pass 
        else:
            self.imagedata, self.header0 = fitsio.read(self.imagefile, 
                                                       header=True)
            self.NXPTS = self.imagedata.shape[1]
            self.NYPTS = self.imagedata.shape[0]
            #print(1, 'self.NXPTS, self.NYPTS', self.NXPTS, self.NYPTS) 
            #sys.exit()
            #print(self.header0)
            P.IMG_HEADER, P.SEx_GAIN = get_header(self.header0)
            #print(P.IMG_HEADER) 
            print("Using large image. imagefile >>> {}".format(self.imagefile))

            if os.path.exists(self.whtfile):
                self.weightexists = True
                
                #XXX

                if re.search("rms", self.whtfile.lower()):
                    self.weightdata = fitsio.read(self.whtfile)
                    print("whtfile >>> {}".format(self.whtfile))
                elif re.search("var", self.whtfile.lower()):
                    self.weightdata = fitsio.read(self.whtfile)
                    self.weightdata = np.sqrt(self.weightdata)
                    print("whtfile >>> {}".format(self.whtfile))
                elif re.search("weight", self.whtfile.lower()):
                    self.weightdata = fitsio.read(self.whtfile)
                    self.weightdata = 1 / self.weightdata
                    print("whtfile >>> {}".format(self.whtfile))
                else:
                    print('Weight file is not understood. Please include ' + \
                          'the word weight/rms to the weight file name. ' + \
                          'If it is weight the rms will be found by 1/sqrt(w)')

                #XXX
                #self.weightdata = wht[0].read()

            elif self.whtfile == 'SDSS':
                self.weightdata = return_sigma(self.imagedata,
                                               gain=4.71, ncombine=1)  
            else:
                #print(4)
                #self.SEx_GAIN = 1.
                self.weightdata = np.ones(self.imagedata.shape)
                self.weightexists = False
                print('No weight image found\n')
            
            #except IOError as e:
            #    print(self.imagefile, "I/O error ({}): {}".format(e.args[0], e.args[1]))
            #    os._exit(0)

            #print(self.sex_cata)


            P.imagedata = self.imagedata
            P.weightdata = self.weightdata
            P.weightexists = self.weightexists
            P.NXPTS = self.NXPTS
            P.NYPTS = self.NYPTS

            if os.path.exists(self.sex_cata):
                pass
            else:
                print('The SExtractor catalogue for your frame is NOT found. ' \
                      'One is being made using the default values. It is always '\
                      'recommended to make SExtractor catalogue by YOURSELF as '\
                      'the pipeline keeps the sky value at the SExtractor value '\
                      'during the decomposition.')
                #print(6)
                PS = PySex(self.SEX_PATH)
                PS.RunSex(self.sex_config,  
                       self.imagefile, self.whtfile,
                       self.sex_cata, P.SEx_GAIN,
                       check_fits=None, sconfig='default')

        #Go to pipeline
        for obj_value in obj_values:
            try:
                P.main(obj_value)
            except:
                pass
        print(123)
        #sys.exit()
        #print(results)

    def find_and_fit_obj_cata(self):
        '''
        Find all the sextractor objects and find fits. It rewrites the
        obj_cata in config.ini
        '''

        #print('find_and_fit 1')
        #fc = open(self.obj_cata, 'w')
        #if self.redshift == 9999:
        #    fc.writelines(['gal_id ra dec mag\n'])
        #else:
        #    fc.writelines(['gal_id ra dec mag z\n'])

        values = np.genfromtxt(self.sex_cata)
        #print(self.sex_cata)

        #print(values.shape)
        con_remove = (values[:, 1] > self.remove_obj_boundary) & (values[:, 2] > self.remove_obj_boundary) & (values[:, 1] < self.NXPTS - self.remove_obj_boundary) & (values[:, 2] < self.NYPTS - self.remove_obj_boundary) 
        con = (values[:, 17] > self.maglim[1]) & (values[:, 17] < self.maglim[0]) & (values[:, 16] < self.stargal_prob_lim) & con_remove  

        #If we want only a few objects then you can use the below line by
        #changing the coordinate of objects
        #con = con & (values[:, 1] > 1839) & (values[:, 1] < 1840) & (values[:, 2] > 1101) & (values[:, 2] < 1102) 

        #con = (values[:, 16] < 0.8) #& (values[:, 16] > 0.65)
        #con = (values[:, 16] < 0.8) #& (values[:, 17] < 23)
         
        values = values[con]
        #print(values.shape)
        values = values[:, [0, 3, 4, 17]]
        values[:, 0] = values[:, 0].astype(int)
        #print(1, values)
        #values[:, 3] = values[:, 3] / 15.0

        if self.redshift != 9999:
            values[:, -1] = self.redshift
            np.savetxt(self.obj_cata, values, fmt='%i %.8f %.8f %.2f %.3f',
                    header='gal_id ra dec mag z', comments='')
        else:
            #print(2, values)
            np.savetxt(self.obj_cata, values, fmt='%i %.8f %.8f %.2f',
                    header='gal_id ra dec mag', comments='')


        searchrad = '0.05arc'

        #print('find_and_fit 2')


        #sys.exit()
    #The old function for psfselect = 2
    #    if psfselect == 2:
    #        center_deviated = 0
    #        starthandle = 0
    #        os.system('ds9 &')
    #        time.sleep(2)
    #        run_psfselect(imagefile, DATADIR, obj_cata, self.galcut)
    #        os.system('xpaset -p ds9 quit')
    #        psflist = '@psflist.list'
    #        find_and_fit()
    #        main()
    #        if crashhandler:
    #            starthandle = 1
    #            os.system('mv restart.cat CRASH.CAT')
    #            obj_cata = 'CRASH.CAT' 
    #            main()
    #New function for psfselect=2 non-interactive for webservice
        #print(self.psfselect)
        #print(self.obj_cata)
    def pymorph(self):


        # now change dir to the 
        thisdir = os.getcwd()
        print("Current directory is ", thisdir)
        print("Output directory is ", self.OUTDIR)

        os.chdir(self.DATADIR)

        #if self.repeat == False & self.galcut == False:

        if self.psfselect == 2:
            self.Interactive = 0
            self.center_deviated = 0
            self.starthandle = 0
            run_psfselect(self.imagefile, 
                          self.DATADIR, 
                          self.obj_cata, 
                          self.galcut)
            psflist = os.path.join(self.DATADIR, '@psflist.list')
            self.find_and_fit()
            self.main()
            if self.crashhandler:
                starthandle = 1
                os.system(f'mv {self.restart_file} CRASH.CAT')
                self.obj_cata = os.path.join(self.OUTDIR, 'CRASH.CAT') 
                self.main()

    #The old function for psfselect=1
        elif self.psfselect == 1:
            self.Interactive = 1
            os.system('ds9 &')
            time.sleep(2)
            run_psfselect(self.imagefile, 
                          self.DATADIR, 
                          self.obj_cata, 
                          self.galcut)
            os.system('xpaset -p ds9 quit')
    #new function for psfselect=1 non-interactive for webservice (abhishek rawat)
    #    elif psfselect == 1:
    #        Interactive = 0
    #        run_psfselect(imagefile, DATADIR, obj_cata, self.galcut)
            



        elif self.psfselect == 0:
            self.center_deviated = 0
            self.starthandle = 0
            #print(111)
            if self.findandfit == True:
                self.find_and_fit_obj_cata()
            #print(111)
            #print(self.sex_cata)
            self.main_thread()
            if self.crashhandler:
                self.starthandle = 1
                os.system(f'mv {self.restart_file} CRASH.CAT')
                self.obj_cata = os.path.join(self.OUTDIR, 'CRASH.CAT')
                self.main()

        os.chdir(thisdir)





if __name__ == '__main__':
    p = PyMorph()
    #p.set_sexparams()
    #p.set_limits()
    p.pymorph()
    #p.main()
