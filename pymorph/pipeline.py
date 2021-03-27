
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
from runsexfunc import PySex
from writehtmlfunc import WriteHtmlFunc
import psffunc as pf  
import mask_or_fit as mf
       
        



class Pipeline(object):


    def __init__(self, params):
        
        self.DATADIR = params[0]
        self.OUTDIR = params[1]
        self.imagefile = params[2]
        self.whtfile = params[3]
        self.sex_cat = params[4]
        self.clus_cat = params[5]
        self.out_cat = params[6]
        self.rootname = params[7]
        self.GALFIT_PATH = params[8]
        self.SEX_PATH = params[9]
        self.repeat = params[10]
        self.galcut = params[11]
        self.decompose = params[12]
        self.detail = params[13]
        self.galfit = params[14]
        self.cas = params[15]
        self.crashhandler = params[16]
        self.components = params[17]
        self.devauc = params[18]
        self.fitting = params[19]
        self.chi2sq = params[20]
        self.goodness = params[21]
        self.center_deviation = params[22]
        self.center_constrain = params[23]
        self.ReSize = params[24]
        self.VarSize = params[25]
        self.FracRad = params[26]
        self.Square = params[27]
        self.FixSize = params[28]
        self.searchrad = params[29]
        self.pixelscale = params[30]
        self.H0 = params[31]
        self.WM = params[32]
        self.WV = params[33]
        self.redshift = params[34]
        self.psfselect = params[35]
        self.stargal_prob = params[36]
        self.star_size = params[37]
        self.psflist = params[38]
        self.which_psf = params[39]
        self.area_obj = params[40]
        self.manual_mask = params[41]
        self.mask_reg = params[42]
        self.thresh_area = params[43]
        self.threshold = params[44]
        self.mag_zero = params[45]
        self.maglim = params[46]
        self.host = params[47]
        self.database = params[48]
        self.table = params[49]
        self.usr = params[50]
        self.pword = params[51]
        self.dbparams = params[52]
        self.galfitv = params[53]
        self.FirstCreateDB = params[54]
        self.FILTER = params[55]
        self.LMag = params[56]
        self.UMag = params[57]
        self.LN = params[58]
        self.UN = params[59]
        self.LRe = params[60]
        self.URe = params[61]
        self.LRd = params[62]
        self.URd = params[63]
        self.bdbox = params[64]
        self.bbox = params[65]
        self.dbox = params[66]
        self.devauc = params[67]
        self.avoidme = params[68]
        self.no_mask = params[69]
        self.norm_mask = params[70]

        self.verbose = False

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
                if self.verbose:
                    print('psflist', self.psflist)
                    print('counter', self.psfcounter)
                
                #pf.UpdatePsfRaDec(self.DATADIR, psffile)
                pf.UpdatePsfRaDec(self.DATADIR, self.psflist)
                distance = 9999
                self.psfcounter += 1
            else:
                if self.verbose:
                    print(self.DATADIR, self.psflist, self.which_psf, ra, dec)
                psffile, distance = pf.getpsf(self.DATADIR,
                                              self.psflist, self.which_psf,
                                              ra, dec)
                distance = distance * 60. * 60.

        return psffile, distance


    def MakeCutOut(self, xcntr, ycntr, 
                   alpha_j, delta_j, 
                   SizeX, SizeY, 
                   TX, TY,
                   gimg, wimg,
                   weightexists, ReSize):
        """

        Make cutout image. The xcntr, ycntr are like iraf.
        SizeX, SizeY are half size

        """
        print(1, 'xcntr, ycntr', xcntr, ycntr)
        ExceedSize = 0
        #All are floor to make the size even number
        xmin = np.floor(xcntr - SizeX)
        ymin = np.floor(ycntr - SizeY)
        xmax = np.floor(xcntr + SizeX)
        ymax = np.floor(ycntr + SizeY)
        if ReSize:
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
        else:
            cut_xcntr = (xmax - xmin) / 2. + np.modf(xcntr)[0]
            cut_ycntr = (ymax - ymin) / 2. + np.modf(ycntr)[0]

        print(2, 'xcntr, ycntr', xcntr, ycntr)
        print(2, 'cut_xcntr, cut_ycntr', cut_xcntr, cut_ycntr)
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

        #print(13, gimg)
        #print(14, data)
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







    def main(self, line_j):
        if self.verbose:
            print('line_j', line_j) 
        pdb = {}                        #The parameter dictionary
        #for line_j in obj_file:
        # declare the flag
        #line_j = obj_file.copy()
        self.flag = 0
        
        try:
            values = line_j.split()
        except Exception as e:
            print('No lines in {}'.format(self.clus_cat)) 
            print('Exiting PyMorph')
            sys.exit()

        k = 0
        for pname in self.pnames:
            pdb[pname] = values[k]
            k += 1
        if self.verbose:
            print('pdb', pdb)
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
        if self.verbose:
            print('fstring {}'.format(self.fstring))


        #gimg = os.path.join(self.DATADIR, 'I{}.fits'.format(self.fstring))
        #wimg = os.path.join(self.DATADIR, 'W{}.fits'.format(self.fstring))
        #oimg = os.path.join(self.OUTDIR, 'O{}.fits'.format(self.fstring))
        #mimg = os.path.join(self.OUTDIR, 'M{}.fits'.format(self.fstring))
        #emimg = os.path.join(self.OUTDIR, 'EM_{}.fits'.format(self.fstring))

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
                gimg = None

        try:
            wimg = pdb["wimg"]   #Weight cut
        except:
            wimg = None
            print('Search for weight image (wimg) in galfit config')


        try:
            alpha_j = pdb["ra"]
            if len(alpha_j.split(':')) == 1:
                alpha_j = float(alpha_j)
                alpha_j *= 15
            else:
                h, m, s = alpha_j.split(':')
                print('alpha_j', h, m, s)
                alpha_j = ut.HMSToDeg(int(h), int(m), float(s)) 
        except:
            print("No ra is given")
            alpha_j = -9999

        try:
            delta_j = pdb["dec"]
            if len(delta_j.split(':')) == 1:
                delta_j = float(delta_j)
            else:
                d, m, s = delta_j.split(':')
                delta_j = ut.DMSToDeg(int(d), int(m), float(s))
        except:
            print("No dec is given")
            delta_j = -9999

        print('alpha_j, delta_j', alpha_j, delta_j)

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
                self.weightexists = True
            else:
                wimg = None
                self.weightexists = False

            print('Using cutouts')

        if self.verbose:
            print('gimg, self.weightexists', gimg, self.weightexists)
        TX = self.imagedata.shape[0]
        TY = self.imagedata.shape[1]

        print('TX TY', TX, TY)
        if self.ReSize:
            gimg = 'I{}.fits'.format(self.fstring)
            wimg = 'W{}.fits'.format(self.fstring)
        print('F1')
        #The sextractor runs on the cutout before resizing to estimate 
        #shallow sky
       
        PS = PySex(c.SEX_PATH)

        shallow_cat = os.path.join(self.DATADIR,
                                  '{}.Shallow'.format(self.sex_cat))
        if self.galcut == True:   
            #Given galaxy cutouts if the user provides sextractor catalogue
            #then it will not run SExtractor else do!
            try:
                self.sex_cat = os.path.join(self.DATADIR, 
                                            'sex_{}.cat'.format(self.fstring))
                check_fits = 'check{}.fits'.format(self.fstring),
                PS.RunSex(self.sex_params, 
                          os.path.join(self.DATADIR, gimg),
                          os.path.join(self.DATADIR, wimg),
                          self.sex_cat, self.SEx_GAIN, 
                          check_fits=check_fits,
                          sconfig='default')


                shallow_cat = os.path.join(self.DATADIR, 
                                      '{}.Shallow'.format(self.fstring))
                check_shallow_fits = 'check{}_shallow.fits'.format(self.fstring)
                PS.RunSex(self.sex_params,
                          os.path.join(self.DATADIR, gimg),
                          os.path.join(self.DATADIR, wimg),
                          shallow_cat, self.SEx_GAIN,
                          check_fits=check_shallow_fits,
                          sconfig='shallow')

            except Exception as e:
                print(type(e))     # the exception instance
                print(e.args)      # arguments stored in\
                                     # .args
                print(e)           # __str__ allows args\
                                     # to printed directly
                print("Something bad happened (Sextractor)!!!!\n\n")
                print(traceback.print_exc())
            #sys.exit() 
        #Lets see whether decorator can be used with ximg and yimg
        try:
            ximg = float(pdb["ximg"])
            if self.ReSize and self.galcut and self.repeat:
                ximg = TX / 2.0
        except:
            if self.verbose:
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
            if self.verbose:
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

        print('F2')
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
        elif self.searchrad.endswith('pix'):
            SeaPix = float(self.searchrad[:-3])
            SeaDeg = pixelscale * SeaPix  / (60.0 * 60.0)
        
        print('F3')
        # first count the number of "potential" targets in the search radius
        SexTargets = 0
        good_objects = []
        bad_objects = []
        good_distance = []
        bad_distance = []
        good_object = ''

        os.system('cp {} sex_{}.txt'.format(self.sex_cat, self.fstring))
        center_distance = 9999 #the distance from the center to the best target

        print('F4', gimg)
        if alpha_j == 9999:
            PS.get_sexobj(self.sex_cat, xpix, ypix, dmin)
        else:
            PS.get_sexobj(self.sex_cat, alpha_j, delta_j, dmin)

            self.flag = SetFlag(self.flag, GetFlag('NO_TARGET'))  

        #sys.exit()
        print('F5')
        if 1:#self.verbose:
            print(good_object)
        # now fit best object            
        #XXX
        if 1:#try:
            sex_id = PS.sex_num
            print("SExtractor ID >>> {}".format(sex_id))
            print(PS.alpha_s, PS.delta_s, PS.xcntr, PS.ycntr) 
            SexHalfRad = PS.radius #Sex halfrad 
            SexPosAng = PS.pos_ang
            pos_ang = PS.pos_ang
            axis_rat = PS.bbya #axis ration b/a
            eg = PS.eg
            if eg <= 0.05:
                eg = 0.07
            major_axis = PS.maj_axis

            SexMagAuto = PS.mag
            SexMagAutoErr = PS.mag_err
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
            # The shallow sky from the shallow run. If no shallow
            # sky it uses the deep sky.
            try:
                SexSky = PS.get_sexvalue(sex_id, shallow_cat, param='sky')
            except:
                SexSky = float(values[10])

            if not np.isnan(UserGivenSky):
                skyval = UserGivenSky #using the sky supplied by the user
            else:
                skyval = SexSky #use the sky value from sextractor
                

            if (position == 1) & (alpha_j == -9999):
                alpha_j = alpha_s
                delta_j = delta_s

            print(alpha_j, delta_j, xcntr, ycntr)

            ut.WriteError('\n\n###########   {} ###########\n'.format(gal_id))
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

            print('TX TY', TX, TY)
            SizeX, SizeY = self.FindCutSize(TX / 2, TY / 2,
                                            SexHalfRad, 
                                            pos_ang, axis_rat) 
            print('Calculated half sizes {}, {}'.format(SizeX, SizeY))

            print('F6', gimg, alpha_j, delta_j)
            # Finding psf and the distance between psf and image
            if self.decompose:
                psffile, distance = self.HandlePsf(UserGivenPsf, 
                                                   alpha_j, delta_j)
            else:
                psffile, distance = 'None', 9999

            print('psffile: {}, distance: {}'.format(psffile, distance))
            print('F7') 
            #sys.exit()
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
                             self.fstring, iseg_fits, oseg_cat, NoMask=False)

            else:
                if gimg is None:
                    gimg = os.path.join(self.DATADIR, 
                                        'I{}.fits'.format(self.fstring))
                    wimg = os.path.join(self.DATADIR, 
                                        'W{}.fits'.format(self.fstring))
                
                if 1: #self.galcut:
                #if self.ReSize: 
                    #Can take a function
                    #Can check whether a decorator can be used
                    if os.path.exists(gimg):
                        ut.WriteError('The file {} exists\n'.format(gimg))
                        run = 0
                        #break #Breaking the sextractor loop
                    try:
                        cut_params = self.MakeCutOut(xcntr, ycntr,
                                                     alpha_j, delta_j,
                                                     SizeX, SizeY,
                                                     TX, TY,
                                                     gimg, wimg,
                                                     self.weightexists,
                                                     self.ReSize)
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
                        #break
                #else:
                #    cut_xcntr = xcntr
                #    cut_ycntr = ycntr 
                #    SizeX = TX
                #    SizeY = TY
                #    ExceedSize = 0

                
                #else:
                #    print('F8', gimg)
                #    #gimg = 'I{}.fits'.format(self.fstring)
                #    #wimg = 'W{}.fits'.format(self.fstring)
                #    gimg = os.path.join(self.DATADIR,
                #                        'I{}.fits'.format(self.fstring))
                #    wimg = os.path.join(self.DATADIR,
                #                        'W{}.fits'.format(self.fstring))
                #    #Can take a function
                #    if os.path.exists(gimg):
                #        ut.WriteError('The file {} exists\n'.format(gimg))
                #        run = 0
                #        #break #Breaking the sextractor loop

                #    # Sizes are total size
                #    print('F9', gimg)
                #    #XXX
                #    if 1:#try: 
                #        #print(xcntr, ycntr, alpha_j, delta_j,SizeX, SizeY,
                #        #                           TX, TY)
                #        #xcntr = SizeX
                #        #ycntr = SizeY
                #        cut_params = self.MakeCutOut(xcntr, ycntr, 
                #                                   alpha_j, delta_j, 
                #                                   SizeX, SizeY, 
                #                                   TX, TY, 
                #                                   gimg, wimg,
                #                                   self.weightexists,
                #                                   self.ReSize)
                #        cut_xcntr = cut_params[0]
                #        cut_ycntr = cut_params[1]
                #        SizeX = cut_params[2]
                #        SizeY = cut_params[3]
                #        ExceedSize = cut_params[4]
                #        print('cut_xcntr, cut_ycntr', cut_xcntr, cut_ycntr)
                #        #sys.exit()
                #    #except  Exception as e:
                #    else:
                #        print(type(e))     # the exception instance
                #        print(e.args)      # arguments stored in\
                #                             # .args
                #        print(e)           # __str__ allows args\
                #                             # to printed directly
                #        ut.WriteError('Cutout exists!')
                #        print(traceback.print_exc())
                #        break

                print('Center of cutimage and exceed size ', \
                      cut_xcntr, cut_ycntr, ExceedSize)
                print('Full Sizes ', SizeX, SizeY)
                good_object[1] = cut_xcntr
                good_object[2] = cut_ycntr
                print('Len good_object', len(good_object))
                print('F10')
                if self.galcut & self.ReSize == 0:
                    pass
                elif ExceedSize:
                    self.flag = SetFlag(self.flag, GetFlag('EXCEED_SIZE'))

                # Runs sextractor to find the segmentation map
                ##RunSegSex(os.path.join(self.DATADIR, gimg))

                # This creates ellipse mask and used for 
                # ellipse fitting and casgm
                                
                
                seg_cat = os.path.join(self.DATADIR,
                                        '{}.cat.seg'.format(self.fstring))
                check_seg_fits = os.path.join(self.OUTDIR, 
                                        'seg_{}.fits'.format(self.fstring))
                if self.verbose:
                    print('gimg, check_seg, fits', gimg, check_seg_fits)

                if os.path.exists(seg_cat):
                    os.remove(seg_cat)

                PS.RunSex(self.sex_params, gimg, wimg, seg_cat,
                          self.SEx_GAIN, check_fits=check_seg_fits,
                          sconfig='seg')

                emimg = os.path.join(self.OUTDIR, 'EM_{}.fits'.format(self.fstring))
                mimg = os.path.join(self.OUTDIR, 'M_{}.fits'.format(self.fstring))
                #sys.exit()
                EM = ElliMaskFunc(emimg, cut_xcntr, cut_ycntr, 1,
                                  center_limit=5., seg_limit=1e-5)
                EM.emask(check_seg_fits, seg_cat, self.fstring)


                # Fitting ellipse task or the manual 1d finder
                print('F11')
                if self.verbose:
                    print('gimg', gimg)
                if self.decompose:
                    FE = ut.FindEllipse(xcntr, ycntr, SexHalfRad,
                                        SexPosAng, axis_rat, skyval,
                                        self.fstring)
                    FE.profile(gimg, output=False, emimg=emimg)

                    #Need XXX
                    print('cut_xcntr, cut_ycnt', cut_xcntr, cut_ycntr)
                    print('Len good_object', len(good_object)) 
                    MF = MaskFunc(mimg, 
                                  cut_xcntr, cut_ycntr, 
                                  SizeX, SizeY, 
                                  good_object)
                    MF.gmask(self.threshold, self.thresh_area, 
                             self.fstring, check_seg_fits, seg_cat,
                             avoidme=self.avoidme, NoMask=self.no_mask,
                             seg_limit=1e-5)
                    #sys.exit()
                    galfit_conf = 'G_{}.in'.format(self.fstring)
                    oimg = 'O_{}.fits'.format(self.fstring)

                    print('F12')

                    print('Len good_object', len(good_object)) 
                    CF = GalfitConfigFunc(self.DATADIR,
                                          gimg, wimg,  
                                          cut_xcntr, cut_ycntr, 
                                          SizeX, SizeY,
                                          self.components, self.fitting,
                                          psffile,
                                          good_object, 
                                          self.sex_cat,
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
                    print('F13') 
            # Estimates sky parameters
            #XXX
            if 1:
                sky_values = FindYetSky(self.fstring, 
                                        self.sex_params, self.SEX_PATH, 
                                        os.path.join(self.DATADIR, gimg),
                                        os.path.join(self.DATADIR, wimg),
                                        self.sex_cat,
                                        cut_xcntr, cut_ycntr,
                                        check_seg_fits,
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
            file_fit = os.path.join(self.OUTDIR, 'fit.log')
            file_fit2 = os.path.join(self.OUTDIR, 'fit2.log')
            if os.access(file_fit, os.F_OK):
                os.remove(file_fit)

            print('F14') 
            if self.decompose:
                try:
                    if self.galfit & self.detail:
                        ConfigIter(gimg, wimg, 
                                   cut_xcntr, cut_ycntr, 
                                   SizeX, SizeY, 
                                   good_object, 
                                   psffile, z)
                    elif self.galfit:
                        cmd = '{} {} '.format(self.GALFIT_PATH, 
                                                         galfit_conf)
                        print(cmd)
                        os.system(cmd)
                        #sys.exit()
                        f_fit2 = open(file_fit2, 'a')

                        if os.path.exists(file_fit):
                            for line in open(file_fit, 'r'):
                                f_fit2.writelines([str(line)])
                        f_fit2.close()
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
                
                print('F15')
                
                oimg = os.path.join(self.OUTDIR, 
                                    'O_{}.fits'.format(self.fstring))
                if os.path.exists(oimg):
                    FE.profile(oimg, output=False, emimg=emimg)
                else:
                    print("Error in finding output elliptical profile\n\n")


            if os.path.exists(mimg):

                try:
                    if os.access('P_' + self.fstring + '.png', os.F_OK):	
                        os.remove('P_' + self.fstring + '.png')
                
                    PF = PlotFunc(oimg, mimg, self.fstring, 
                                  SexSky, SkySig, self.mag_zero,
                                  save_name=None)
                    PF.plot_profile()
                    goodness = 9999 #PF.plot_profile().goodness
                except Exception as e:
                    print(type(e))     # the exception instance
                    print(e.args)      # arguments stored in\
                                             # .args
                    print(e)           # __str__ allows args\
                                           # to printed directly

                    goodness = -9999
                    ut.WriteError('Error in plotting \n')
                    run = 0	
                    self.flag = SetFlag(self.flag, GetFlag('PLOT_FAIL'))
                
            else:
                ut.WriteError('Could not find Mask image for plottong \n')

            print('F16')

            galfit_failed = []
            f_cat_results = []

            if isset(self.flag, GetFlag("GALFIT_FAIL")): #or \
               #isset(flag, GetFlag("LARGE_CHISQ")) or \
               #isset(flag, GetFlag("FAKE_CNTR")) or \
               #isset(flag, GetFlag("BULGE_AT_LIMIT")) or \
               #isset(flag, GetFlag("DISK_AT_LIMIT")):
                failedvalues = line_j.split()
                ###for fv in failedvalues:
                ###    self.f_failed.writelines([' '.format(fv)])
                ###self.f_failed.writelines(['{}\n'.format(self.flag)])
                galfit_failed = failedvalues.append(self.flag)
            
            #print(1, gal_id)
            #print(2, good_object) 
            #for i in (gal_id + good_object).split(' '):
            #    if len(i) == 0:
            #        pass
            #    else:
            #        f_cat_results.append(i)

            f_cat_results.extend(good_object)

            #self.f_cat.writelines(['{} '.format(gal_id)])
            #self.f_cat.write(good_object)


            print('F17')

            try:
                WF = WriteHtmlFunc(self.fstring,
                                   cut_xcntr, cut_ycntr, 
                                   alpha_j, delta_j, 
                                   distance, z,
                                   SexMagAuto, SexMagAutoErr, 
                                   SexTargets, SexSky,
                                   self.flag, SexHalfRad, self.mag_zero,
                                   C, C_err, A, A_err, S, S_err, G, M, 
                                   self.components, self.decompose, 
                                   self.repeat,
                                   goodness, self.EXPTIME,
                                   self.H0, self.WM, self.WV, self.pixelscale)
                WF.writeparams(self.params_to_write)
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
            if os.access(self.sex_cat, os.F_OK):
                os.remove(self.sex_cat)

        #except Exception as e:
        #    print(type(e))     # the exception instance
        #    print(e.args)      # arguments stored in\
            # .args
        #    print(e)           # __str__ allows args\
            # to printed directly
        #    print("something bad happened (general object search)!!!!\n\n")
        #    print(traceback.print_exc())


        print('F18')
    
        return f_cat_results, galfit_failed





