
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
from runsexfunc import PySex
#from writehtmlfunc import WriteParams
import psffunc as pf  
import mask_or_fit as mf
       
        


class ReturnClass(object):
    
    def __init__(self):
        pass
    
    def _find_cutout_size(self): 
        """Return the size of the cutout. SizeXX, SizeYY are the size of 
          cutout if c.galcut is True"""

        print(self.pos_ang)
        self.pos_ang = self.pos_ang * np.pi / 180. #pa in radian
        SizeX = self.SexHalfRad * self.FracRad * np.abs(np.cos(self.pos_ang)) + \
                self.axis_rat * self.SexHalfRad * self.FracRad * np.abs(np.sin(self.pos_ang))
        SizeY = self.SexHalfRad * self.FracRad * np.abs(np.sin(self.pos_ang)) + \
                self.axis_rat * self.SexHalfRad * self.FracRad * np.abs(np.cos(self.pos_ang))
        SizeX = int(SizeX)
        SizeY = int(SizeY)

        if self.Square:
            SizeX = max(SizeX, SizeY)
            SizeY = SizeX
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

    def _distance_psf_obj(self):
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

            distance = 3600.0 * np.sqrt((dec - dec_p)**2.0 + \
                       ((ra - ra_p) * np.cos(dec * ut.Get_R()))**2.0)
        except:
            distance = 9999
            print('Distance between psf and image is 9999')

        return distance

    def _handle_psf(self, UserGivenPsf, ra, dec):

        """Determine the psf used for fitting"""

        if self.repeat:
            if np.abs(ra) == 9999 or np.abs(dec) == 9999:
                distance = 9999
                psffile = c.pfile
            else:
                psffile = c.pfile
                pf.UpdatePsfRaDec(c.pfile)
                distance = super()._distance_psf_obj(c.pfile, ra, dec)

        else:
            if UserGivenPsf != 'None':
                psffile = UserGivenPsf
                pf.UpdatePsfRaDec(self.DATADIR, psffile)
                distance = super()._distance_psf_obj(psffile, ra, dec)
            elif np.abs(ra) == 9999 or np.abs(dec) == 9999:
                psffile = self.psflist[self.psfcounter]
                pf.UpdatePsfRaDec(self.DATADIR, psffile)
                distance = 9999
                self.psfcounter += 1
            else:
                psffile, distance = pf.getpsf(self.DATADIR,
                                              self.psflist, self.which_psf,
                                              ra, dec)
                distance = distance * 60. * 60.

        return psffile, distance


    def _make_cutout(self):
        """

        Make cutout image. The xcntr, ycntr are like iraf.
        SizeX, SizeY are half size

        """
        print(9, self.xcntr, self.SizeX, self.ycntr, self.SizeY)
        ExceedSize = 0
        #All are floor to make the size even number
        xmin = np.floor(self.xcntr - self.SizeX)
        ymin = np.floor(self.ycntr - self.SizeY)
        xmax = np.floor(self.xcntr + self.SizeX)
        ymax = np.floor(self.ycntr + self.SizeY)
        if xmin < 0:
            xmin = 0
            cut_xcntr = self.xcntr
            self.ExceedSize = 1
        else:
            cut_xcntr = self.SizeX + np.modf(self.xcntr)[0]
        if ymin < 0:
            cut_ycntr = self.ycntr
            ymin = 0
            self.ExceedSize = 1
        else:
            cut_ycntr = self.SizeY + np.modf(self.ycntr)[0]
        if xmax > self.TX - 1:
            xmax = self.TX
            self.ExceedSize = 1
        if ymax > self.TY - 1:
            ymax = self.TY
            self.ExceedSize = 1

        ymin = int(ymin)
        ymax = int(ymax)
        xmin = int(xmin)
        xmax = int(xmax)
        # FIX check c.imagedata or c.ggimage
        print(10, ymin, ymax, xmin, xmax)
        data = self.imagedata[ymin:ymax, xmin:xmax]
        print(100, self.imagedata.shape)
        # END
        self.SizeY, self.SizeX = data.shape

        hdict = {}
        try:
            hdict['RA_TARG'] = self.alpha_j
            hdict['DEC_TARG'] = self.delta_j
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

        print(13, self.gimg)
        print(14, data)
        fits = fitsio.FITS(os.path.join(self.DATADIR, self.gimg), 'rw')
        fits.write(data, header=hdict)
        fits.close()


        # FIX
        #Making weight image cut
        if self.weightexists:
            data_wht = self.weightdata[ymin:ymax,xmin:xmax].copy()

            fits = fitsio.FITS(os.path.join(self.DATADIR, self.wimg), 'rw')
            fits.write(data_wht)
            fits.close()
        else:
            print('Cannot creat weight image. If you supply weight image \
                   please check whether it exists or report a bug')
        #END
        return cut_xcntr, cut_ycntr




class Pipeline(ReturnClass):
    
    __version__ = 3.0

    def __init__(self, config_file='config.ini'):

        pass
        #self.SP = SexParams(self.sex_params, self.mag_zero)


    def _target_initialize(self, obj_value, pdb, pnames):


        k = 0
        for pname in pnames:
            print(pname)
            pdb[pname] = obj_value[k]
            k += 1

        try:
            self.gal_id = pdb["gal_id"]
        except:
            print("No gal_id. Try to use gimg")

            try:
                self.gal_id = pdb["gimg"].split('.')[0] #gal_id will be 
                                                   #filename without .fits
            except:
                print("No image or gal_id found in the object catalogue." \
                      "Exiting")
                os._exit(0)

        self.fstring = '{}_{}'.format(self.rootname, self.gal_id)
        print('fstring {}'.format(self.fstring))


        try:
            self.gimg = pdb["gimg"]    #Galaxy cutout
        except:
            print("No gimg given.")
            if os.path.exists(os.path.join(self.DATADIR, 
                             'I{}.fits'.format(self.fstring))):
                self.gimg = 'I{}.fits'.format(self.fstring)
            elif os.path.exists(os.path.join(self.DATADIR, 
                               '{}.fits'.format(self.gal_id))): 
                self.gimg = '{}.fits'.format(self.gal_id)
            else:
                print("No possible gimg found")
                self.gimg = 'None'

        try:
            self.wimg = pdb["wimg"]   #Weight cut
        except:
            self.wimg = 'None'
            print('Search for weight image (wimg) in galfit config')


        if isinstance(pdb["ra"], float): 
            self.alpha_j = pdb["ra"]     
        elif isinstance(pdb["ra"], str):
            h, m, s = pdb["ra"].split(':')
            self.alpha_j = ut.HMSToDeg(int(h), int(m), float(s)) 
        else:
            print("No ra is given")
            self.alpha_j = -9999

        if isinstance(pdb["dec"], float):
            self.delta_j = pdb["dec"]
        elif isinstance(pdb["dec"], str):
            d, m, s = delta.split(':')
            self.delta_j = ut.DMSToDeg(int(d), int(m), float(s))
        else:
            print("No dec is given")
            self.delta_j = -9999

        #print(self.alpha_j, self.delta_j)
        try:
            self.z = float(pdb["z"])
        except:
            print("No z is given")
            self.z = 9999
        
        #1000
        #Paste from modes.py
        #2000
    #3000
    #Paste from modes.py
    #4000
            
    def _search_radius(self):
        
        if self.searchrad is None:
            if self.position:
                self.searchrad = '1arc'
                print('No search radius found. Setting to 1 arc sec')
            else:
                self.searchrad = '10pix'
                print('No search radius found. Setting to 10 pix')

        if self.searchrad.endswith('arc'):
            self.SeaDeg = float(self.searchrad[:-3]) / (60.0 * 60.0)
            self.SeaPix = 10.0
            print('Arcsec')
        elif self.searchrad.endswith('pix'):
            self.SeaPix = float(self.searchrad[:-3])
            self.SeaDeg = pixelscale * self.SeaPix  / (60.0 * 60.0)
            
        print('self.SeaPix, self.SeaDeg', self.SeaPix, self.SeaDeg)

            
    def _get_shallow_sky(self, target_sky):
        
        # The shallow sky from the shallow run. If no shallow
        # sky it uses the deep sky.

        file_shallow = '{}.Shallow'.format(self.sex_cata)
        if os.path.exists(file_shallow):
            v_shallow = np.genfromtxt(file_shallow)
            try:
                con = (v_shallow[:, 19] == sex_id)
                
                self.SexSky = v_shallow[:, 10][con]
            except:
                pass
        else:
            self.SexSky = target_sky
            
        print('SexSky', self.SexSky)

        if self.UserGivenSky is not None:
            self.SexSky = self.UserGivenSky #using the sky supplied by the user
        #else:
        #    self.SexSky = self.SexSky #use the sky value from sextractor



    def _flags_initialization(self):
        # Adding initial setup to the flag
        if self.repeat:
            self.flag = SetFlag(self.flag, GetFlag('REPEAT'))
        if self.fitting[0] and 'bulge' in self.components:
            self.flag = SetFlag(self.flag, GetFlag('FIT_BULGE_CNTR'))
        if self.fitting[1] and 'disk' in self.components:
            self.flag = SetFlag(self.flag, GetFlag('FIT_DISK_CNTR'))
        if self.fitting[2]:
            self.flag = SetFlag(self.flag, GetFlag('FIT_SKY'))


    def _galcut_images(self):
        '''
        Useful for galcut mode. Here I am planning to give gimg in obj_cata 
        If ReSize is True then it will create the usual file and weight names
        '''

        print('Image is >>> {}'.format(self.gimg))

        gfits = fitsio.FITS(os.path.join(self.DATADIR, self.gimg), 'r')
        self.imagedata = gfits[0].read()
        self.header0 = gfits[0].read_header()
        gfits.close()

        gheader = ut.CheckHeader(self.header0) #Will set up global header parameters
        self.EXPTIME = gheader[0]
        self.RDNOISE = gheader[1], 
        self.GAIN = gheader[2]
        self.SEx_GAIN = gheader[3]
        self.NCOMBINE = gheader[4]

        if os.path.exists(os.path.join(self.DATADIR, self.wimg)):
            wfits = fitsio.FITS(os.path.join(self.DATADIR, self.wimg))
            self.weightdata = wfits[0].read()
            wfits.close()
            self.weightexists = True
        else:
            self.wimg = None
            self.weightexists = False

        print('Using cutouts')

        print(1, self.gimg)
        self.TX = self.imagedata.shape[0]
        self.TY = self.imagedata.shape[1]

        if self.ReSize:
            self.gimg = 'I{}.fits'.format(self.fstring)
            self.wimg = 'W{}.fits'.format(self.fstring)

        #The sextractor runs on the cutout before resizing to estimate 
        #shallow sky
    ##if self.galcut == True:   
        #Given galaxy cutouts if the user provides sextractor catalogue
        #then it will not run SExtractor else do!
        PS = PySex(c.SEX_PATH)
        if not os.path.exists(self.sex_cata):
            if 1:

                PS.RunSex(self.sex_params, 
                          os.path.join(self.DATADIR, self.gimg),
                          os.path.join(self.DATADIR, self.wimg),
                          os.path.join(self.DATADIR, self.sex_cata),
                          self.SEx_GAIN, check_fits=None, 
                          sconfig='default')

                PS.RunSex(self.sex_params,
                          os.path.join(self.DATADIR, self.gimg),
                          os.path.join(self.DATADIR, self.wimg),
                          os.path.join(self.DATADIR, 
                                      '.Shallow'.format(self.sex_cata)),
                          self.SEx_GAIN, check_fits=None, 
                          sconfig='shallow')

#             except Exception as e:
#                 print(type(e))     # the exception instance
#                 print(e.args)      # arguments stored in\
#                                      # .args
#                 print(e)           # __str__ allows args\
#                                      # to printed directly
#                 print("Something bad happened (Sextractor)!!!!\n\n")
#                 print(traceback.print_exc())
#         sys.exit() 
    
    def _galcut_cutout(self):
        
        if self.ReSize: 
            #Can take a function
            #Can check whether a decorator can be used
            if os.path.exists(self.gimg):
                ut.WriteError('The file {} exists\n'.format(self.gimg))
                run = 0
                #break #Breaking the sextractor loop
            try:
                cut_params = self._make_cutout(self.TX, self.TY,
                                               self.gimg, self.wimg,
                                               self.weightexists)
                cut_xcntr = cut_params[0]
                cut_ycntr = cut_params[1]
                self.ExceedSize = cut_params[4]
            except  Exception as e:
                print(type(e))     # the exception instance
                print(e.args)      # arguments stored in\
                                     # .args
                print(e)           # __str__ allows args\
                                     # to printed directly
                ut.WriteError('Cutout exists!')
                print(traceback.print_exc())
                #break
        else:
            cut_xcntr = xcntr
            cut_ycntr = ycntr 
            SizeX = self.TX
            SizeY = self.TY
            self.ExceedSize = 0
        
        return cut_xcntr, cut_ycntr, SizeX, SizeY
        ##return run
                            
            
    def _not_galcut_cutout(self):
        
        print(11, self.gimg)
        self.gimg = 'I{}.fits'.format(self.fstring)
        self.wimg = 'W{}.fits'.format(self.fstring)
        #Can take a function
        if os.path.exists(self.gimg):
            ut.WriteError('The file {} exists\n'.format(self.gimg))
            run = 0
            #break #Breaking the sextractor loop
            print(3)

        # Sizes are total size
        print(12, self.gimg)
        #XXX
        if 1:#try: 
            cut_params = self._make_cutout()
            cut_xcntr = cut_params[0]
            cut_ycntr = cut_params[1]

        #except  Exception as e:
        else:
            print(type(e))     # the exception instance
            print(e.args)      # arguments stored in\
                                 # .args
            print(e)           # __str__ allows args\
                                 # to printed directly
            ut.WriteError('Cutout exists!')
            print(traceback.print_exc())
            #break
        
        return cut_xcntr, cut_ycntr

    def _get_ximg_yimg(self, pdb):
           
        '''
        It will be useful for mode galcut
        '''
        try:
            ximg = float(pdb["ximg"])
            if self.ReSize and self.galcut and self.repeat:
                ximg = self.TX / 2.0
        except:
            print('No ximg is given or either ReSize or galcut or repeat '\
                  'keywords are False. Trying to find from the cutout if '\
                  'no ra dec in the image header')
            if self.galcut == True & self.position == 0:
                ximg = self.TX / 2.0
            else:
                ximg = -9999

        try:
            yimg = float(pdb["yimg"])
            if self.ReSize and self.galcut and self.repeat:
                yimg = self.TY / 2.0
        except:
            print('No yimg is given or either ReSize or galcut or repeat '\
                  'keywords are False. Trying to find from the cutout if '\
                  'no ra dec in the image header')
            if self.galcut == True & self.position == 0:
                yimg = self.TY / 2.0
            else:
                yimg = -9999
                
        return ximg, yimg

    
    def _get_bkg_center(self, pdb):
            
        try:
            bxcntr = float(pdb["bxcntr"])
        except:
            bxcntr = -9999

        try:
            bycntr = float(pdb["bycntr"])
        except:
            bycntr = -9999
                
        return bxcntr, bycntr
 
    def _potential_target(self):
           
        '''
        It finds the potential target. It uses sextractor catalog and first 
        uses the sky coordinate to find any objects at a user given distance.
        Then it uses the image coordinates to find the objects. If there is 
        sky objects in the obj_cata it will find by comparing the sextractor
        catalog. If there is only image coordinate then it will find the
        objects in the sextractor image coordinate system

        This is used in the case of find and fit, a few objects in large frame
        and galcut images
        '''

        SexTargets = 0
        good_objects = []
        bad_objects = []
        good_distance = []
        bad_distance = []
        good_object = None

        subprocess.call(['cp', '{}'.format(self.sex_cata), 
                         'sex_{}.txt'.format(self.fstring)])

        print(2, self.gimg)

        values_sex = np.genfromtxt(self.sex_cata)
        print(self.sex_cata)
        alpha_s = values_sex[:, 3]
        delta_s = values_sex[:, 4]
        if self.position == 0:
            alpha_s = np.full_like(values_sex.shape[0], 9999)
            delta_s = np.full_like(values_sex.shape[0], 9999)
        sex_id = values_sex[0]
        xcntr  = values_sex[:, 1]
        ycntr  = values_sex[:, 2]

        print(self.alpha_j, alpha_s)
        print(self.SeaDeg, self.SeaPix)
        
        if self.alpha_j is not None:
            con = (abs(self.alpha_j - alpha_s) < self.SeaDeg) & \
            (abs(delta_s - self.delta_j) < self.SeaDeg)
            curr_distance = np.sqrt((self.alpha_j - alpha_s)**2 + (delta_s - self.delta_j)**2)
            
        else:
            con = (abs(xcntr - self.ximg) < self.SeaPix) & (abs(ycntr - self.yimg) < self.SeaPix)
            curr_distance = np.sqrt((xcntr - self.ximg)**2 + (ycntr - self.yimg)**2)

        curr_distance = curr_distance[con]
        values_sex = values_sex[con] 
        
        #print(curr_distance)
        #print(values_sex)
        
#         print(curr_distance_img, curr_distance_sky)
#         curr_distance = np.min(curr_distance_img, curr_distance_sky)
#         print(curr_distance.shape)
#         curr_distance = curr_distance[con]
        for cud in curr_distance:
            print("Candidate distance: {:3f}".format(cud)) 
        SexTargets = curr_distance.shape[0]

        if SexTargets > 0:
            center_distance = np.min(curr_distance)
            print("New Preferred target!!")
            good_object = values_sex[np.argmin(curr_distance)]
            print('good_object', good_object.shape)

            print("Target distance: {:.3f}".format(center_distance))
        else:
            # No suitable target found
            print("No Target Found!!!!")
            good_object = np.full(19, 9999)
            good_object[11] = 0
            self.flag = SetFlag(self.flag, GetFlag('NO_TARGET'))                 #except:

        #print(good_object)

        
        target = mf.GetSExObj(values=good_object)

        sex_id = target.sex_num
        print("SExtractor ID >>> {}".format(sex_id))


        self.xcntr, self.ycntr = target.sex_center
        print('sex_xcntr', self.xcntr)
        self.SexHalfRad = target.radius #Sex halfrad 
        self.SexPosAng = target.pos_ang
        self.pos_ang = target.pos_ang
        self.axis_rat = target.bbya #axis ration b/a
        self.eg = target.eg
        if self.eg <= 0.05:
            self.eg = 0.07
        self.major_axis = target.maj_axis

        self.SexMagAuto = target.mag
        self.SexMagAutoErr = target.mag_err
        self.UMag = -1e5
        self.LMag = -1e5
        self.URe = 1e0
        self.URd = 1e0
        if self.UMag > -9999.0:
            self.UMag = self.SexMagAuto - 7.0
        if self.LMag < 9999.0:
            self.LMag = self.SexMagAuto + 7.0

        # FIX
        if self.URe < 9999.0:
            self.URe = self.SexHalfRad * 50.0
        if self.URd < 9999.0:
            self.URd = self.SexHalfRad * 50.0
        # END
        print(1)

        return good_object
                
    def main(self, obj_value):
        
        pdb = {}                        #The parameter dictionary
        #for obj_value in obj_values:
        print(1, obj_value)
        # declare the flag
        self.flag = 0

        print(11)
        self._target_initialize(obj_value, pdb, self.pnames)
        print(22)
        print(2, self.alpha_j, self.delta_j)
        if self.alpha_j == -9999 or self.delta_j == -9999:
            self.position = 0
        else:
            self.position = 1 #Understood position is given
        
        print(3)
        
        self.ximg, self.yimg = self._get_ximg_yimg(pdb)
        
        print(self.ximg, self.yimg)
        self.bxcntr, self.bycntr = self._get_bkg_center(pdb)

        try:
            self.mag_zero = float(pdb["mzero"])
        except:
            self.mag_zero = self.mag_zero

        try:
            self.UserGivenPsf = pdb["star"]
            print('PSF is assigned individually to galaxies')
        except:
            self.UserGivenPsf = 'None'

        try:
            self.UserGivenSky = float(pdb['sky'])
            print('Sky is assigned individually to galaxies')
        except:
            self.UserGivenSky = None
        print(4)

        
        
        # Crashhandling starts
        if self.crashhandler:# & starthandle:
            print("CrashHandler is invoked")

            ut.CrashHandlerToRemove(gal_id, self.fstring, self.OUTDIR)

            self._crash_handler(pdb)

        # Determine Search Radius 
        self._search_radius()
        
        print(5)
        
        # first count the number of "potential" targets in the search radius
        
        good_object = self._potential_target()
        print(6)
        #print(2, imgsize, whtsize)
        ##XXX, xcntr or ximg?
        xcntr = good_object[1]
        ycntr = good_object[2]
        # now fit best object            
        #XXX
                
        if 1:#try:
            
            print('Shallow', self.sex_cata)
            self._get_shallow_sky(good_object[10])
           

            #XXX Need to check whether a few line should be needed  
            if self.position == 0:
                alpha_s = 9999
                delta_s = 9999

             #Assume that when alpha_j is -9999 then delta_j will also be -9999
            if self.position & (self.alpha_j == -9999):
                self.alpha_j = alpha_s
                self.delta_j = delta_s

            print('alpha', self.alpha_j)
            print(0, self.gimg)

            self._flags_initialization()
            
            
            # Calculating the cutout size (half size).
            # SizeX, SizeY return from make_cutout are full sizes
            self.SizeX, self.SizeY = super()._find_cutout_size()
            print('Calculated half sizes {}, {}'.format(self.SizeX, self.SizeY))

            # Finding psf and the distance between psf and image
            if self.decompose:
                psffile, distance = super()._handle_psf(self.UserGivenPsf, 
                                                   self.alpha_j, self.delta_j)
            else:
                psffile, distance = 'None', 9999

            print('psffile: {}, distance: {}'.format(psffile, distance))

            
            run = 1 #run =1 if pipeline runs sucessfuly

            ut.WriteError('\n\n###########   {} ###########\n'.format(self.gal_id))

            #print(3, imgsize, whtsize)
            # For the new run
            if self.repeat:
                cut_xcntr = self.TX / 2.0
                cut_ycntr = self.TY / 2.0
                if os.path.exists(mimg):
                    pass
                else:
                    MF = MaskFunc(mimg, 
                                  cut_xcntr, cut_ycntr, 
                                  self.SizeX, self.SizeY, 
                                  good_object)
                    MF.gmask(self.threshold, self.thresh_area,
                             0, NoMask=False)

            else:
                if self.galcut:
                    #Handling image cutout names
                    #XXX
                    self._galcut_images()
                    cut_xcntr, cut_ycntr = self._galcut_cutout()
                else:
                    cut_xcntr, cut_ycntr = self._not_galcut_cutout()

                print(2, self.gimg)

                print('Center of cutimage and exceed size ', \
                      cut_xcntr, cut_ycntr, self.ExceedSize)
                print('Full Sizes ', self.SizeX, self.SizeY)

                if self.galcut & self.ReSize == 0:
                    pass
                if self.ExceedSize:
                    self.flag = SetFlag(self.flag, GetFlag('EXCEED_SIZE'))

                # Runs sextractor to find the segmentation map
                ##RunSegSex(os.path.join(self.DATADIR, self.gimg))

                # This creates ellipse mask and used for 
                # ellipse fitting and casgm
                seg_cata = os.path.join(self.DATADIR,
                                        '{}.seg'.format(self.sex_cata))

                
                self.gimg = os.path.join(self.DATADIR, self.gimg)
                self.wimg = os.path.join(self.DATADIR, self.wimg)

                seg_file = os.path.join(self.OUTDIR, 'check.fits')

                #if os.path.exists(seg_file):
                #    os.remove(seg_file)
                if os.path.exists(seg_cata):
                    os.remove(seg_cata)

                #sys.exit()
                PS = PySex(self.SEX_PATH)
                PS.RunSex(self.sex_params, self.gimg, self.wimg, 
                          seg_cata, self.SEx_GAIN, check_fits=seg_file, 
                          sconfig='seg')
                #sys.exit()
                EM = ElliMaskFunc(cut_xcntr, cut_ycntr,
                                  center_limit=5., seg_limit=1e-5)
                EM.emask(seg_file, seg_cata, self.fstring)

                print(11)
                #sys.exit()
                # Fitting ellipse task or the manual 1d finder
                mimg = 'M_{}.fits'.format(self.fstring)
                print(3, self.gimg)
                if self.decompose:
                    print(4, self.gimg, self.SexSky)
                    ut.FindEllipse(xcntr, ycntr, self.SexHalfRad,
                                  self.SexPosAng, self.axis_rat, self.SexSky, 
                                  self.fstring)

                    MaskFunc(self.gimg, 
                             cut_xcntr, cut_ycntr, 
                             self.SizeX, self.SizeY, 
                             good_object)

                    galfit_conf = 'G_{}.in'.format(self.fstring)
                    oimg = 'O_{}.fits'.format(self.fstring)

                    print(4, self.gimg)

                    print(good_object.shape)
                    CF = GalfitConfigFunc(self.DATADIR,
                                          self.gimg, self.wimg,  
                                          cut_xcntr, cut_ycntr, 
                                          self.SizeX, self.SizeY,
                                          self.components, self.fitting,
                                          psffile,
                                          good_object, 
                                          self.sex_cata,
                                          self.SexSky, self.mag_zero, self.flag)
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
                                        self.gimg, self.wimg,
                                        seg_file, 
                                        cut_xcntr, cut_ycntr,
                                        seg_cata,
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

                print('SexySky, self.SexSky', SexySky, self.SexSky)
                #sys.exit()
                if SkyMin != 9999:
                    SkyMin = SkyYet * 1.0
                    self.SkySig = SkySig * 1.0
                else:
                    SkyMin = self.SexSky * 1.0 
                    self.SkySig = np.sqrt(np.abs(self.SexSky))
                print('Sky Sigma {} '.format(self.SkySig))
                print('SkyMin {} SexSky {}'.format(SkyMin, self.SexSky))
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
                cas_values = ut.HandleCasgm(self.gimg, 
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
                        ConfigIter(self.gimg, self.wimg, 
                                   cut_xcntr, cut_ycntr, 
                                   SizeX, SizeY, 
                                   good_object, 
                                   psffile, z)
                    elif self.galfit:
                        #cmd = '{} {}'.format(self.GALFIT_PATH, galfit_conf)
                        #os.system(cmd)
                        subprocess.call(['{}'.format(self.GALFIT_PATH), '{}'.format(galfit_conf)])

                        f_fit = open('fit2.log','a')
                        if os.path.exists('fit.log'):
                            for line in open('fit.log','r'):
                                f_fit.writelines([str(line)])
                        f_fit.close()
                        # FIX
                        # Mainly the FindEllipse problem
                        ut.OImgFindEllipse(oimg, 
                                           cut_xcntr, cut_ycntr,
                                           self.SexHalfRad,
                                           self.SexPosAng, self.axis_rat, 
                                           self.SexSky, self.fstring, 
                                           self.repeat,
                                           output=False)
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
                                    self.SexSky, SkySig)
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
                            self.gimg, 
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
    #XXX
    #f_cat.close()
    #f_failed.close()






    def find_and_fit(self):
        print('find_and_fit 1')
        #fc = open(self.obj_cata, 'w')
        #if self.redshift == 9999:
        #    fc.writelines(['gal_id ra dec mag\n'])
        #else:
        #    fc.writelines(['gal_id ra dec mag z\n'])

        values = np.genfromtxt(self.sex_cata)
        
        con = (values[17] > self.UMag) & (values[17] < self.LMag) & (values[16] < self.stargal_prob)
        values = values[con]
        values = values[:, [0, 3, 4, 17]]
        values[3] = values[3] / 15.0

        if self.redshift != 9999:
            values[:, -1] = self.redshift
            np.savetxt(self.obj_cata, values, fmt='%i %.8f %.8f %.2f %.3f',
                    header='gal_id ra dec mag z', comments='')
        else:
            np.savetxt(self.obj_cata, values, fmt='%i %.8f %.8f %.2f',
                    header='gal_id ra dec mag', comments='')


        searchrad = '0.05arc'

        print('find_and_fit 2')

    def pymorph(self):


        # now change dir to the 
        thisdir = os.getcwd()
        print("Current directory is ", thisdir)
        print("Output directory is ", self.OUTDIR)

        os.chdir(self.OUTDIR)



        if 1:
            print(1)
            
            if self.repeat == False & self.galcut == False:
                print(2)
                fimg = fitsio.FITS(self.imagefile)
                self.imagedata = fimg[0].read()
                self.header0 = fimg[0].read_header()
                fimg.close()
                self.TX = self.imagedata.shape[0]
                self.TY = self.imagedata.shape[1]
                
                gheader = ut.CheckHeader(self.header0)
                self.EXPTIME = gheader[0]
                self.RDNOISE = gheader[1],
                self.GAIN = gheader[2]
                self.SEx_GAIN = gheader[3]
                self.NCOMBINE = gheader[4]

                print("Using large image. imagefile >>> {}".format(self.imagefile))
            else:
                print(4)
                self.SEx_GAIN = 1.

            if self.repeat == False & self.galcut == False & os.path.exists(self.whtfile):
                self.weightexists = True
                
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
                print('No weight image found\n')


            else:
                self.weightexists = False
            
        #except IOError as e:
        #    print(self.imagefile, "I/O error ({}): {}".format(e.args[0], e.args[1]))
        #    os._exit(0)
        
        print(self.sex_cata)
        if os.path.exists(self.sex_cata):
            pass
            print(5)
        elif self.galcut == False:
            print('The SExtractor catalogue for your frame is NOT found. ' \
                  'One is being made using the default values. It is always '\
                  'recommended to make SExtractor catalogue by YOURSELF as '\
                  'the pipeline keeps the sky value at the SExtractor value '\
                  'during the decomposition.')
            print(6)
            RunSex(self.sex_params, self.SEX_PATH, 
                   self.imagefile, self.whtfile,
                   self.sex_cata, self.SEx_GAIN,
                   seg_file=None, sconfig='default')

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
        print(self.psfselect)
        print(self.obj_cata)
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
                os.system('mv restart.cat CRASH.CAT')
                self.obj_cata = 'CRASH.CAT' 
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
            self.find_and_fit()
            print(self.sex_cata)
            self.main()
            if self.crashhandler:
                self.starthandle = 1
                os.system('mv restart.cat CRASH.CAT')
                self.obj_cata = 'CRASH.CAT' 
                self.main()

        os.chdir(thisdir)





if __name__ == '__main__':
    p = PyMorph()
    p.pymorph()
    #p.main()
