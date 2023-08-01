
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
from .pymorphutils import Get_R, HMSToDeg, DMSToDeg, check_header, write_error, CrashHandlerToRemove, FindEllipse, HandleCasgm, OImgFindEllipse
from .flagfunc import GetFlag, Get_FitFlag, isset, SetFlag

from .ellimaskfunc_easy import ElliMaskFunc
#from ellimaskfunc import ElliMaskFunc

from .maskfunc_easy import MaskFunc
#from maskfunc import MaskFunc

from .configfunc import GalfitConfigFunc
#from configtwostep import ConfigIter
from .yetbackfunc import FindYetSky
from .plotfunc import PlotFunc
from .runsexfunc import PySex
from .writehtmlfunc import WriteHtmlCSV
from .psffunc import update_psf_ra_dec, getpsf
from .mask_or_fit import GetSExObj
       
        


class ReturnClass(object):
    
    def __init__(self):
        pass

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
                       ((ra - ra_p) * np.cos(dec * Get_R()))**2.0)
        except:
            distance = 9999
            print('Distance between psf and image is 9999')

        return distance

    def _handle_psf(self, UserGivenPsf, ra, dec):

        """Determine the psf used for fitting"""

        #print('psflist', self.psflist)
        if self.repeat:
            if np.abs(ra) == 9999 or np.abs(dec) == 9999:
                distance = 9999
                psffile = c.pfile
            else:
                psffile = c.pfile
                update_psf_ra_dec(c.pfile)
                distance_psf_gal = super()._distance_psf_obj(c.pfile, ra, dec)

        else:
            if UserGivenPsf != 'None':
                psffile = UserGivenPsf
                #print('UserGivenPsf', UserGivenPsf)
                update_psf_ra_dec(self.DATADIR, psffile)
                #print('U2')
                #sys.exit()
                psffile, distance_psf_gal = getpsf(self.DATADIR,
                                              UserGivenPsf, self.which_psf,
                                              ra, dec)
                distance_psf_gal = distance_psf_gal * 60. * 60.
                #distance_psf_gal = super()._distance_psf_obj(psffile, ra, dec)
            elif np.abs(ra) == 9999 or np.abs(dec) == 9999:
                psffile = self.psflist[self.psfcounter]
                update_psf_ra_dec(self.DATADIR, psffile)
                distance_psf_gal = 9999
                self.psfcounter += 1
            else:
                psffile, distance_psf_gal = getpsf(self.DATADIR,
                                              self.psflist, self.which_psf,
                                              ra, dec)
                distance_psf_gal = distance_psf_gal * 60. * 60.
        #print(f'psf {psffile} dis {distance_psf_gal}')
        #print('U2')
        return psffile, round(distance_psf_gal, 1)

    
    def _find_cutout_size(self): 
        """Return the size of the cutout. half_sizeX, y_half_sizeY are the size of 
          cutout if c.galcut is True"""

        #print(self.SexPosAng)
        self.SexPosAng = self.SexPosAng * np.pi / 180. #pa in radian
        half1 = self.SexHalfRad * self.FracRad * np.abs(np.cos(self.SexPosAng))
        half2 = self.SexAxisRatio * self.SexHalfRad * self.FracRad * np.abs(np.sin(self.SexPosAng))
        half_size = half1 + half2
        #y_half_size = self.SexHalfRad * self.FracRad * np.abs(np.sin(self.SexPosAng)) + \
        #        self.SexAxisRatio * self.SexHalfRad * self.FracRad * np.abs(np.cos(self.SexPosAng))
        half_size = int(half_size)

        #if self.Square:
        #    half_size = max(half_size, y_half_size)
        #    y_half_size = half_size
        if self.galcut:
            if self.ReSize:
                if self.VarSize:
                    pass
                else:
                    half_size = self.FixSize
            else:
                # FIX
                pass
                # END
        else:
            if self.VarSize:
                pass
            else:
                half_size = self.FixSize
        #print('half_size', half_size)
        return half_size



    def _make_img_square(self, half_size):

        #print(9, self.NXPTS, self.sex_xcntr, self.sex_ycntr, self.half_size)
        ExceedSize = 0
        #All are floor to make the size even number
        xmin = np.floor(self.sex_xcntr - half_size)
        ymin = np.floor(self.sex_ycntr - half_size)
        xmax = np.floor(self.sex_xcntr + half_size)
        ymax = np.floor(self.sex_ycntr + half_size)
        if xmin < 0:
            xmin = 0
            xcntr_img = self.sex_xcntr
            ExceedSize = 1
        else:
            xcntr_img = half_size + np.modf(self.sex_xcntr)[0]
        #print(ymin, ExceedSize)
        if ymin < 0:
            ycntr_img = self.sex_ycntr
            ymin = 0
            ExceedSize = 1
        else:
            ycntr_img = half_size + np.modf(self.sex_ycntr)[0]
        #print(ExceedSize)
        if xmax > self.NXPTS - 1:
            xmax = self.NXPTS
            ExceedSize = 1
        if ymax > self.NYPTS - 1:
            ymax = self.NYPTS
            ExceedSize = 1
        
        if ExceedSize:
            half_size = min(half_size, half_size)
            #half_size = self.half_size
        else:
            half_size = max(half_size, half_size)
        return half_size

    def _make_cutout(self, half_size):
        """

        Make cutout image. The sex_xcntr, sex_ycntr are like iraf.
        half_size, y_half_size are half size

        """
        #print(9, self.sex_xcntr, self.sex_ycntr, half_size)
        ExceedSize = 0
        #All are floor to make the size even number
        xmin = np.floor(self.sex_xcntr - half_size)
        xmax = np.floor(self.sex_xcntr + half_size)
        ymin = np.floor(self.sex_ycntr - half_size)
        ymax = np.floor(self.sex_ycntr + half_size)
        #print('self.NXPTS, self.NYPTS', self.NXPTS, self.NYPTS)
        #print(10, xmin, xmax, self.NXPTS - xmax, ymin, ymax, self.NYPTS - ymax, half_size * 2 - 1)
        #print(11, xmin < 0, xmax > (self.NXPTS - 2), ymin < 0, ymax > (self.NYPTS - 2))
        if xmin < 0 or xmax > (self.NXPTS - 2) or ymin < 0 or ymax > (self.NYPTS - 2):
            max_diff = np.array([xmin, self.NXPTS - 2 - xmax, ymin, self.NYPTS - 2 - ymax])
            con = (max_diff < 0)
            max_diff = abs(max_diff[con]).max() + 1
            half_size = half_size - max_diff #np.max([abs(xmin), abs(ymin), xmax, ymax])
            xcntr_img = half_size + np.modf(self.sex_xcntr)[0]
            ycntr_img = half_size + np.modf(self.sex_ycntr)[0]
            ExceedSize = 1

            #print('inside', ymin, ymax, xmin, xmax)
            xmin = int(self.sex_xcntr - half_size)
            xmax = int(self.sex_xcntr + half_size)
            ymin = int(self.sex_ycntr - half_size)
            ymax = int(self.sex_ycntr + half_size)
            data = self.imagedata[ymin:ymax, xmin:xmax]
        else:
            xcntr_img = half_size + np.modf(self.sex_xcntr)[0]
            ycntr_img = half_size + np.modf(self.sex_ycntr)[0]
            #print('outside', ymin, ymax, xmin, xmax)
            xmin = int(xmin)
            xmax = int(xmax)
            ymin = int(ymin)
            ymax = int(ymax)
            # FIX check c.imagedata or c.ggimage
            #print(10, xmin, xmax, ymin, ymax)
            data = self.imagedata[ymin:ymax, xmin:xmax]
            #print(100, self.imagedata.shape)
            # END
        self.img_NPTS = data.shape[0]

        print('2 xcntr_img ', xcntr_img, ycntr_img)
        #print('half_size', half_size)
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

        #print(13, self.gimg)
        print(13, self.DATADIR, self.cutimage_file)
        #print(14, data)
        print(14)
        #fits = fitsio.FITS(os.path.join(self.DATADIR, self.gimg), 'rw')
        fits = fitsio.FITS(os.path.join(self.DATADIR, self.cutimage_file), 'rw')
        fits.write(data, header=hdict)
        fits.close()
        print(15)
        #sys.exit()

        # FIX
        #Making weight image cut
        if self.weightexists:
            data_wht = self.weightdata[ymin:ymax,xmin:xmax].copy()

            fits = fitsio.FITS(os.path.join(self.DATADIR, self.wimg_file), 'rw')
            fits.write(data_wht)
            fits.close()
        else:
            print('Cannot creat weight image. If you supply weight image \
                   please check whether it exists or report a bug')
        #END
        print('1 xcntr_img', xcntr_img)
        return [xcntr_img, ycntr_img, ExceedSize, half_size]




class Pipeline(ReturnClass):
    
    #__version__ = 3.0

    def __init__(self, config_file='config.ini'):

        pass
        #self.SP = SexParams(self.sex_config, self.mag_zero)


    def _target_initialize(self, pdb):

        if "gal_id" in pdb.keys():
            self.gal_id = int(pdb["gal_id"])
        elif "gimg" in pdb.keys():
            self.gal_id = int(pdb["gimg"].split('.')[0]) #gal_id will be 
                                                   #filename without .fits
        else:
            print("No image or gal_id found in the object catalogue." \
                  "Exiting")
            os._exit(0)

        self.fstring = '{}_{}'.format(self.rootname, self.gal_id)
        print('fstring {}'.format(self.fstring))

        #If we don't set galcut then the below ones are the cutimage_file and
        #wimg_file. Otherwise, this will change in the next condition, that is
        #galcut=True 
        self.cutimage_file = 'I{}.fits'.format(self.fstring)
        self.wimg_file = 'W{}.fits'.format(self.fstring)

        if self.galcut:
            if "gimg" in pdb.keys():
                self.cutimage_file = pdb["gimg"]    #Galaxy cutout
            elif os.path.exists(os.path.join(self.DATADIR, 
                             'I{}.fits'.format(self.fstring))):
                self.cutimage_file = 'I{}.fits'.format(self.fstring)
            elif os.path.exists(os.path.join(self.DATADIR, 
                               '{}.fits'.format(self.gal_id))): 
                self.cutimage_file = '{}.fits'.format(self.gal_id)
            else:
                print("No gimg key is given. No possible gimg found")
                self.cutimage_file = 'None'
                #sys.exit()

            if "wimg" in pdb.keys():
                self.wimg_file = pdb["wimg"]   #Weight cut
            elif os.path.exists(os.path.join(self.DATADIR,
                             'W{}.fits'.format(self.fstring))):
                self.wimg_file = 'W{}.fits'.format(self.fstring)
            else:
                self.wimg_file = 'None'
                print('Search for weight image (wimg) in galfit config')


        if isinstance(pdb["ra"], float): 
            self.alpha_j = pdb["ra"]     
        elif isinstance(pdb["ra"], str):
            h, m, s = pdb["ra"].split(':')
            self.alpha_j = HMSToDeg(int(h), int(m), float(s)) 
        else:
            print("No ra is given")
            self.alpha_j = -9999

        if isinstance(pdb["dec"], float):
            self.delta_j = pdb["dec"]
        elif isinstance(pdb["dec"], str):
            d, m, s = delta.split(':')
            self.delta_j = DMSToDeg(int(d), int(m), float(s))
        else:
            print("No dec is given")
            self.delta_j = -9999

        #print(self.alpha_j, self.delta_j)
        if "z" in pdb.keys():
            self.z = float(pdb["z"])
        else:
            print("No z is given")
            self.z = 9999
        
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
            
        #print('self.SeaPix, self.SeaDeg', self.SeaPix, self.SeaDeg)

            
    def _get_shallow_sky(self, target_sex_xcntr, target_sex_ycntr, target_sex_sky):
        
        # The shallow sky from the shallow run. If no shallow
        # sky it uses the deep sky.

        file_shallow = '{}.Shallow'.format(self.sex_cata)
        if os.path.exists(file_shallow):
            v_shallow = np.genfromtxt(file_shallow)
            try:
                distance_diff = (v_shallow[1] - target_sex_xcntr)**2 + (v_shallow[2] - target_sex_ycntr)**2
                self.SexSky = v_shallow[:, 10][distance_diff.argmin()]
            except:
                pass
        else:
            self.SexSky = target_sex_sky
            
        #print('SexSky', self.SexSky)

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

        #print('Image is >>> {}'.format(self.gimg))
        print('Image is >>> {}'.format(self.cutimage_file))

        #gfits = fitsio.FITS(os.path.join(self.DATADIR, self.gimg), 'r')
        gfits = fitsio.FITS(os.path.join(self.DATADIR, self.cutimage_file), 'r')
        self.imagedata = gfits[0].read()
        self.header0 = gfits[0].read_header()
        gfits.close()

        gheader = check_header(self.header0) #Will set up global header parameters
        self.EXPTIME = gheader[0]
        self.RDNOISE = gheader[1], 
        self.GAIN = gheader[2]
        self.SEx_GAIN = gheader[3]
        self.NCOMBINE = gheader[4]

        if os.path.exists(os.path.join(self.DATADIR, self.wimg_file)):
            wfits = fitsio.FITS(os.path.join(self.DATADIR, self.wimg_file))
            self.weightdata = wfits[0].read()
            wfits.close()
            self.weightexists = True
        else:
            self.wimg_file = None
            self.weightexists = False

        print('Using cutouts')

        #print(1, self.gimg)
        #print(1, self.cutimage_file)
        self.NXPTS = self.imagedata.shape[1]
        self.NYPTS = self.imagedata.shape[0]

        if self.ReSize:
            #self.gimg = 'I{}.fits'.format(self.fstring)
            self.cutimage_file = 'I{}.fits'.format(self.fstring)
            self.wimg_file = 'W{}.fits'.format(self.fstring)

        #The sextractor runs on the cutout before resizing to estimate 
        #shallow sky
    ##if self.galcut == True:   
        #Given galaxy cutouts if the user provides sextractor catalogue
        #then it will not run SExtractor else do!
        PS = PySex(c.SEX_PATH)
        if not os.path.exists(self.sex_cata):
            if 1:

                PS.RunSex(self.sex_config, 
                          os.path.join(self.DATADIR, self.cutimage_file),
                          os.path.join(self.DATADIR, self.wimg_file),
                          os.path.join(self.DATADIR, self.sex_cata),
                          self.SEx_GAIN, check_fits=None, 
                          sconfig='default')

                PS.RunSex(self.sex_config,
                          os.path.join(self.DATADIR, self.cutimage_file),
                          os.path.join(self.DATADIR, self.wimg_file),
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
    
    def _galcut_cutout(self, half_size):
        
        if self.ReSize: 
            #Can take a function
            #Can check whether a decorator can be used
            #if os.path.exists(self.gimg):
            if os.path.exists(self.cutimage_file):
                #write_error('The file {} exists\n'.format(self.gimg))
                write_error('The file {} exists\n'.format(self.cutimage_file))
                run = 0
                #break #Breaking the sextractor loop
            try:
                cntr_half = self._make_cutout(half_size)
            except  Exception as e:
                print(type(e))     # the exception instance
                print(e.args)      # arguments stored in\
                                     # .args
                print(e)           # __str__ allows args\
                                     # to printed directly
                write_error('Cutout exists!')
                print(traceback.print_exc())
                #break
        else:
            xcntr_img = sex_xcntr
            ycntr_img = sex_ycntr 
            ExceedSize = 0
            cntr_half = [xcntr_img, ycntr_img, ExceedSize, half_size]

        cntr_half.append(self.cutimage_file)
        cntr_half.append(self.wimg_file) 
        return cntr_half
        ##return run
                            
            
    def _not_galcut_cutout(self, half_size):
        
        #print(11, self.gimg)
        #print('self.cutimage_file', self.cutimage_file)
        #self.gimg = 'I{}.fits'.format(self.fstring)
        cutimage_file = 'I{}.fits'.format(self.fstring)
        wimg_file = 'W{}.fits'.format(self.fstring)
        #Can take a function
        #if os.path.exists(self.gimg):
        if os.path.exists(self.cutimage_file):
            #write_error('The file {} exists\n'.format(self.gimg))
            write_error('The file {} exists\n'.format(self.cutimage_file))
            run = 0
            #break #Breaking the sextractor loop
            #print(3)

        # Sizes are total size
        #print(12, self.gimg)
        #print(12, self.cutimage_file)
        #XXX
        if 1:#try: 
            cntr_half = self._make_cutout(half_size)
        #except  Exception as e:
        else:
            print(type(e))     # the exception instance
            print(e.args)      # arguments stored in\
                                 # .args
            print(e)           # __str__ allows args\
                                 # to printed directly
            write_error('Cutout exists!')
            print(traceback.print_exc())
            #break

        cntr_half.append(cutimage_file)
        cntr_half.append(wimg_file)
        print('cntr_half ', cntr_half)
        return cntr_half
 

    def _get_ximg_yimg(self, pdb):
           
        '''
        It will be useful for mode galcut
        '''
        if "ximg" in pdb.keys():
            ximg = float(pdb["ximg"])
            if self.ReSize and self.galcut and self.repeat:
                ximg = self.NXPTS / 2.0
        else:
            print('No ximg is given or either ReSize or galcut or repeat '\
                  'keywords are False. Trying to find from the cutout if '\
                  'no ra dec in the image header')
            if self.galcut == True & self.position == 0:
                ximg = self.NXPTS / 2.0
            else:
                ximg = -9999

        if "yimg" in pdb.keys():
            yimg = float(pdb["yimg"])
            if self.ReSize and self.galcut and self.repeat:
                yimg = self.NYPTS / 2.0
        else:
            print('No yimg is given or either ReSize or galcut or repeat '\
                  'keywords are False. Trying to find from the cutout if '\
                  'no ra dec in the image header')
            if self.galcut == True & self.position == 0:
                yimg = self.NYPTS / 2.0
            else:
                yimg = -9999
                
        return ximg, yimg

    
    def _get_bkg_center(self, pdb):
            
        if "bxcntr" in pdb.keys():
            bxcntr = float(pdb["bxcntr"])
        else:
            bxcntr = -9999

        if "bycntr" in pdb.keys():
            bycntr = float(pdb["bycntr"])
        else:
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
        good_object = None

        subprocess.call(['cp', f'{self.sex_cata}', 
                         f'sex_{self.fstring}.txt'])

        #print(2, self.gimg)
        #print(2, self.cutimage_file)

        values_sex = np.genfromtxt(self.sex_cata)
        sex_id = values_sex[0]
        sex_xcntr  = values_sex[:, 1]
        sex_ycntr  = values_sex[:, 2]
        #print(self.sex_cata)
        if self.position == 0:
            alpha_s = np.full_like(values_sex.shape[0], 9999)
            delta_s = np.full_like(values_sex.shape[0], 9999)
            con = (abs(sex_xcntr - self.ximg) < self.SeaPix) & (abs(sex_ycntr - self.yimg) < self.SeaPix)
            curr_distance = np.sqrt((sex_xcntr - self.ximg)**2 + (sex_ycntr - self.yimg)**2)

        else:
            alpha_s = values_sex[:, 3]
            delta_s = values_sex[:, 4]
            con = (abs(self.alpha_j - alpha_s) < self.SeaDeg) & \
            (abs(delta_s - self.delta_j) < self.SeaDeg)
            curr_distance = np.sqrt((self.alpha_j - alpha_s)**2 + (delta_s - self.delta_j)**2)
         
        
        #print(self.alpha_j, alpha_s)
        #print(self.SeaDeg, self.SeaPix)
        
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
            #print('good_object', good_object.shape)

            print("Target distance: {:.3f}".format(center_distance))
        else:
            # No suitable target found
            print("No Target Found!!!!")
            good_object = np.full(19, 9999)
            good_object[11] = 0
            self.flag = SetFlag(self.flag, GetFlag('NO_TARGET'))                 #except:

        #print(good_object)

        
        target = GetSExObj(values=good_object)

        sex_id = target.sex_num
        print("SExtractor ID >>> {}".format(sex_id))


        self.sex_xcntr, self.sex_ycntr = target.sex_center
        #print('sex_xcntr', self.sex_xcntr)
        self.SexHalfRad = target.radius #Sex halfrad 
        self.SexPosAng = target.pos_ang
        self.SexAxisRatio = target.bbya #axis ration b/a
        self.eg = target.eg
        if self.eg <= 0.05:
            self.eg = 0.05

        self.SexMagAuto = target.mag
        self.SexMagAutoErr = target.mag_err
        
        self.LMag = self.SexMagAuto + self.NMag
        self.UMag = self.SexMagAuto - self.NMag

        # FIX
        self.LRadius = 0.2
        self.URadius = self.SexHalfRad * self.NRadius
            # END
        #print(1)

        return good_object, SexTargets
                
    def main(self, obj_value):
        
        pdb = {}                        #The parameter dictionary
        print(obj_value)
        for k, pname in enumerate(self.pnames):
            print(pname)
            pdb[pname] = obj_value[k]
            k += 1

        #print(1, obj_value)
        # declare the flag
        self.flag = 0
        
        
        #print('pnames', self.pnames)
        #print(11)
        self._target_initialize(pdb)
        #print(22)
        #print(2, self.alpha_j, self.delta_j)
        if self.alpha_j == -9999 or self.delta_j == -9999:
            self.position = 0
        else:
            self.position = 1 #Understood position is given
        
        #print(3)
       
        #XXX Need to check whether we use this anywhere 
        self.ximg, self.yimg = self._get_ximg_yimg(pdb)
        
        #print(self.ximg, self.yimg)
        self.bxcntr, self.bycntr = self._get_bkg_center(pdb)

        if "mzero" in pdb.keys():
            self.mag_zero = float(pdb["mzero"])
        else:
            self.mag_zero = self.mag_zero

        #print('pnames', self.pnames)
        if "star" in pdb.keys():
            print(pdb["star"])
            self.UserGivenPsf = pdb["star"]
            print('PSF is assigned individually to galaxies')
        else:
            self.UserGivenPsf = 'None'
        #print(self.UserGivenPsf)
        if "sky" in pdb.keys():
            self.UserGivenSky = float(pdb['sky'])
            print('Sky is assigned individually to galaxies')
        else:
            self.UserGivenSky = None
        #print(4)

        
        
        # Crashhandling starts
        if self.crashhandler:# & starthandle:
            print("CrashHandler is invoked")

            CrashHandlerToRemove(gal_id, self.fstring, self.OUTDIR)

            self._crash_handler(pdb)

        # Determine Search Radius 
        self._search_radius()
        
        #print(5)
        
                #print(6)
        #print(2, imgsize, whtsize)
        # now fit best object            
        #XXX
                
        if 1:#try:
            
                       

            #XXX Need to check whether a few line should be needed  
            if self.position == 0:
                alpha_s = 9999
                delta_s = 9999

             #Assume that when alpha_j is -9999 then delta_j will also be -9999
            if self.position & (self.alpha_j == -9999):
                self.alpha_j = alpha_s
                self.delta_j = delta_s

            #print('alpha', self.alpha_j)
            #print(0, self.gimg)
            #print(0, self.cutimage_file)

            self._flags_initialization()
            
            
            
            # Finding psf and the distance between psf and image
            if self.decompose:
                psffile, distance_psf_gal = super()._handle_psf(self.UserGivenPsf, 
                                                   self.alpha_j, self.delta_j)
            else:
                psffile, distance_psf_gal = 'None', 9999

            print('psffile: {}, distance: {}'.format(psffile, distance_psf_gal))

            if self.galcut:
                #Handling image cutout names
                #XXX
                self._galcut_images()

            good_object, SexTargets = self._potential_target()
            #print('good_object', good_object)

            print('Shallow', self.sex_cata)
            #good_object[10] is the sky value
            self._get_shallow_sky(self.sex_xcntr, self.sex_ycntr, 
                                  good_object[10])

            # Calculating the cutout size (half size).
            # half_size, y_half_size return from make_cutout are full sizes
            self.half_size = super()._find_cutout_size()
 
            #sys.exit() 
            run = 1 #run =1 if pipeline runs sucessfuly

            write_error('\n\n###########   {} ###########\n'.format(self.gal_id))

            #print(3, imgsize, whtsize)
            # For the new run
            if self.repeat:
                xcntr_img = self.NXPTS / 2.0
                ycntr_img = self.NYPTS / 2.0
                if os.path.exists(mimg):
                    pass
                else:
                    MaskF = MaskFunc(mimg, 
                                     xcntr_img, ycntr_img, 
                                     self.NXPTS, self.NYPTS, 
                                     good_object)
                    MaskF.gmask(self.threshold, self.thresh_area,
                             0, NoMask=False)

            else:
                if self.galcut:
                    #Handling image cutout names
                    #XXX
                    #self._galcut_images()
                    xcntr_img, ycntr_img = self._galcut_cutout(self.half_size)
                else:
                    cntr_half = self._not_galcut_cutout(self.half_size)
                    xcntr_img = cntr_half[0] 
                    ycntr_img = cntr_half[1]
                    ExceedSize = cntr_half[2]
                    self.half_size = cntr_half[3]
                    self.cutimage_file = cntr_half[4]
                    self.wimg_file = cntr_half[5]
                #print(2, self.gimg)
                #print(2, self.cutimage_file)

                print('Center of cutimage and exceed size ', \
                      xcntr_img, ycntr_img, ExceedSize)
                #print('Full Sizes ', self.half_size)

                if self.galcut & self.ReSize == 0:
                    pass
                if ExceedSize:
                    self.flag = SetFlag(self.flag, GetFlag('EXCEED_SIZE'))

                # first count the number of "potential" targets in the search radius
        
              
                if self.Square:
                    self.half_size = self._make_img_square(self.half_size)

                #print(self.half_size)
                #sys.exit()
                print('Calculated half sizes {}'.format(self.half_size))

                # Runs sextractor to find the segmentation map
                ##RunSegSex(os.path.join(self.DATADIR, self.gimg))

                # This creates ellipse mask and used for 
                # ellipse fitting and casgm
                seg_cata = os.path.join(self.DATADIR,
                                        '{}.seg'.format(self.sex_cata))

                
                #self.gimg = os.path.join(self.DATADIR, self.gimg)
                self.cutimage_file = os.path.join(self.DATADIR, 
                                                  self.cutimage_file)
                self.wimg_file = os.path.join(self.DATADIR, self.wimg_file)
                
                #seg_fits is used to generate elliptical mask and galfit mask 
                seg_fits = os.path.join(self.OUTDIR, 'check.fits')

                #if os.path.exists(seg_fits):
                #    os.remove(seg_fits)
                if os.path.exists(seg_cata):
                    os.remove(seg_cata)

                #sys.exit()
                # The following is to generate segmentation image for mask
                PS = PySex(self.SEX_PATH)
                #PS.RunSex(self.sex_config, self.gimg, self.wimg_file, 
                #          seg_cata, self.SEx_GAIN, check_fits=seg_fits, 
                #          sconfig='seg')
                PS.RunSex(self.sex_config, self.cutimage_file, self.wimg_file, 
                          seg_cata, self.SEx_GAIN, check_fits=seg_fits, 
                          sconfig='seg')
                #sys.exit()
                EM = ElliMaskFunc(self.DATADIR, xcntr_img, ycntr_img,
                                  center_limit=5., seg_limit=1e-5)
                EM.emask(seg_fits, seg_cata, self.fstring)

                #print(11)
                #sys.exit()
                # Fitting ellipse task or the manual 1d finder
                mimg = 'M_{}.fits'.format(self.fstring)
                mimg = os.path.join(self.DATADIR, mimg)
                #print(sex_xcntr, xcntr_img, sex_ycntr, ycntr_img)
                #sys.exit()
                #print(3, self.gimg)
                #print(3, self.cutimage_file)
                if self.decompose:
                    #print(4, self.gimg, self.SexSky)
                    #print(4, self.cutimage_file, self.SexSky)
                    #print(5, self.gimg, self.SexSky)
                    #print(5, self.cutimage_file, self.SexSky)
                    FE_gimg = FindEllipse(xcntr_img, ycntr_img, 
                                             self.SexHalfRad, self.SexPosAng, 
                                             self.SexAxisRatio, self.SexSky, 
                                             self.fstring)
                    #FE_gimg.profile(self.gimg, output=False)
                    FE_gimg.profile(self.cutimage_file, output=False)
                    #print('Pipeline repeat', self.repeat)
                    #MaskF = MaskFunc(mimg, xcntr_img, ycntr_img, self.NXPTS, self.NYPTS, good_object)
                    #MaskF.gmask(self.threshold, self.thresh_area, seg_fits, seg_cata,         self.avoidme, NoMask=False)
                    #print(4, self.gimg)
                    #print(4, self.cutimage_file)

                    #print(good_object.shape)
                    #print('MZ', self.mag_zero)
                    #XXX Here self.gimg was instead of self.cutimage_file
                    CF = GalfitConfigFunc(self.DATADIR,
                                          self.cutimage_file, self.wimg_file,  
                                          xcntr_img, ycntr_img,
                                          good_object,
                                          self.half_size,
                                          self.components, self.fitting,
                                          psffile,
                                          self.sex_cata,
                                          self.SexSky, 
                                          self.fstring, 
                                          self.threshold, self.thresh_area, 
                                          self.center_deviated, 
                                          self.center_constrain,
                                          self.avoidme,
                                          self.LMag, self.UMag,
                                          self.LN, self.UN,
                                          self.LRadius, self.URadius,
                                          self.bdbox, self.bbox,
                                          self.dbox, self.devauc,
                                          self.galfitv, 
                                          self.mag_zero, self.flag)
                    CF.write_config()
                    #print('CF done')
                    galfit_conf = CF.config_file
                    self.flag = CF.flag
                    oimg = CF.oimg
                    #continue
                    
                    MaskF = MaskFunc(mimg) 
                    #print('seg_fits', seg_fits)
                    MaskF.gmask(self.threshold, self.thresh_area, 
                                seg_fits, CF.fit_neighbor_cutimage,         
                                self.avoidme, NoMask=False)



            # Estimates sky parameters
            #XXX
            if 1:
                #XXX It was self.gimg was instead of self.cutimage_file
                sky_values = FindYetSky(self.sex_config, self.SEX_PATH, 
                                        self.cutimage_file, self.wimg_file,
                                        seg_fits, 
                                        xcntr_img, ycntr_img,
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
                write_error('YetSky estimation failed\n')
                print(traceback.print_exc())

            # Estimate CASGM  
            if self.cas:
                #XXX self.gimg was instead of self.cutimage_file
                cas_values = HandleCasgm(self.cutimage_file, 
                                            xcntr_img, ycntr_img, 
                                            alpha_j, delta_j, 
                                            z, 
                                            half_size, y_half_size, 
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


            # Decomposition
            if os.access('fit.log', os.F_OK):
                os.remove('fit.log')

            if self.decompose:
                if 1:#try:
                    if self.galfit & self.detail:
                        #XXX self.gimg was instead of self.cutimage_file
                        ConfigIter(self.cutimage_file, self.wimg_file, 
                                   xcntr_img, ycntr_img, 
                                   half_size, y_half_size, 
                                   good_object, 
                                   psffile, z)
                    elif self.galfit:
                        cmd = '{} {}'.format(self.GALFIT_PATH, galfit_conf)
                        #print(cmd)
                        #os.system(cmd)
                        subprocess.call(['{}'.format(self.GALFIT_PATH), '{}'.format(galfit_conf)])

                        f_fit = open('fit2.log','a')
                        if os.path.exists('fit.log'):
                            for line in open('fit.log','r'):
                                f_fit.writelines([str(line)])
                        f_fit.close()
                        write_error('((((( Decomposition Successful )))))\n')
                        
                        # FIX
                        # Should be GALFIT sky instead of self.SexSky 
                        self.flag = OImgFindEllipse(oimg, 
                                           xcntr_img, ycntr_img,
                                           self.SexHalfRad,
                                           self.SexPosAng, self.SexAxisRatio, 
                                           self.SexSky, self.fstring, 
                                           self.repeat, self.flag,
                                           self.components,
                                           output=False)
                        #END
                        #print(self.flag)
                        #print(11)
                #except Exception as e:
                #    print(type(e))     # the exception instance
                #    print(e.args)      # arguments stored in\
                                         # .args
                #    print(e)           # __str__ allows args\
                                         # to printed directly
                #    print("Something bad happened (GALFIT)\n\n")
                #    print(traceback.print_exc())

            
                        
            #sys.exit()
            t1 = time.time()
            #print(self.do_plot)
            if self.do_plot == True:
                if os.access('P_' + self.fstring + '.png', os.F_OK):	
                    os.remove('P_' + self.fstring + '.png')
                if not os.path.exists(oimg):
                    print('No output image exists. GALFIT might have been failed')
                    print('Exiting plotting')
                else:
                    try:
                        PlotF = PlotFunc(oimg, mimg, self.fstring, 
                                            self.SexSky, SkySig, self.mag_zero)
                        PlotF.plot_profile()
                    except: 
                        write_error('Error in plotting \n')
                        if mimg == 'None':
                            write_error('Could not find Mask image for plottong \n')
                        run = 0	
                        self.flag = SetFlag(self.flag, GetFlag('PLOT_FAIL'))
            else:
                write_error('No plot \n')
                run = 0	
                self.flag = SetFlag(self.flag, GetFlag('PLOT_NOT_SELECTED'))

            t2 = time.time()
            #print('PF time >>> ', t2 - t1)
            #FIX
            #fout.writelines(['{} '.format(gal_id)])
            #fout.write(good_object)


            #sys.exit()
            if 1:
                WF = WriteHtmlCSV(self.fstring, 
                                  xcntr_img, ycntr_img, 
                                  self.alpha_j, self.delta_j,
                                  self.SexMagAuto, self.SexMagAutoErr, 
                                  SexTargets, self.SexSky,
                                  self.flag, self.SexHalfRad, self.mag_zero,
                                  C, C_err, A, A_err, S, S_err, G, M,
                                  self.components, self.decompose, 
                                  self.repeat, self.detail, 
                                  self.final_result_file,
                                  self.H0, self.WM, self.WV, self.pixelscale,
                                  self.pymorph_config)

                WF.UMag = self.UMag
                WF.LMag = self.LMag
                WF.LRe = self.LRadius
                WF.URe = self.URadius
                WF.LN = self.LN
                WF.UN = self.UN
                WF.LRd = self.LRadius 
                WF.URd = self.URadius

                #print('WRITE HTML', self.params_to_write)
                # self.params_to_write is defined in pymorph.py
                WF.writeparams(self.params_to_write, distance_psf_gal, 
                               self.z)
                self.flag = WF.flag
           # except Exception as e:
           #     print(type(e))     # the exception instance
           #     print(e.args)      # arguments stored in\
           #     # .args
           #     print(e)        # __str__ allows args\
           #     # to printed directly
           #     print("Something bad happened (Writing)!!!!\n\n")
           #     print(traceback.print_exc())

                #except:
                #    write_error('Error in writing html\n')
                #    run = 0
            #print('WF time >>> ', time.time() - t2)
            #The following removes all the temporary files 
            # after every fit
            if isset(self.flag, GetFlag("GALFIT_FAIL")) or \
               isset(WF.fit_flag, Get_FitFlag("LARGE_CHISQ")) or \
               isset(WF.fit_flag, Get_FitFlag("FAKE_CNTR")):# or \
               #isset(WF.fit_flag, Get_FitFlag("BULGE_AT_LIMIT")) or \
               #isset(WF.fit_flag, Get_FitFlag("DISK_AT_LIMIT")):
                f_failed = open(self.restart_file, 'a')
                failedvalues = obj_value.copy()
                for fv in failedvalues:
                    f_failed.writelines(['{} '.format(fv)])
                f_failed.writelines(['{}\n'.format(self.flag)])
                f_failed.close()

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

        #if values[0].strip().isdigit():
        #    print('Something happend in the pipeline. ' + \
        #          'Check error.log 2')
        #else:
        #    pass

        if self.galcut:
            if os.access(self.sex_cata, os.F_OK):
                os.remove(self.sex_cata)
        #sys.exit()
    #except Exception as e:
    #    print(type(e))     # the exception instance
    #    print(e.args)      # arguments stored in\
        # .args
    #    print(e)           # __str__ allows args\
        # to printed directly
    #    print("something bad happened (general object search)!!!!\n\n")
    #    print(traceback.print_exc())
    #XXX
    #fout.close()
    #f_failed.close()







if __name__ == '__main__':
    p = PyMorph()
    p.pymorph()
    #p.main()
