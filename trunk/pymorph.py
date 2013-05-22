#!/data2/home/ameert/python/bin/python2.5

"""PyMorph [Py MOrphological Parameters' Hunter], is a pipeline to find the Morphological parameters of galaxy. Authors: Vinu Vikram , Yogesh Wadadekar, Ajit K. Kembhavi. 2008 Feb, Alan Meert 2010"""

import os
from os.path import exists
import sys
import time
import csv
import traceback
from optparse import OptionParser, OptParseError
import pyfits
import numpy as np
import re

configdir = '.'

print 'configdir is : ', configdir
sys.path.append(configdir)

import config as c
import pymorphutils as ut
from flagfunc import GetFlag, isset, SetFlag

from ellimaskfunc_easy import ElliMaskFunc
#from ellimaskfunc import ElliMaskFunc

from maskfunc_easy import MaskFunc
#from maskfunc import MaskFunc

from configfunc import ConfigFunc
from config_twostep import ConfigIter
from yetbackfunc import FindYetSky, RunSegSex
from plotfunc import PlotFunc
from runsexfunc import RunSex, SexShallow
from writehtmlfunc import WriteParams
    
#from configiter import *
#from configbarpoint import *
#from configbulgedisk import *



def main():
    imagefile = c.imagefile
    whtfile = c.whtfile
    sex_cata = c.sex_cata
    clus_cata = c.clus_cata
    out_cata = c.out_cata
    c.weightexists = 0

    ReSize = c.size[0]
    try:
        VarSize = c.size[1]
    except:
        print "c.size[1] undefined"
        if ReSize:
            VarSize = 1
        else:
            VarSize = 0
    try:
        Square = c.size[3]
    except:
        print "c.size[3] undefined"
        Square = 1
    try:
        FracRad = c.size[2]  
    except:
        print "c.size[2] undefined"
        FracRad = 20
    if VarSize == 0:
        try:
            FixSize = c.size[4]
        except:
            print "c.size[4] undefined"
            FixSize = 120
    threshold = c.threshold
    thresh_area = c.thresh_area
    try:
        os.system('rm -f TmpElliMask.fits TmpElliMask1.fits')
    except:
        pass

    #Initialize index.html
    if exists('index.html'):
        pass
    else:
        indexfile = open('index.html', 'w')
        indexfile.writelines(['<HTML>\n<BODY>\n'])
        indexfile.writelines(['</BODY></HTML>'])
        indexfile.close()

    #Reading image and weigh files
    if(c.repeat == False and c.galcut == False):
            print "Using large image. c.imagefile >>> ", imagefile
            TX = c.imagedata.shape[1]
            TY = c.imagedata.shape[0]
            if exists(c.datadir + whtfile):
                wht = pyfits.open(c.datadir + whtfile)
                if re.search("rms", whtfile.lower()):
                    c.weightdata = wht[0].data
                    print "whtfile >>> ", whtfile
                    c.weightexists = 1
                elif re.search("weight", whtfile.lower()):
                    c.weightdata = 1 / np.sqrt(wht[0].data)
                    print "whtfile >>> ", whtfile
                    c.weightexists = 1
                else:
                    print 'Weight file is not understood. Please include ' + \
                          'the word weight/rms to the weight file name. ' + \
                          'If it is weight the rms will be found by 1/sqrt(w)' 
                wht.close()
            else:
               print 'No weight image found\n'
    
    #Initializing psf array. ie. creating c.psflist from file 
    ut.PsfArr()
    if c.decompose:
        for psfelement in c.psflist:
            ut.UpdatePsfRaDec(psfelement)
    try:
        ComP = c.components
    except:
        print "c.components undefined. Asuming bulge+disk model"
        ComP = ['bulge', 'disk']
    if len(ComP) == 0:
        print "No model specified. Asuming bulge+disk model"
        ComP = ['bulge', 'disk']
    # Writing csvn header and findeing the output parameters
    ParamToWrite = ut.PyMorphOutputParams(c.dbparams, c.decompose)    
    if exists('result.csv'):
        pass
    else:
        f_res = open("result.csv", "ab")
        csvlist = ['%s_%d'%(ParamToWrite[par_key][0], par_key)
                   for par_key in ParamToWrite.keys()]
        print csvlist
        writer = csv.writer(f_res)
        writer.writerow(csvlist)
        f_res.close()

    f_cat = open(out_cata, 'w')
    f_failed = open('restart.cat', 'w')
    obj_file = open(c.datadir + clus_cata, 'r')  # The file contains the 
                                                 # objects of interest
    pnames = obj_file.readline().split() #The names of the parameters given 
                                         #in the first line in the clus_cata

    # writing a input catalogue (restart.cat) for failed objects
    for FailedParam in pnames:
        f_failed.writelines([str(FailedParam), ' '])
    f_failed.writelines(['flag \n'])

    pdb = {}                        #The parameter dictionary
    c.psfcounter = 0                  #For getting psf in the case of unknown ra
    for line_j in obj_file:
        # declare the flag
        c.Flag = 0
        
        try:
            values = line_j.split()
            k = 0
            for pname in pnames:
                pdb[pname] = values[k]
                k += 1
            try:
                gal_id = pdb["gal_id"]
            except:
                print "no gal_id using gimg"
                try:
                    gal_id = pdb["gimg"].split('.')[0] #gal_id will be 
                                                       #filename without .fits
                except:
                    print "No image or gal_id found in the object catalogue." \
                          "Exiting"
                    os._exit(0)
            c.fstring = str(c.rootname) + '_' + str(gal_id)
            try:
                alpha1 = float(pdb["ra1"])
            except:
                print "No ra1 (hour) is given"
                alpha1 = -9999
            try:
                alpha2 = float(pdb["ra2"])
            except:
                print "No ra2 (minuite) is given. Asuming ra is in deg"
                alpha2 = 0
            try:
                alpha3 = float(pdb["ra3"])
            except:
                print "No ra3 (second) is given"
                alpha3 = 0
            try:
                delta1 = float(pdb["dec1"])
            except:
                print "No dec1 (deg) is given"
                delta1 = -9999
            try:
                delta2 = float(pdb["dec2"])
            except:
                print "No dec2 (min) is given"
                delta2 = 0
            try:
                delta3 = float(pdb["dec3"])
            except:
                print "No dec3 (sec) is given"
                delta3 = 0
            if alpha1 == -9999 or delta1 == -9999:
		RaDecInfo = 0
            else:
                RaDecInfo = 1 #Understood position is given
            try:
                z = float(pdb["z"])
            except:
                print "No z is given"
                z = 9999
            try:
                gimg = pdb["gimg"]    #Galaxy cutout
            except:
                print "No gimg given."
                if exists(c.datadir + 'I' + c.fstring + '.fits'):
                    gimg = 'I' + c.fstring + '.fits'
                elif exists(c.datadir + str(gal_id) + '.fits'): 
                    gimg = str(gal_id) + '.fits'
                else:
                    print "No possible gimg found"
                    gimg = 'None'
            print 'fstring ', c.fstring
            try:
                wimg = pdb["wimg"]   #Weight cut
            except:
                if c.galcut:
                    print "No wimg given"
                    if exists(c.datadir + 'W' + c.fstring + '.fits'):
                        wimg = 'W' + c.fstring + '.fits'
                    else:
                        print "No possible wimg found"
                        wimg = 'None'
            try:
                cfile = pdb["cfile"]  #GALFIT configuration file
            except:
                print "No cfile given"
                if c.repeat == True:
                    cfile = 'G_' + c.fstring + '.in'
                else:
                    print "Repeat is false. No possible cfile (galfit " + \
                          "config file) is using"
                    cfile = 'None'
            #Reading galfit config file, if it exits, to know the gimg etc.
            if exists(cfile):
                gimg, oimg, wimg, c.pfile, mimg, cofile = \
                                                     ut.ReadGalfitConfig(cfile)
                gimg = gimg.split('/')[-1]
                oimg = oimg.split('/')[-1]
                wimg = wimg.split('/')[-1]
                c.pfile = c.pfile.split('/')[-1]
                mimg = mimg.split('/')[-1]
                cofile = cofile.split('/')[-1]
            else:
                cfile = 'None'
            #Handling image cutout names
            if c.galcut == True:
                if ReSize:
                    cutimage = 'I' + c.fstring + '.fits'
                    whtimage = 'W' + c.fstring + '.fits'
                else:
                    cutimage = gimg
                    whtimage = wimg                
            else:
                cutimage = 'I' + c.fstring + '.fits'
                whtimage = 'W' + c.fstring + '.fits'
            if c.galcut == True:
                print 'Image is >>> ', gimg
                ggimg = pyfits.open(c.datadir + gimg)
                c.imagedata = ggimg[0].data
                header0 = ggimg[0].header
                ggimg.close()
                ut.CheckHeader(header0) #Will set up global header parameters
                TX = c.imagedata.shape[1]
                TY = c.imagedata.shape[0]
                if exists(c.datadir + whtimage):
                    gwimg = pyfits.open(c.datadir + wimg)
                    c.weightdata = gwimg[0].data
                    gwimg.close()
                    c.weightexists = 1
                else:
                    whtimage = 'None'
                print 'Using cutouts'
            try:
                ximg = float(pdb["ximg"])
                if ReSize and c.galcut and c.repeat:
                    ximg = TX / 2.0
            except:
                print 'No ximg is given. Trying to find from the cutout if '\
                      'no ra dec in the image header'
                if(c.galcut == True and RaDecInfo == 0):
                    ximg = TX / 2.0
                else:
                    ximg = -9999
            try:
                yimg = float(pdb["yimg"])
                if ReSize and c.galcut and c.repeat:
                    yimg = TY / 2.0
            except:
                print 'No yimg is given. Trying to find from the cutout if '\
                      'no ra dec in the image header'
                if(c.galcut == True and RaDecInfo == 0):
                    yimg = TY / 2.0
                else:
                    yimg = -9999
            try:
                bxcntr = float(pdb["bxcntr"])
            except:
                bxcntr = 9999
            try:
                bycntr = float(pdb["bycntr"])
            except:
                bycntr = 9999
            try:
                c.mag_zero = float(pdb["mzero"])
            except:
                c.mag_zero = c.mag_zero
            try:
                UserGivenPsf = pdb["star"]
                print 'Psf is assigned individually to galaxies'
            except:
                UserGivenPsf = 'None'
            try:
                UserGivenSky = pdf['sky']
                print 'Sky is assigned individually to galaxies'
            except:
                UserGivenSky = -9999

            # Crashhandling starts
            if c.crashhandler and c.starthandle:
                print "CrashHandler!!!"
                CrashHandlerToRemove(gal_id)
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
                        c.fitting[2] = 0
                    else:
                        c.fitting[2] = 1
                    if isset(CrashFlag, GetFlag("FIT_BULGE_CNTR")) and\
                       isset(CrashFlag, GetFlag("FIT_DISK_CNTR")):
                        pass
                    else:
                        c.fitting[0] = 1
                        c.fitting[1] = 1
                if isset(CrashFitFlag, Get_FitFlag("LARGE_CHISQ")):
                    if isset(CrashFlag, GetFlag("FIT_BULGE_CNTR")) and\
                       isset(CrashFlag, GetFlag("FIT_DISK_CNTR")):
                        pass
                    else:
                        c.fitting[0] = 1
                        c.fitting[1] = 1
                if isset(CrashFitFlag, Get_FitFlag("FAKE_CNTR")):
                    c.center_deviated = 1

            #The sextractor runs on the cutout before resizing to estimate 
            #shallow sky
            if(c.galcut == True):   #Given galaxy cutouts
                if exists(sex_cata): #If the user provides sextractor catalogue
                                     #then it will not run SExtractor else do!
                    pass
                else:
                    try:
                        RunSex(c.datadir + gimg, c.datadir + wimg, 'None', \
                               9999, 9999, 0)
                        SexShallow(c.datadir + gimg, c.datadir + wimg, \
                                   'None', 9999, 9999, 0)
                    except Exception, inst:
                        print type(inst)     # the exception instance
                        print inst.args      # arguments stored in\
                                             # .args
                        print inst           # __str__ allows args\
                                             # to printed directly
                        print "something bad happened!!!!\n\n"
                        print traceback.print_exc()
                        
                        #print 'Problem running Sextractor (line no. 342)' 
            if(alpha1 == -9999 or delta1 == -9999):
                alpha_j = -9999
                delta_j = -9999
            else:
                alpha_j = ut.HMSToDeg(alpha1, alpha2, alpha3)
                delta_j = ut.DMSToDeg(delta1, delta2, delta3)
            # Determine Search Radius 
            try:
                SearchRad = c.searchrad
            except:
                if RaDecInfo:
                    SearchRad = '1arc'
                    print 'No search radius found. Setting to 1 arc sec'
                else:
                    SearchRad = '10pix'
                    print 'No search radius found. Setting to 10 pix'
            if SearchRad.endswith('arc'):
                SeaDeg = float(SearchRad[:-3]) / (60.0 * 60.0)
                SeaPix = 10.0
            elif SearchRad.endswith('pix'):
                SeaPix = float(SearchRad[:-3])
                SeaDeg = c.pixelscale * SeaPix  / (60.0 * 60.0)

            # first count the number of "potential" targets in the search radius
            c.SexTargets = 0
            good_objects = []
            bad_objects = []
            good_distance = []
            bad_distance = []
            good_object = ''

            os.system('cp %s sex_%s.txt' %(sex_cata, c.fstring))
            for line_s in open(sex_cata, 'r'):
                try:
                    values = line_s.split()
                    alpha_s = float(values[3])
                    delta_s = float(values[4])
                    if RaDecInfo == 0:
                        alpha_s = 9999
                        delta_s = 9999
                    sex_id = values[0]
                    xcntr  = float(values[1])
                    ycntr  = float(values[2])
                    
                    if(abs(alpha_j - alpha_s) < SeaDeg and \
                       abs(delta_s - delta_j) < SeaDeg or \
                       abs(xcntr - ximg) < SeaPix and \
                       abs(ycntr - yimg) < SeaPix):
                        c.SexTargets +=1
                        if c.SexTargets == 1:
                            good_object = line_s
                except:
                    if values[0].strip().isdigit():
                        print 'Something happend in the pipeline. ' + \
                              'Check error.log'
                    else:
                        pass
            if len(good_object) <1:
                # No suitable target found
                print "NO TARGET FOUND!!!!"
                good_object = ' 9999  9999 9999  9999 9999  9999 9999  9999 9999  9999 9999  0 9999  9999 9999  9999 9999 9999 9999\n'
                c.Flag = SetFlag(c.Flag,GetFlag('NO_TARGET'))  
                

            # now fit best object            
            try:
                values = good_object.split()
                alpha_s = float(values[3])
                delta_s = float(values[4])
                if RaDecInfo == 0:
                    alpha_s = 9999
                    delta_s = 9999
                sex_id = values[0]
                xcntr  = float(values[1])
                ycntr  = float(values[2])
                print "SExtractor ID >>> ", values[0]
                c.SexMagAuto = float(values[17])
                c.SexMagAutoErr = float(values[18])
                c.SexHalfRad = float(values[9]) #Sex halfrad 
                c.SexPosAng = float(values[11])
                c.pos_ang = ut.pa(float(values[11]))
                c.axis_rat = 1.0 / float(values[12]) #axis ration b/a
                c.eg = 1 - c.axis_rat
                if c.eg <= 0.05:
                    c.eg = 0.07
                c.major_axis = float(values[14])
                if c.UMag > -9999.0:
                    c.UMag = c.SexMagAuto - 7.0
                if c.LMag < 9999.0:
                    c.LMag = c.SexMagAuto + 7.0

                # FIX
                if c.URe < 9999.0:
                    c.URe = c.SexHalfRad * 50.0
                if c.URd < 9999.0:
                    c.URd = c.SexHalfRad * 50.0
                # END
                # The shallow sky from the shallow run. If no shallow
                # sky it uses the deep sky
                ShallowSky = 9999
                if exists(str(c.sex_cata) + '.Shallow'):
                    f_sex_shallow = open(str(c.sex_cata) + \
                                    '.Shallow', 'r')
                    for line_shallow in f_sex_shallow:
                        v_shallow = line_shallow.split()
                        try:
                            if str(v_shallow[19]) == str(values[0]):
                                ShallowSky = float(v_shallow[10])
                        except:
                            pass
                    f_sex_shallow.close()
                if ShallowSky == 9999:
                    c.SexSky = float(values[10])
                    c.GalSky = 9999
                else:
                    c.SexSky = ShallowSky
                    c.GalSky = 9999
                if(alpha_j == -9999 or delta_j == -9999):
                    if RaDecInfo: 
                        alpha_j = alpha_s
                        delta_j = delta_s
                ut.WriteError('\n\n###########   ' + str(gal_id) + \
                                  '   ###########\n')
                c.run = 1 #run =1 if pipeline runs sucessfuly
                # Adding initial setup to the flag
                if c.repeat:
                    c.Flag = SetFlag(c.Flag,GetFlag('REPEAT'))
                if c.fitting[0] and 'bulge' in ComP:
                    c.Flag = SetFlag(c.Flag,GetFlag('FIT_BULGE_CNTR'))
                if c.fitting[1] and 'disk' in ComP:
                    c.Flag = SetFlag(c.Flag,GetFlag('FIT_DISK_CNTR'))
                if c.fitting[2]:
                    c.Flag = SetFlag(c.Flag,GetFlag('FIT_SKY'))
                
                # Calculating the cutout size (half size).
                # SizeX, SizeY return from MakeCutOut are full sizes
                SizeX, SizeY = ut.FindCutSize(ReSize, VarSize, \
                               Square, FracRad, c.size[4], TX/2, TY/2) 
                print 'Calculated half sizes ', SizeX, SizeY
                # Finding psf and the distance between psf and image
                if c.decompose:
                    psffile, distance = ut.HandlePsf(cfile, \
                                       UserGivenPsf, alpha_j, delta_j)
                    print 'psffile, distance > ', psffile, distance
                else:
                    psffile, distance = 'None', 9999
                # For the new run
                if c.repeat == False:
                    if c.galcut == False:
                        if exists(cutimage):
                            ut.WriteError('The file ' + cutimage +\
                                          ' exists\n')
                            c.run = 0
                            break #Breaking the sextractor loop
                        # Sizes are total size
                        try: 
                            cut_xcntr, cut_ycntr, SizeX, SizeY, \
                            ExceedSize = \
                            ut.MakeCutOut(xcntr, ycntr, alpha_j, \
                                    delta_j, SizeX, SizeY, \
                                    TX, TY, cutimage, whtimage, ReSize) 
                        except:
                            ut.WriteError('Cutout exists!')
                            break
                    if c.galcut:
                        if ReSize: 
                            if exists(cutimage):
                                ut.WriteError('The file ' + cutimage +\
                                              ' exists\n')
                                c.run = 0
                                break #Breaking the sextractor loop
                            try:
                                cut_xcntr, cut_ycntr, SizeX, SizeY, \
                                ExceedSize = \
                                ut.MakeCutOut(xcntr, ycntr, alpha_j, \
                                              delta_j, SizeX, SizeY, \
                                    TX, TY, cutimage, whtimage, ReSize)
                            except:
                                ut.WriteError('Cutout exists!')
                                break
                        else:
                            cut_xcntr, cut_ycntr, SizeX, SizeY, \
                            ExceedSize = xcntr, ycntr, TX, TY, 0
                    print 'Center of cutimage and exceed size ', \
                          cut_xcntr, cut_ycntr, ExceedSize
                    print 'Full Sizes ', SizeX, SizeY
                    if c.galcut and ReSize == 0:
                        pass
                    elif ExceedSize:
                        c.Flag = SetFlag(c.Flag,GetFlag('EXCEED_SIZE'))
                    # Runs sextractor to find the segmentation map
                    RunSegSex(c.datadir + cutimage)

                    # This creates ellipse mask and used for 
                    # ellipse fitting and casgm
                    ElliMaskFunc(cutimage, cut_xcntr, cut_ycntr, \
                                 SizeX, SizeY, good_object, 1)
                    # Fitting ellipse task or the manual 1d finder
                    if c.decompose:
                        ut.HandleEllipseTask(cutimage, cut_xcntr, \
                                       cut_ycntr, \
                                       SizeX, SizeY, c.SexSky, 0)
                        MaskFunc(cutimage, cut_xcntr, cut_ycntr, \
                                         SizeX, SizeY, good_object)
                        maskimage = 'M_' + c.fstring  + '.fits'
                        config_file = 'G_' + c.fstring + '.in'
                        outimage = 'O_' + c.fstring + '.fits'
                        ConfigFunc(cutimage, whtimage,  cut_xcntr,\
                                   cut_ycntr, SizeX, SizeY, good_object, \
                                   psffile, 'SegCat.cat')
                        #continue
                else:
                    cut_xcntr, cut_ycntr = TX / 2.0, TY / 2.0
                    maskimage = mimg
                    if exists(mimg):
                        pass
                    else:
                        MaskFunc(cutimage, cut_xcntr, cut_ycntr, \
                                 SizeX, SizeY, good_object)
                    config_file = cfile
                    outimage = str(oimg)
                # Estimates sky parameters
                try:
                    SexySky, SkyYet, SkyMed, SkyMin, SkyQua, \
                    SkySig = \
                    FindYetSky(c.datadir + cutimage, \
                              cut_xcntr, cut_ycntr)
                    if SkyMin != 9999:
                        c.SkyMin = SkyYet * 1.0
                        c.SkySig = SkySig * 1.0
                    else:
                        c.SkyMin = c.SexSky * 1.0 
                        c.SkySig = np.sqrt(np.abs(c.SexSky))
                    print 'Sky Sigma >>> ', c.SkySig
                    print 'SkyMin SexSky > ', c.SkyMin, c.SexSky
                except:
                    ut.WriteError('Sky estimation failed\n')
                # Estimate CASGM  
                if(c.cas):
                    C, C_err, A, A_err, S, S_err, G, M = \
                    ut.HandleCasgm(cutimage, cut_xcntr, cut_ycntr, \
                                   alpha_j, delta_j, z,
                                   SizeX, SizeY, good_object, bxcntr, \
                                   bycntr)
                else:
                    C, C_err, A, A_err, S, S_err, G, M = \
                    9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999
                # Removing CASGM temp files
                for ff in ['BMask.fits', 'MRotated.fits', \
                          'MaskedGalaxy.fits', 'Rotated.fits']:
                    if os.access(ff, os.F_OK):
                        os.remove(ff)
                # Decomposition
                if os.access('fit.log', os.F_OK):
                    os.remove('fit.log')
                if c.decompose:
                    try:
                        try:
                            DetailFit = c.detail
                        except:
                            DetailFit = 0
                        if c.galfit and DetailFit:
                            ConfigIter(cutimage, whtimage, cut_xcntr,\
                                       cut_ycntr, SizeX, \
                                       SizeY, good_object, psffile, z)
                        elif c.galfit:
                            cmd = str(c.GALFIT_PATH) + ' ' + \
                                     config_file
                            os.system(cmd)
                            f_fit = open('fit2.log','a')
                            if exists('fit.log'):
                                for line in open('fit.log','r'):
                                    f_fit.writelines([str(line)])
                            f_fit.close()
                            # FIX
                            # Mainly the ellipse fit problem
                            ut.HandleGalfitOutput(cutimage, outimage, \
                                                  cut_xcntr, cut_ycntr, \
                                                  SizeX, SizeY, good_object)
                            #END
                    except Exception, inst:
                        print type(inst)     # the exception instance
                        print inst.args      # arguments stored in\
                                             # .args
                        print inst           # __str__ allows args\
                                             # to printed directly
                        print "something bad happened!!!!\n\n"
                        print traceback.print_exc()

                
                            
                    if(c.run == 1):
                        ut.WriteError('((((( Decomposition '\
                                          'Successful )))))\n')

                try:
                    Goodness = -9999
                    if os.access('P_' + c.fstring + '.png', \
                                 os.F_OK):	
                        os.remove('P_' + c.fstring + '.png')
                        GoodNess = PlotFunc(outimage, \
                                 maskimage, cut_xcntr, cut_ycntr, \
                                            c.SexSky, c.SkySig)
                        Goodness = GoodNess.plot_profile
                except:
                    ut.WriteError('Error in plotting \n')
                    if maskimage == 'None':
                        ut.WriteError('Could not find Mask image\n')
                    c.run = 0	
                    c.Flag = SetFlag(c.Flag,GetFlag('PLOT_FAIL'))
                
                if isset(c.Flag, GetFlag("GALFIT_FAIL")): #or \
                   #isset(c.Flag, GetFlag("LARGE_CHISQ")) or \
                   #isset(c.Flag, GetFlag("FAKE_CNTR")) or \
                   #isset(c.Flag, GetFlag("BULGE_AT_LIMIT")) or \
                   #isset(c.Flag, GetFlag("DISK_AT_LIMIT")):
                    FailedValues = line_j.split()
                    for FailedValue in FailedValues:
                        f_failed.writelines([str(FailedValue), ' '])
                    f_failed.writelines([str(c.Flag), '\n'])
                f_cat.writelines([str(gal_id), ' '])
                f_cat.write(good_object)
                #The following removes all the temporary files 
                # after every fit
                ToClean = 0


                try:
                    WriteParams(ParamToWrite, cutimage, cut_xcntr, cut_ycntr, \
                                distance, alpha_j, \
                                delta_j, z, Goodness, \
                                C, C_err, A, A_err, S, S_err, \
                                G, M, c.EXPTIME)
                except Exception, inst:
                    print type(inst)     # the exception instance
                    print inst.args      # arguments stored in\
                    # .args
                    print inst           # __str__ allows args\
                    # to printed directly
                    print "something bad happened!!!!\n\n"
                    print traceback.print_exc()

                    #except:
                    #    ut.WriteError('Error in writing html\n')
                    #    c.run = 0

                if ToClean:
                    ClaCli = 'E*fits E*txt G_* I*fits *.pl \
                         *.con M_* \
                         O*fits O*txt P*png R*html error.log  \
                         galfit.*  \
                         Tmp* SO* agm_r* \
                         BMask.fits MaskedGalaxy.fits \
                         MRotated.fits   \
                         B.fits GalEllFit.fits AResidual.fits \
                         ellip err BackMask.fits'
                    vClaCli = ClaCli.split()
                    for v1ClaCli in vClaCli:
                        if os.access(v1ClaCli, os.F_OK):
                            os.remove(v1ClaCli)
                for myfile in ['ellip','err','test.tab']:
                    if os.access(myfile,os.F_OK):
                            os.remove(myfile)
            except Exception, inst:
                print type(inst)     # the exception instance
                print inst.args      # arguments stored in\
                                            # .args
                print inst           # __str__ allows args\
                                             # to printed directly
                print "something bad happened!!!!\n\n"
                print traceback.print_exc()

                #NOTE CHANGE THIS ADD ADDITIONAL IF len > 0 HERE

		if values[0].strip().isdigit():
                    print 'Something happend in the pipeline. ' + \
                          'Check error.log 2'
                else:
                    pass
            if(c.galcut == True):
                if os.access(sex_cata, os.F_OK):
                    os.remove(sex_cata)
        except Exception, inst:
            print type(inst)     # the exception instance
            print inst.args      # arguments stored in\
            # .args
            print inst           # __str__ allows args\
            # to printed directly
            print "something bad happened!!!!\n\n"
            print traceback.print_exc()
    f_cat.close()
    f_failed.close()
def selectpsf(ImG, CaT):
    c.psff = []
    im1 = pyfits.open(ImG)
    image = im1[0].data
    im1.close()
    AreaOfObj = c.AreaOfObj
    def FindPsf(AreaOfObj, CaT):
        for line in open(CaT,'r'):
            values = line.split()
            try:
                BaKgR  = float(values[10])
                if float(values[16]) >= c.StarGalProb and \
		   float(values[13]) > AreaOfObj \
                   and float(values[14]) < 50.0:
                    xcntr = float(values[1]) - 1
                    ycntr = float(values[2]) - 1
                    #size of the psf is 8 times the sigma assuming the star has 
                    #Gaussian profile
                    PsfSize = np.floor(float(values[14])) * c.starsize 
                    x1 = int(xcntr) + (PsfSize/2)
                    x2 = int(xcntr) - (PsfSize/2)
                    y1 = int(ycntr) + (PsfSize/2)
                    y2 = int(ycntr) - (PsfSize/2)
                    ra1 = int(float(values[3]) / 15.0)
                    ra2 = int((float(values[3]) / 15.0 - int(float(values[3])\
                          / 15.0))*60.0)
                    ra3 = (((float(values[3]) / 15.0 - int(float(values[3]) / \
                          15.0))*60.0) - ra2) * 60.0
                    dec1 = int(float(values[4]))
                    dec2 = abs(int((float(values[4]) - dec1) * 60.0))
                    dec3 = (abs(((float(values[4]) - dec1) * 60.0)) - dec2) \
                           * 60.0
                    if ra1 < 10:
                        ra11 = '0' + str(ra1)
                    else:
                        ra11 = str(ra1)
                    if ra2 < 10:
                        ra22 = '0' + str(ra2)
                    else:
                        ra22 = str(ra2)
                    if ra3 < 10:
                        ra33 = '0' + (str(np.round(ra3, 1))[:3]).split('.')[0]\
                                      + \
                                     (str(np.round(ra3, 1))[:3]).split('.')[1] 
                    else:
                        ra33 = (str(np.round(ra3, 1))[:4]).split('.')[0] + \
                               (str(np.round(ra3, 1))[:4]).split('.')[1]
                    if abs(dec1) < 10:
                        if dec1 < 0.0:
                            dec11 = '-0' + str(abs(dec1))
                        else:
                            dec11 = '+0' + str(abs(dec1))
                    else:
                        if dec1 < 0.0:
                            dec11 = str(dec1)
                        else:
                            dec11 = '+' + str(dec1)
                    if dec2 < 10:
                        dec22 = '0' + str(dec2)
                    else:
                        dec22 = str(dec2)
                    if dec3 < 10:
                        dec33 = '0' + (str(np.round(dec3, 1))[:3]).split('.')[0]\
                                 + (str(np.round(dec3, 1))[:3]).split('.')[1]
                    else:
                        dec33 = (str(np.round(dec3, 1))[:4]).split('.')[0] +\
                                (str(np.round(dec3, 1))[:4]).split('.')[1]
                    psffile = 'psf_' + str(ra11) + str(ra22) + str(ra33) + str(dec11) +str(dec22) + str(dec33) + '.fits'
                    if psffile in c.psff:
                        pass
                    else:
                        c.psff.append(psffile)
                        psf = image[y2:y1, x2:x1]
                        psf = psf - BaKgR
                        if os.access(psffile, os.F_OK):
                            os.remove(psffile)
                        hdu = pyfits.PrimaryHDU(psf.astype(np.float32))
                        hdu.header.update('XCNTR', int(xcntr))
                        hdu.header.update('YCNTR', int(ycntr))
                        hdu.writeto(psffile)
            except:
                pass
    while len(c.psff) < 5 and AreaOfObj > 10:
        FindPsf(AreaOfObj, CaT)
        AreaOfObj -= 5
    if len(c.psff) == 0:
        manualpsf = raw_input("Unfortunately, NO psf is found. Please enter "\
                              "the psf name >>> ")
        c.psff.append(manualpsf)
    if c.Interactive:
        PsfList = []
        TmPLST = []
        for element in c.psff:
            TmPLST.append(element)
        print 'Checking Started. You can just visually check the psf. You can'\
               ' do the thorough checking later'
        for element in c.psff:
            os.system('xpaset -p ds9 frame clear')
            if exists(element):
                os.system('cat ' + str(element) + ' | xpaset ds9 fits')
                os.system('xpaset -p ds9 scale mode zscale')
                os.system('xpaset -p ds9 zoom to fit')
            else:
                print 'The psf you have given does NOT exists!!!'
                pass
            write = raw_input("Do you need this psf? ('y' if yes, 'c'"\
                              " to cancel psf checking) " )
            if write == 'y':
                PsfList.append(element)
                TmPLST.remove(element)
#                c.psff.remove(element)
            elif write == 'c':
                for element1 in TmPLST:
                    try:
                        os.remove(element1)
                    except:
                        pass
                break
            else:
                try:
                    os.remove(element)
                except:
                    pass
        print '\nFinal Checking Started. If you are using psfselect = 2, be '\
              'carefull!!! This is your last chance for selecting psf. '\
              'Select ONLY GOOD psf. Press "y" to accept the previous psf. '\
             'ENTER to go to the next one. This will continue until you press '\
             '"1" when it asked Finished? ALL THE BEST! \n'
        UrPsfChk = raw_input("Do you want to use your own psf? Enter 'y' or " \
                             "'n' >>> ")
        if UrPsfChk == 'y':
            UrPsf = raw_input("Enter your psf >>> ")
            fff = open('psflist.list', 'ab')
            fff.writelines([str(UrPsf), '\n'])
            fff.close()
            if(c.galcut):
                c.ValueS.append(UrPsf)
                fwithpsf = open('CatWithPsf.cat', 'ab')
                for v in c.ValueS:
                    fwithpsf.writelines([str(v), ' '])
                fwithpsf.writelines(['\n'])
                fwithpsf.close()
            finish = 1
            for element in PsfList:
                if os.access(element, os.F_OK):
                    os.remove(element)
        else:
            finish = 0
            TmpPsfList = []
        UpdateCounter = 1
        while finish == 0:
            for element in PsfList:
                if exists(element):
                    os.system('xpaset -p ds9 frame clear')
                    os.system('cat ' + str(element) + ' | xpaset ds9 fits')
                    os.system('xpaset -p ds9 scale mode zscale')
                    os.system('xpaset -p ds9 zoom to fit')
                else:
                    print 'The psf you have given is NOT exists!!!'
                    pass
                write = raw_input("Do you REALLY need this psf? ('y' or," \
                                  "'n' or press any key to continue) ") 
                if write == 'y':
                    TmpPsfList.append(element)
                    fff = open('psflist.list', 'ab')
                    fff.writelines([str(element), '\n'])
                    if(c.galcut):
                        if UpdateCounter:
                            c.ValueS.append(element)
                            fwithpsf = open('CatWithPsf.cat', 'ab')
                            for v in c.ValueS:
                                fwithpsf.writelines([str(v), ' '])
                            fwithpsf.writelines(['\n'])
                            fwithpsf.close()
                            UpdateCounter = 0
                    fff.close()
                if write == 'n':
                    TmpPsfList.append(element)
                    try:
                        os.remove(element)
                    except:
                        pass
                else:
                    pass
            for element in TmpPsfList:
                try:
                    PsfList.remove(element)
                except:
                    pass
            TmpPsfList = []
            fi = raw_input("Finished? ('1' to finish, any other key to continue) ") 
            if fi == '0' or fi == '1':
                finish = int(fi)
                if finish == 1:
                    for element in PsfList:
                        try:
                            os.remove(element)
                        except:
                            pass
            else:
                finish = 0
    else:
        for element in c.psff:
            fff = open('psflist.list', 'ab')
            fff.writelines([str(element), '\n'])
        fff.close()
def SExtractorConf():
    SEx_DETECT_MINAREA = raw_input("DETECT_MINAREA (6) >>> ")
    try:
        c.SEx_DETECT_MINAREA = float(SEx_DETECT_MINAREA)
        c.SEx_DETECT_MINAREA = int(c.SEx_DETECT_MINAREA)
    except:
        c.SEx_DETECT_MINAREA = 6
    SEx_DETECT_THRESH = raw_input('DETECT_THRESH (1.5) >>> ')
    try:
        c.SEx_DETECT_THRESH = float(SEx_DETECT_THRESH)
    except:
        c.SEx_DETECT_THRESH = 1.5
    SEx_ANALYSIS_THRESH = raw_input('ANALYSIS_THRESH (1.5) >>> ')
    try:
        c.SEx_ANALYSIS_THRESH = float(SEx_ANALYSIS_THRESH)
    except:
        c.SEx_ANALYSIS_THRESH = 1.5
    SEx_FILTER = raw_input('FILTER (Y/N) >>> ')
    while SEx_FILTER != 'Y' and SEx_FILTER != 'N' and SEx_FILTER != '':
        SEx_FILTER = raw_input('FILTER (Y/N) >>> ')
    if len(SEx_FILTER) == 0:
        c.SEx_FILTER = 'Y'
    else:
        c.SEx_FILTER = SEx_FILTER
    print 'Available options for convolve filter are gauss_1.5_3x3.conv(1) '\
          'gauss_2.0_3x3.conv(2) gauss_2.0_5x5.conv(3) gauss_2.5_5x5.conv(4) '\
          'gauss_3.0_5x5.conv(5) gauss_3.0_7x7.conv(6) gauss_4.0_7x7.conv(7) '\
          'gauss_5.0_9x9.conv(8) default(0)'
    SEx_FILTER_NAME  = raw_input('FILTER_NAME (default.conv) >>> ')
    if len(SEx_FILTER_NAME) == 0 or SEx_FILTER_NAME == '0':
        c.SEx_FILTER_NAME = 'default.conv'
    elif SEx_FILTER_NAME == '1':
        c.SEx_FILTER_NAME = 'gauss_1.5_3x3.conv'
    elif SEx_FILTER_NAME == '2':
        c.SEx_FILTER_NAME = 'gauss_2.0_3x3.conv'
    elif SEx_FILTER_NAME == '3':
        c.SEx_FILTER_NAME = 'gauss_2.0_5x5.conv'
    elif SEx_FILTER_NAME == '4':
        c.SEx_FILTER_NAME = 'gauss_2.5_5x5.conv'
    elif SEx_FILTER_NAME == '5':
        c.SEx_FILTER_NAME = 'gauss_3.0_5x5.conv'
    elif SEx_FILTER_NAME == '6':
        c.SEx_FILTER_NAME = 'gauss_3.0_7x7.conv'
    elif SEx_FILTER_NAME == '7':
        c.SEx_FILTER_NAME = 'gauss_4.0_7x7.conv'
    elif SEx_FILTER_NAME == '8':
        c.SEx_FILTER_NAME = 'gauss_5.0_9x9.conv'
    SEx_DEBLEND_NTHRESH = raw_input('DEBLEND_NTHRESH (32) >>> ')
    try:
        c.SEx_DEBLEND_NTHRESH = float(SEx_DEBLEND_NTHRESH)
        c.SEx_DEBLEND_NTHRESH = int(c.SEx_DEBLEND_NTHRESH)
    except:
        c.SEx_DEBLEND_NTHRESH = 32
    SEx_DEBLEND_MINCONT = raw_input('DEBLEND_MINCONT (0.005) >>> ')
    try:
        c.SEx_DEBLEND_MINCONT = float(SEx_DEBLEND_MINCONT)
    except:
        c.SEx_DEBLEND_MINCONT = 0.005
    SEx_PHOT_FLUXFRAC = raw_input('PHOT_FLUXFRAC (0.5) >>> ')
    try:
        c.SEx_PHOT_FLUXFRAC = float(SEx_PHOT_FLUXFRAC)
    except:
        c.SEx_PHOT_FLUXFRAC = 0.5
    SEx_pix_scale_disp = 'PIXEL_SCALE (' + str(c.SEx_PIXEL_SCALE) + ') >>> '
    SEx_PIXEL_SCALE = raw_input(SEx_pix_scale_disp)
    try:
        c.SEx_PIXEL_SCALE = float(SEx_PIXEL_SCALE)
    except:
        c.SEx_PIXEL_SCALE = c.pixelscale
    SEx_SEEING_FWHM = raw_input('SEEING_FWHM (0.11) >>> ')
    try:
        c.SEx_SEEING_FWHM = float(SEx_SEEING_FWHM )
    except:
        c.SEx_SEEING_FWHM = c.pixelscale * 3.37
    SEx_BACK_SIZE = raw_input('BACK_SIZE (64) >>> ')
    try:
        c.SEx_BACK_SIZE = float(SEx_BACK_SIZE)
        c.SEx_BACK_SIZE = int(c.SEx_BACK_SIZE)
    except:
        c.SEx_BACK_SIZE = 64
    SEx_BACK_FILTERSIZE = raw_input('BACK_FILTERSIZE (3) >>> ')
    try:
        c.SEx_BACK_FILTERSIZE = float(SEx_BACK_FILTERSIZE)
        c.SEx_BACK_FILTERSIZE = int(c.SEx_BACK_FILTERSIZE)
    except:
        c.SEx_BACK_FILTERSIZE = 3
    SEx_BACKPHOTO_TYPE = raw_input('BACKPHOTO_TYPE (G)LOBAL/(L)OCAL) >>> ')
    while SEx_BACKPHOTO_TYPE != 'G' and SEx_BACKPHOTO_TYPE != 'L' \
          and SEx_BACKPHOTO_TYPE != '':
        SEx_BACKPHOTO_TYPE = raw_input('BACKPHOTO_TYPE (G)LOBAL/(L)OCAL) >>> ')
    if len(SEx_BACKPHOTO_TYPE) == 0:
        c.SEx_BACKPHOTO_TYPE = 'GLOBAL'
    elif SEx_BACKPHOTO_TYPE == 'G':
        c.SEx_BACKPHOTO_TYPE = 'GLOBAL'
    elif SEx_BACKPHOTO_TYPE == 'L':
        c.SEx_BACKPHOTO_TYPE = 'LOCAL'
    if c.SEx_BACKPHOTO_TYPE == 'LOCAL':
        SEx_BACKPHOTO_THICK = raw_input('BACKPHOTO_THICK (24) >>> ')
    try:
        c.SEx_BACKPHOTO_THICK = float(SEx_BACKPHOTO_THICK)
	c.SEx_BACKPHOTO_THICK = int(c.SEx_BACKPHOTO_THICK)
    except:
        c.SEx_BACKPHOTO_THICK = 24
    SEx_WEIGHT_TYPE = raw_input('WEIGHT_TYPE (MAP_RMS) >>> ')
    c.SEx_WEIGHT_TYPE = SEx_WEIGHT_TYPE

def run_SExtractorConf(option, opt, value, parser):
    try:
        SExtractorConf()
    except:
        raise OptionValueError("failure in SextractorConf()")
    return

def rm_sex_cata(option, opt, value, parser):
    if exists(sex_cata):
        os.remove(sex_cata)
    else:
        pass
    return

def run_test(option, opt, value, parser):
    print "using directory " + value + "for output\n"
    c.center_deviated = 0
    c.starthandle = 0
    FindAndFit()
    main()
    if c.crashhandler:
        c.starthandle = 1
        os.system('mv restart.cat CRASH.CAT')
        c.clus_cata = 'CRASH.CAT' 
        main()
    
    sys.exit(0)
    return

if __name__ == '__main__':
    c.FirstCreateDB = 1 #Won't create table when c.FirstCreateDB=0
    c.VERSION = 3.0
    c.FILTER = 'UNKNOWN'
    c.Filter = 'UNKNOWN'
    c.SEx_DETECT_MINAREA = 6
    c.SEx_DETECT_THRESH = 1.5
    c.SEx_ANALYSIS_THRESH = 1.5
    c.SEx_FILTER = 'Y'
    c.SEx_FILTER_NAME = 'default.conv'
    c.SEx_DEBLEND_NTHRESH = 32
    c.SEx_DEBLEND_MINCONT = 0.005
    c.SEx_PHOT_FLUXFRAC = 0.5
    c.SEx_PIXEL_SCALE = c.pixelscale
    c.SEx_SEEING_FWHM = c.pixelscale * 3.37 
    c.SEx_BACK_SIZE = 64
    c.SEx_BACK_FILTERSIZE = 3
    c.SEx_BACKPHOTO_TYPE = 'GLOBAL'
    c.SEx_BACKPHOTO_THICK = 24
    c.SEx_WEIGHT_TYPE = 'DECIDE'
    sex_cata = c.sex_cata

    def FindAndFit():
        if c.findandfit == 1:
            if c.psfselect > 2:
                #magmin = raw_input("Enter Minimum Magnitude >>> ")
                magmin = raw_input("What is faintest magnitude galaxy that you want to fit? >>> ")     #modified by abhishek
                try:
                    magmin = float(magmin) * 1.0
                except:
                    magmin = 9999
                #magmax = raw_input("Enter Maximum Magnitude >>> ")
                magmax = raw_input("What is brightest magnitude galaxy that you want to fit? >>> ")    #modified by abhishek
                try:
                    magmax = float(magmax) * 1.0
                except:
                    magmax = -9999
                stargal = raw_input("Enter star-galaxy classification "\
                                    "(1 for star and 0 is galaxy (0.8 is "\
                                    " a good number) >>> ")
                redshift = raw_input("Enter redshift, if you know >>> ")
                try:
                    redshift = float(redshift)*1.0
                except:
                    redshift = 9999
            elif c.psfselect <= 2:   
                try: 
                    magmin = c.maglim[0]
                    magmax = c.maglim[1] 
                    stargal = c.stargal
                    redshift = 9999
                except:
                    magmin = 9999
                    magmax = -9999
                    stargal = 0.8
                    redshift = 9999                  
            NewClusCata = open(c.clus_cata,'w')
            if redshift == 9999:
                NewClusCata.writelines(['gal_id ra1 dec1 mag\n'])
            else:
                NewClusCata.writelines(['gal_id ra1 dec1 mag z\n'])
            for lineSex in open(sex_cata, 'r'):
                ValueSex = lineSex.split()
                try:
                    if float(ValueSex[17]) > magmax and \
                       float(ValueSex[17]) < magmin and \
                       float(ValueSex[16]) < float(stargal):
                        if redshift == 9999:
                            NewClusCata.writelines([str(ValueSex[0]), ' ', \
                                       str(float(ValueSex[3]) / 15.0), ' ', \
                                       str(float(ValueSex[4])), ' ',\
                                       str(ValueSex[17]), '\n'])
                        else:
                            NewClusCata.writelines([str(ValueSex[0]), ' ', \
                                str(float(ValueSex[3]) / 15.0), ' ', \
                                str(float(ValueSex[4])), \
                                ' ', str(redshift), ' ', \
                                str(ValueSex[17]), '\n'])
                except:
                    pass
            NewClusCata.close()
            c.searchrad = '0.05arc'
        else:
            pass

    # Note I set defaults here and call them in creating the options. This is so 
    # I can use them in determining whether to retain the default or not.

    usage = "Usage: pymorph [--edit-conf[-e]] [--with-psf [-p]] [--force[-f]] "\
        "[--help[-h]] [--lmag] [--umag] [--ln] [--un] [--lre] [--ure] "\
        "[--lrd] [--urd] [--with-in] [--with-filter] [--with-db] "\
        "[--with-area]  [--no-mask] [--norm-mask] [--with-sg] [--bdbox] "\
        "[--bbox] [--dbox] [--test [-t]] [--devauc] [--outdir] [--datadir]"
    parser = OptionParser(usage=usage)
    parser.add_option("-e", "--edit_conf", action="callback", 
                      callback=run_SExtractorConf, 
                      help="runs SExtractor configuration")
    parser.add_option("-f", "--force", action="callback", 
                      callback=rm_sex_cata,
                      help="removes SExtractor catalog")
    parser.add_option("-p", "--with-psf", action="store", type="int",
                      dest="WhichPsf",default = 0,
                      help="Nearest/farthest PSF")
    parser.add_option("-t", "--test", action="callback", callback=run_test, 
                      type="string", help="runs the test instance-OVERRIDES ALL OTHER INPUT-User must supply a directory for output")
    parser.add_option("--lmag", action="store", type="float",
                      dest="LMag",default = 500.0, 
                      help="lower magnitude cutoff")
    parser.add_option("--umag", action="store", type="float",
                      dest="UMag", default = -500.0,
                      help="upper magnitude cutoff")
    parser.add_option("--ln", action="store", type="float",
                      dest="LN", default = 0.1, help="Lower Sersic")
    parser.add_option("--un", action="store", type="float",
                      dest="UN", default = 20.0, help="Upper Sersic")
    parser.add_option("--lre", action="store", type="float",
                      dest="LRe", default = 0.0,
                      help="Lower Bulge Radius")
    parser.add_option("--ure", action="store", type="float",
                      dest="URe", default = 500.0,
                      help="Upper Bulge Radius")
    parser.add_option("--lrd", action="store", type="float",
                      dest="LRd", default = 0.0, help="Lower Disk Radius")
    parser.add_option("--urd", action="store", type="float",
                      dest="URd", default = 500.0, help="Upper Disk Radius")
    parser.add_option("--with-in", action="store", type="float",
                      dest="avoidme", default = 50.0, help="avoid me!")
    parser.add_option("--with-filter", action="store", type="string", 
                      dest="Filter", default = 'UNKNOWN', help="Filter used")
    parser.add_option("--with-db", action="store", type="string",
                       dest="database", default = 'UNKNOWN', 
                      help="database used")
    parser.add_option("--with-area", action="store", type="float",
                      dest="AreaOfObj", default = 40.0, 
                      help="min Area of psf for selection")
    parser.add_option("--no_mask", action="store_true", dest="NoMask",
                      default = False, help="turns off masking")
    parser.add_option("--norm_mask", action="store_true", dest="NormMask",
                      default = False, help="turns on Normal masking")
    parser.add_option("--with-sg", action="store", type="float",
                      dest="StarGalProb", default = 0.9, 
                      help="for psf identification")
    parser.add_option("--bdbox", action="store_true", dest="bdbox",
                      default = False, help="turns on bdbox")
    parser.add_option("--bbox", action="store_true", dest="bbox",
                      default = False, help="turns on bbox")
    parser.add_option("--dbox", action="store_true", dest="dbox",
                      default = False, help="turns on dbox")
    parser.add_option("--devauc", action="store_true", dest="devauc",
                      default = False, 
                      help="turns on DeVacouleur's bulge fitting (sersic index of bulge frozen at 4.0 for fit)")
    parser.add_option("-o","--outdir", action="store", type="string",
                      dest="outdir", default = os.getcwd()+'/',
                      help="path to directory that will contain all output. MUST end in '/'")
    parser.add_option("-d","--datadir", action="store", dest="datadir", 
                      type="string", default = os.getcwd()+'/',
                      help="path to directory containing all input. MUST end in '/'")
    parser.add_option("--with-host", action="store", type="string",
                       dest="host", default = 'localhost', 
                      help="mysql host used")
    
    # parses command line aguments for pymorph
    (options, args) = parser.parse_args()

    # compares command line arguments to defaults,
    # if default, only stores value if value is not set in config file
    # if not default, then stores to config regardless of config file value.
    for key, value in options.__dict__.items():
        if  parser.defaults[key] ==  value:
            if not hasattr(c, key):
                setattr(c,key, value)
        else:
            setattr(c,key, value)

    if c.Filter == 'UNKNOWN':
        pass
    else:
        c.FILTER = c.Filter


    # now change dir to the 
    thisdir = os.getcwd()
    print "thisdir is ", thisdir
    print "outdir is ", c.outdir
    os.chdir(c.outdir)
    try:
        if(c.repeat == False and c.galcut == False):
            img = pyfits.open(c.datadir + c.imagefile)
            c.imagedata = img[0].data
            c.HeAdEr0 = img[0].header
            img.close()
            ut.CheckHeader(c.HeAdEr0)
    except IOError, (errno, strerror):
        print imagefile, "I/O error(%s): %s" % (errno, strerror)
        os._exit(0)

    if exists(sex_cata):
        pass
    elif(c.galcut == False):
        print 'The SExtractor catalogue for your frame is NOT found. ' \
              'One is being made using the default values. It is always '\
              'recommended to make SExtractor catalogue by YOURSELF as '\
              'the pipeline keeps the sky value at the SExtractor value '\
              'during the decomposition.'
        if exists(c.datadir + c.whtfile):
            RunSex(c.datadir + c.imagefile, c.datadir + c.whtfile, 'None', 9999, 9999, 0)
	    SexShallow(c.datadir + c.imagefile, c.datadir + c.whtfile, 'None', 9999, 9999, 0)
        else:
            RunSex(c.datadir + c.imagefile, 'None', 'None', 9999, 9999, 0)
	    SexShallow(c.datadir + c.imagefile, 'None', 'None', 9999, 9999, 0)
    def runpsfselect():
        if(c.galcut):   #Given galaxy cutouts
            obj_file = open(c.datadir + c.clus_cata,'r') 
            pnames = obj_file.readline().split() 
            c.ValueS = []
            for v in pnames:
                c.ValueS.append(v)
            c.ValueS.append('star')
            fwithpsf = open('CatWithPsf.cat', 'ab')
            for v in c.ValueS:
                fwithpsf.writelines([str(v), ' '])
            fwithpsf.writelines(['\n'])
            fwithpsf.close()
            pdb = {}                        #The parameter dictionary
            for line_j in obj_file:
                try:
                    values = line_j.split()
                    c.ValueS = []
                    for v in values:
                        c.ValueS.append(v)
                    k = 0
                    for pname in pnames:
                        pdb[pname] = values[k]
                        k += 1
                    try:
                        gal_id = pdb["gal_id"]
                    except:
                        try:
                            gal_id = pdb["gimg"][:-5]
                        except:
                            print "No image or gal_id found in the object" \
                                  "catalogue. Exiting"
                            os._exit(0)
                    try:
                        gimg = pdb["gimg"]    #Galaxy cutout
                    except:
                        if exists(c.datadir + 'I' + str(c.rootname) + '_' \
                                   + str(gal_id) + '.fits'):
                            gimg = 'I' + str(c.rootname) + '_' \
                                    + str(gal_id) + '.fits'
                        elif exists(c.datadir + str(gal_id) + '.fits'):
                            gimg = str(gal_id) + '.fits'
                        else:
                            print "No image found. Exiting"
                            os._exit(0)
                    try:
                        wimg = pdb["wimg"]   #Weight cut
                    except:
                        if exists(c.datadir + 'W' + str(c.rootname) + '_' + \
                                  str(gal_id) + '.fits'):
                            wimg = 'W' + str(c.rootname) + '_' + \
                                   str(gal_id) + '.fits'
                        else:
                            wimg = 'None'
                    GiMg = pyfits.open(c.datadir + gimg)
                    headerGiMg = GiMg[0].header
                    if (headerGiMg.has_key('GAIN')):
                        c.SEx_GAIN = headerGiMg['GAIN']
                    else:
                        c.SEx_GAIN = 1
                    GiMg.close()
                    if exists(sex_cata): 
                        pass
                    else:
                        RunSex(c.datadir + gimg,c.datadir + wimg, 'None', 9999, 9999, 0)
                    try:
                        selectpsf(c.datadir + gimg,sex_cata)
                    except:
                        pass
                    if os.access(sex_cata, os.F_OK):
                        os.remove(sex_cata)
                except:
                    pass
            obj_file.close()  
            AskForUpdate = raw_input("Do you want to update the clus_cata? " \
                            "('y' for yes) ")
            if AskForUpdate == 'y':
                cmd = 'mv CatWithPsf.cat ' +str(c.clus_cata)
                os.system(cmd)
            else:
                pass
        else:
            selectpsf(c.datadir + c.imagefile,sex_cata)
#The old function for psfselect = 2
#    if c.psfselect == 2:
#        c.center_deviated = 0
#        c.starthandle = 0
#        os.system('ds9 &')
#        time.sleep(2)
#        runpsfselect()
#        os.system('xpaset -p ds9 quit')
#        c.psflist = '@psflist.list'
#        FindAndFit()
#        main()
#        if c.crashhandler:
#            c.starthandle = 1
#            os.system('mv restart.cat CRASH.CAT')
#            c.clus_cata = 'CRASH.CAT' 
#            main()
#New function for psfselect=2 non-interactive for webservice
    if c.psfselect == 2:
        c.Interactive = 0
        c.center_deviated = 0
        c.starthandle = 0
        runpsfselect()
        c.psflist = c.datadir + '@psflist.list'
        FindAndFit()
        main()
        if c.crashhandler:
            c.starthandle = 1
            os.system('mv restart.cat CRASH.CAT')
            c.clus_cata = 'CRASH.CAT' 
            main()

#The old function for psfselect=1
    elif c.psfselect == 1:
        c.Interactive = 1
        os.system('ds9 &')
        time.sleep(2)
        runpsfselect()
        os.system('xpaset -p ds9 quit')
#new function for psfselect=1 non-interactive for webservice (abhishek rawat)
#    elif c.psfselect == 1:
#        c.Interactive = 0
#        runpsfselect()
        



    elif c.psfselect == 0:
        c.center_deviated = 0
        c.starthandle = 0
        FindAndFit()
        main()
        if c.crashhandler:
            c.starthandle = 1
            os.system('mv restart.cat CRASH.CAT')
            c.clus_cata = 'CRASH.CAT' 
            main()

    os.chdir(thisdir)
