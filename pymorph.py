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
import pymorphutils_rewriting as ut
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
    sex_cata = c.sex_cata

    #Size parameters
    ReSize, VarSize, Square, FracRad, FixSize = ut.GetSizeInfo()

    #Initialize index.html
    if exists('index.html'):
        pass
    else:
        indexfile = open('index.html', 'w')
        indexfile.writelines(['<HTML>\n<BODY>\n'])
        indexfile.writelines(['</BODY></HTML>'])
        indexfile.close()

    # Opens files to write output
    f_cat = open(c.out_cata, 'w')
    f_failed = open('restart.cat', 'w')

    # Remove some intermediate files
    for xf in ['TmpElliMask.fits', 'TmpElliMask1.fits']:
        if os.path.exists(xf):
            os.system('rm -f %s'%xf)

    #Initializing psf array. ie. creating c.psflist from file 
    ut.PsfArr() #Now c.psflist has a list

    # Updating PSF header
    if c.decompose:
        for psfelement in c.psflist:
            ut.UpdatePsfRaDec(psfelement)

    # Model for fitting
    ComP = ut.GetModel()

    # Writing csv header and finding the output parameters
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


    # The file contains the objects of interest 
    obj_file = open(os.path.join(c.datadir, clus_cata), 'r')  
    pnames = obj_file.readline().split() #The names of the parameters given 
                                         #in the first line in the clus_cata

    # writing a input catalogue (restart.cat) for failed objects
    for p in pnames:
        f_failed.writelines(['%s '%p])
    f_failed.writelines(['flag \n'])

    #The parameter dictionary
    pdb = {} 

    #For getting psf in the case of unknown ra
    c.psfcounter = 0 

    for line_j in obj_file:
        # declare the flag
        c.Flag = 0
       
        # Assign parameter values
        alpha1, alpha2, alpha3, delta1, delta2, delta3, z, bxcntr, bycntr,\
        UserGivenPsf, UserGivenSky, gimg, wimg, ximg, yimg, RaDecInfo = \
        ut. GetInputParams(pnames, line_j)
 
        # Find what is input images
        gimg, outimage, wimg, config_file, c.pfile, maskimage, confile = \
                                  DecisionMaker(gimg, wimg, ximg, yimg)

	# Crashhandling starts
	SetCrashHandler()

        # Set coordinates               
	if(alpha1 == -9999 or delta1 == -9999):
	    alpha_j = -9999
	    delta_j = -9999
	else:
	    alpha_j = ut.HMSToDeg(alpha1, alpha2, alpha3)
	    delta_j = ut.DMSToDeg(delta1, delta2, delta3)
	# Determine Search Radius 
        SeaDeg, SeaPix = ut.GetSearchRad(RaDecInfo)

        # Get data
        ut.GetImageData()

        #Reading weigh files
        if(c.repeat == False and c.galcut == False):
            TY, TX = c.imagedata.shape
            ut.GetWhtImage()


	# Getting the nearest object from the sextractor catalog
        good_object = ut.FindSexObj(sex_cata, RaDecInfo, SeaDeg, SeaPix,
                                    alpha_j, delta_j, ximg, yimg)

	if len(good_object) <1:
	    # No suitable target found
	    print "NO TARGET FOUND!!!!"
	    good_object = ' 9999  9999 9999  9999 9999  9999 9999  9999 9999  9999 9999  0 9999  9999 9999  9999 9999 9999 9999\n'
	    c.Flag = SetFlag(c.Flag,GetFlag('NO_TARGET'))  

	# now fit best object            
	try:
            # Since we found the object making some initial setup
            c.run = 1 #run =1 if pipeline runs sucessfuly
	    ut.WriteError('\n\n###########   ' + str(gal_id) + \
			      '   ###########\n')
            # Adding initial setup to the flag
	    if c.repeat:
		c.Flag = SetFlag(c.Flag,GetFlag('REPEAT'))
	    if c.fitting[0] and 'bulge' in ComP:
		c.Flag = SetFlag(c.Flag,GetFlag('FIT_BULGE_CNTR'))
	    if c.fitting[1] and 'disk' in ComP:
		c.Flag = SetFlag(c.Flag,GetFlag('FIT_DISK_CNTR'))
	    if c.fitting[2]:
		c.Flag = SetFlag(c.Flag,GetFlag('FIT_SKY'))
	    
            # Getting sextractor parameters
	    values = good_object.split()
	    print "SExtractor ID >>> ", values[0]
	    alpha_s = float(values[3])
	    delta_s = float(values[4])
	    sex_id = values[0]
	    xcntr  = float(values[1])
	    ycntr  = float(values[2])

            # Bunch of sextractor global parameters 
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
          
            # Some fit limits
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

            # Set sky values (c.SexSky & c.GalSky)
            ut.SetSkies(values[0], values[10])

	    
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
	    if not c.repeat:
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
		else:
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
		print 'Center of cutimage and exceed size ', \
		      cut_xcntr, cut_ycntr, ExceedSize
		print 'Full Sizes ', SizeX, SizeY
		if c.galcut and ReSize == 0:
		    pass
		elif ExceedSize:
		    c.Flag = SetFlag(c.Flag,GetFlag('EXCEED_SIZE'))

		# Runs sextractor to find the segmentation map
		RunSegSex(os.path.join(c.datadir, cutimage))

		# This creates ellipse mask and used for 
		# ellipse fitting and casgm
		ElliMaskFunc(cutimage, cut_xcntr, cut_ycntr, \
			     SizeX, SizeY, good_object, 1)
		if c.decompose:
		    ut.HandleEllipseTask(cutimage, cut_xcntr, \
				   cut_ycntr, \
				   SizeX, SizeY, c.SexSky, 0)
		    MaskFunc(cutimage, cut_xcntr, cut_ycntr, \
				     SizeX, SizeY, good_object)
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
		    pass
		    MaskFunc(cutimage, cut_xcntr, cut_ycntr, \
			     SizeX, SizeY, good_object)
		config_file = cfile
		outimage = str(oimg)

	    # Estimates sky parameters
	    try:
		SexySky, SkyYet, SkyMed, SkyMin, SkyQua, \
		SkySig = \
		FindYetSky(os.path.join(c.datadir, cutimage), \
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
		    print "something bad happened (GALFIT)!!!!\n\n"
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
		print "something bad happened (Writing)!!!!\n\n"
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
	    print "something bad happened (general sextractor)!!!!\n\n"
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
            print "something bad happened (general object search)!!!!\n\n"
            print traceback.print_exc()
    f_cat.close()
    f_failed.close()

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
                setattr(c, key, value)
        else:
            setattr(c, key, value)

    # now change dir to the 
    thisdir = os.getcwd()
    print "This dir is > ", thisdir
    print "Output dir is > ", c.outdir
    os.chdir(c.outdir)

    # Read images 
    try:
        if not c.repeat and not c.galcut:
            img = pyfits.open(os.path.join(c.datadir, c.imagefile))
            c.imagedata = img[0].data
            c.HeAdEr0 = img[0].header
            img.close()
            ut.CheckHeader(c.HeAdEr0)
    except IOError, (errno, strerror):
        print imagefile, "I/O error(%s): %s" % (errno, strerror)
        os._exit(0)

    # Generate sextractor catalogs if not exists
    if exists(sex_cata):
        pass
    elif not c.galcut:
        print 'The SExtractor catalogue for your frame is NOT found. ' \
              'One is being made using the default values. It is always '\
              'recommended to make SExtractor catalogue by YOURSELF as '\
              'the pipeline keeps the sky value at the SExtractor value '\
              'during the decomposition.'
        if exists(os.path.join(c.datadir, c.whtfile)):
            RunSex(os.path.join(c.datadir, c.imagefile), 
                   os.path.join(c.datadir , c.whtfile), 
                   'None', 9999, 9999, 0)
	    SexShallow(os.path.join(c.datadir, c.imagefile), 
                       os.path.join(c.datadir, c.whtfile), 
                       'None', 9999, 9999, 0)
        else:
            RunSex(os.path.join(c.datadir, c.imagefile), 
                   'None', 'None', 9999, 9999, 0)
	    SexShallow(os.path.join(c.datadir, c.imagefile), 
                       'None', 'None', 9999, 9999, 0)

    #New function for psfselect=2 non-interactive for webservice
    if c.psfselect == 2:
        c.Interactive = 0
        c.center_deviated = 0
        c.starthandle = 0
        runpsfselect()
        c.psflist = os.path.join(c.datadir, '@psflist.list')
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
