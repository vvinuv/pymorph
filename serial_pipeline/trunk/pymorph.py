#!/home/vinu/software/Python2.5/bin/python
"""PyMorph [Py MOrphological Parameters' Hunter], is a pipeline to find the Morphological parameters of galaxy. Authors: Vinu Vikram , Yogesh Wadadekar, Ajit K. Kembhavi. 2008 Feb"""

import os
import time
from os.path import exists
import sys
from getopt import getopt, GetoptError
import csv
import pyfits
import numpy as n
from ndimage import center_of_mass
sys.path.append('.')
import config as c
from maskfunc import *
from configfunc import *
from ellimaskfunc import *
from outmaskfunc import *
from plotfunc import *
from writehtmlfunc import *
from runsexfunc import *
from casgm import *
from bkgdfunc import *
#from configiter import *
from configbarpoint import *

try:
    from pyraf import iraf
    from fitellifunc import *
except:
    pass

def main():
    imagefile = c.imagefile
    whtfile = c.whtfile
    sex_cata = c.sex_cata
    clus_cata = c.clus_cata
    out_cata = c.out_cata
    try:
        if c.psflist.startswith('@'):
            psffi = open(c.psflist.split('@')[1], 'r')
            c.psflist = []
            for pline in psffi: 
                c.psflist.append(pline.split()[0])
    except:
        pass     
    ReSize = c.size[0]
    try:
        VarSize = c.size[1]
    except:
        if ReSize:
            VarSize = 1
        else:
            VarSize = 0
    try:
        Square = c.size[3]
    except:
        Square = 1
    try:
        FracRad = c.size[2]  
    except:
        FracRad = 20
    if VarSize == 0:
        try:
            FixSize = c.size[4]
        except:
            FixSize = 120
    threshold = c.threshold
    thresh_area = c.thresh_area
    try:
        os.system('rm -f TmpElliMask.fits TmpElliMask1.fits')
    except:
        pass
    if exists('index.html'):
        pass
    else:
        indexfile = open('index.html', 'w')
        indexfile.writelines(['<HTML>\n<BODY>\n'])
        indexfile.writelines(['</BODY></HTML>'])
        indexfile.close()
    try:
        if(c.repeat == False and c.galcut == False):
            image = c.ImAgE
            header0 = c.HeAdEr0
            if (header0.has_key('EXPTIME')):
                EXPTIME = header0['EXPTIME']
            else:
                EXPTIME = -9999
            if (header0.has_key('RDNOISE')):
                RDNOISE= header0['RDNOISE']
            else:
                RDNOISE = -9999
            if (header0.has_key('GAIN')):
                GAIN = header0['GAIN']
            else:
                GAIN = -9999
            if (header0.has_key('NCOMBINE')):
                NCOMBINE= header0['NCOMBINE']
            else:
                NCOMBINE = -9999
            print "imagefile >>> ", imagefile
            TX = image.shape[1]
            TY = image.shape[0]
    except IOError, (errno, strerror):
        print imagefile, "I/O error(%s): %s" % (errno, strerror)
        os._exit(0)
    try:
        if exists(whtfile):
            if(c.repeat == False and c.galcut == False):
                wht = pyfits.open(whtfile)
                weight = wht[0].data
                wht.close()
                print "whtfile >>> ", whtfile
        else:
           print 'No weight image found\n'
    except IOError, (errno, strerror):
        print whtfile, "I/O error(%s): %s" % (errno, strerror)
        pass
    psflist = c.psflist
    def psfradec(element):
        """The function which will update the psf header if the psf files
           are the specified format"""
        try:
            ra1 = float(str(element)[4:6])
            ra2 = float(str(element)[6:8])
            ra3 = float(str(element)[8:10]) + float(str(element)[10]) / 10.0
            dec1 = float(str(element)[11:-10])
            dec2 = float(str(element)[-10:-8])
            dec3 = float(str(element)[-8:-6]) + float(str(element)[-6]) / 10.0
            ra = (ra1 + (ra2 + ra3 / 60.0) / 60.0) * 15.0
            if dec1 < 0.0:
                dec = (dec1 - (dec2 + dec3 / 60.0) / 60.0)
            else:
                dec = (dec1 + (dec2 + dec3 / 60.0) / 60.0)
            try:
                iraf.hedit(element, 'RA_TARG', ra, add= 'yes', verify= 'no', \
                               show='no', update='yes')
                iraf.flpr()
                iraf.hedit(element, 'DEC_TARG', dec, add= 'yes', verify= 'no',\
                       show='no', update='yes')
                iraf.flpr()
            except:
                pass
        except:
            pass
    def failedgalfit(WhichGalaxy):
        f_fail = open("fit.log", "w")
        f_fail.writelines(['-----------------------------------------------'\
                            '------------------------------\n\n'])
        f_fail.writelines(['Input image     : ', str(WhichGalaxy), '\n'])
        f_fail.writelines(['Init. par. file : Failed! :(', '\n'])
        f_fail.writelines(['Restart file    : Failed! :(', '\n'])
        f_fail.writelines(['Output image    : Failed! :(', '\n\n'])
        if 'bulge' in c.components:
            f_fail.writelines([' sersic   : (9999, 9999)   9999   9999   9999'\
                           '    9999   9999   9999', '\n'])
            f_fail.writelines(['              (9999, 9999)   9999    9999'\
                               '    9999    9999   9999   9999', '\n'])
        if 'disk' in c.components:
            f_fail.writelines([' expdisk   : (9999, 9999)   9999   9999'\
                               '    9999   9999   9999', '\n'])
            f_fail.writelines(['              (9999, 9999)   9999    9999'\
                               '    9999   9999   9999', '\n'])
        if 'point' in c.components:
            f_fail.writelines([' gaussian   : (9999, 9999)   9999   9999'\
                               '    9999   9999   9999', '\n'])
            f_fail.writelines(['              (9999, 9999)   9999    9999'\
                               '    9999   9999   9999', '\n'])
        f_fail.writelines([' sky      : [9999, 9999]   9999   9999   '\
                           '9999', '\n'])
        f_fail.writelines(['                             9999   9999   9999\n'])
        f_fail.writelines([' Chi^2 = 9999,  ndof = 9999\n'])
        f_fail.writelines([' Chi^2/nu = 9999\n\n'])
        f_fail.writelines(['-----------------------------------------------'\
                            '------------------------------'])
        f_fail.close()
    def pa(x):
        """ The function which will bring position angle 
         measured by sextrator in the range -90 and 90"""		
        if(float(x)>=0 and float(x)<=180.0): 
            pos_ang = float(x) - 90.0 #position angle
        if(float(x)<0 and float(x)>=-180.0):
            pos_ang = 90.0 - abs(float(x))  #position angle
        if(float(x)>180 and float(x)<=360.0):
            pos_ang = float(x) - 360.0 + 90.0 #position angle
        if(float(x)>=-360 and float(x)<-180.0):
            pos_ang = float(x) + 360.0 - 90.0 #position angle	
        return pos_ang

    def psf_select(alpha_j, delta_j):					
        """This function will select the nearest psf from the psflist.
           The distance is calculated by using the following equation
           d = Sqrt((dec_a - dec_b) ^ 2 + ((ra_a - ra_b) * sin(0.5) * 
           (dec_a - dec_b)) ^ 2.0 )"""
        PsfDistanceDict = {}
        distance = 9999.0
        psffile = 'test.fits'
        psflist = c.psflist
        r = 3.14159265 / 180.0
        for element in psflist:
            p=pyfits.open(element)
            header = p[0].header
            if (header.has_key('RA_TARG')):
                ra = header['RA_TARG']
            else:
                ra = 9999
            if (header.has_key('DEC_TARG')):
                dec= header['DEC_TARG']
            else:
                dec= 9999
            p.close()
#		d = sqrt((ra - alpha_j) ** 2.0 + (dec - delta_j) ** 2.0)
#            d = n.arccos(n.cos((90.0 - delta_j) * r) * n.cos((90.0 - dec) *\
#                r) + n.sin((90.0 - delta_j) * r) *  n.sin((90.0 - dec) * r) * \
#                n.cos((alpha_j - ra) * r))
#            d = n.sqrt((delta_j - dec)**2.0 + ((alpha_j-ra)*n.sin((0.5) *\
#                (delta_j+dec)))**2.0)
            d = n.sqrt((delta_j - dec)**2.0 + ((alpha_j - ra) * \
                n.cos(delta_j * r))**2.0)
            PsfDistanceDict[element] = d
            #print 'alp dec alpsf decpsf d', alpha_j, delta_j, ra, dec, d
#            if(d < distance):
#                psffile = element
#                distance = d
        ItemS = PsfDistanceDict.items()
        ItemS = [(v, k) for (k, v) in ItemS]
        ItemS.sort()
        ItemS = [(k, v) for (v, k) in ItemS]
        psffile = ItemS[c.WhichPsf][0]
        distance = ItemS[c.WhichPsf][1]
        return psffile, distance 
    def CrashHandlerToRemove(gal_id):
        RemoveMe = 'I' + str(c.rootname) + '_' + str(gal_id) + '.fits ' +\
                   'I' + str(c.rootname) + '_' + str(gal_id) + '.con ' +\
                   'W' + str(c.rootname) + '_' + str(gal_id) + '.fits ' +\
                   'E_I' + str(c.rootname) + '_' + str(gal_id) + '.txt ' +\
                   'EM_I' + str(c.rootname) + '_' + str(gal_id) + '.fits ' +\
                   'G_I' + str(c.rootname) + '_' + str(gal_id) + '.in ' + \
                   'M_I' + str(c.rootname) + '_' + str(gal_id) + '.fits ' +\
                   'OE_I' + str(c.rootname) + '_' + str(gal_id) + '.txt ' +\
                 'OEM_O_I' + str(c.rootname) + '_' + str(gal_id) + '.fits ' +\
                 'O_I' + str(c.rootname) + '_' + str(gal_id) + '.fits ' +\
                 'Tmp* ' + \
                 'SO_I' + str(c.rootname) + '_' + str(gal_id) + '.fits '
        f_R_crash =  'R_I' + str(c.rootname) + '_' + str(gal_id) + '_1.html'
        f_R_cra = open(f_R_crash, 'w')
	P_cra = 'P_I' + str(c.rootname) + '_' + str(gal_id) + '.png'
        P_new = 'P_I' + str(c.rootname) + '_' + str(gal_id) + '_1.png'
        for line_crash in open('R_I' + str(c.rootname) + '_' + str(gal_id) + \
                               '.html', 'r'):
            line_crash2wri = line_crash.replace(P_cra, P_new)
	    f_R_cra.write(line_crash2wri)
        f_R_cra.close() 
	CmdToRename = 'mv ' + 'P_I' + str(c.rootname) + '_' + str(gal_id) + \
	               '.png ' + \
		       'P_I' + str(c.rootname) + '_' + str(gal_id) + '_1.png '
	os.system(CmdToRename)
        LinuxCommand = 'rm -f ' + str(RemoveMe) + 'R_I' + str(c.rootname) + \
                       '_' + str(gal_id) + '.html'
        os.system(LinuxCommand)
    def GetFlag(flagname):
        FlagDict = dict([('REPEAT', 0),
                         ('FIT_BULGE_CNTR', 1),
                         ('FIT_DISK_CNTR', 2),
                         ('FIT_SKY', 3),
                         ('EXCEED_SIZE', 4),
                         ('ELLIPSE_FAIL', 5),
                         ('CASGM_FAIL', 6),
                         ('GALFIT_FAIL', 7),
                         ('PLOT_FAIL', 8),
                         ('FIT_BULGE', 9),
                         ('FIT_DISK', 10),
                         ('FIT_POINT', 11),
                         ('NEIGHBOUR_FIT', 12),
                         ('LARGE_CHISQ', 13),
                         ('SMALL_GOODNESS', 14),
                         ('FAKE_CNTR', 15),
                         ('BULGE_AT_LIMIT', 16),
                         ('DISK_AT_LIMIT', 17),
                         ('ASYM_NOT_CONV', 18),
                         ('ASYM_OUT_FRAME', 19),
                         ('BACK_FAILED', 20)])
        return FlagDict[flagname]
    def isset(flag, bit):
        """Return True if the specified bit is set in the given bit mask"""
        return (flag & (1 << bit)) != 0
    try:
        ComP = c.components
    except:
        ComP = ['bulge', 'disk']
    if len(ComP) == 0:
        ComP = ['bulge', 'disk']
#weight = where(weight1 > 0, 1.0 / sqrt(weight1), 0.0)
    if exists('result.csv'):
        pass
    else:
        f_res = open("result.csv", "ab")
        writer = csv.writer(f_res)
        if(c.decompose):
            ParamToWrite = ['Name_0', 'ra__1', 'dec__2', 'z_3', 'mag_auto_4',\
                            'magerr_auto_5', 'Ie_6', 'Ie_err_7', 're_pix_8',\
                            're_err_pix_9', 're_kpc_10', 're_err_kpc_11', \
                            'n_12', 'n_err_13', 'AvgIe_14', 'AvgIe_err_15', \
                            'eb_16', 'eb_err_17', 'bboxy_18', 'bboxy_err_19',\
                            'Id_20', 'Id_err_21', 'rd_pix_22', \
                            'rd_err_pix_23', 'rd_kpc_24', 'rd_err_kpc_25', \
                            'ed_26', 'ed_err_27', 'dboxy_28', 'dboxy_err_29',\
                            'BD_30', 'BT_31', 'Point_32', 'Point_err_33', \
                            'Pfwhm_34', 'Pfwhm_kpc_35', 'chi2nu_36', \
                            'Goodness_37', 'run_38', 'C_39', 'C_err_40',\
                            'A_41', 'A_err_42', 'S_43', 'S_err_44', 'G_45', \
                            'M_46', 'SexSky_47', 'GalSky_48', 'dis_modu_49', \
                            'distance_50', 'fit_51', 'flag_52', \
                            'HalfRadius_53', 'BarMag_54', 'BarMagErr_55', \
                            'BarRePix_56', 'BarRePixErr_57', 'BarReKpc_58',\
                            'BarReKpcErr_59', 'BarIndex_60', 'BarIndexErr_61',\
                            'BarEll_62', 'BarEllErr_63', 'BarBoxy_64', \
                            'Comments_65']
#            if 'bulge' in ComP:
#                for bulgecomp in ['Ie','Ie_err','re(pixels)','re_err(pixels)',\
#                                  're(kpc)', 're_err(kpc)' ,'n', 'n_err']
#                    ParamToWrite.append(bulgecomp)
#            if 'disk' in ComP:
#                for diskcomp in ['Id','Id_err','rd(pixels)','rd_err(pixels)', \
#                                 'rd(kpc)', 'rd_err(kpc)']: 
#                    ParamToWrite.append(diskcomp)
#            if 'bulge' in ComP and 'disk' in ComP:
#                ParamToWrite.append('BD')
#                ParamToWrite.append('BT')
#            if 'point' in ComP:
#                ParamToWrite.append('Point')
#                ParamToWrite.append('Point_err')
#            for otherparam in ['chi2nu', 'Goodness', 'run', 'C', 'C_err', 'A',\
#                               'A_err', 'S', 'S_err', 'G', 'M', 'distance', \
#                               'fit', 'flag', 'Comments']:
#                ParamToWrite.append(otherparam)
            writer.writerow(ParamToWrite)
        else:
            writer.writerow(['Name','ra_','dec_','z', 'mag_auto', \
	                 'magerr_auto', 'C', \
                         'C_err', 'A', 'A_err', 'S', 'S_err', 'G', 'M', \
                         'flag', 'HalfRadius', 'Comments'])
        f_res.close()
    f_cat = open(out_cata,'w')
    f_failed = open('restart.cat', 'w')
    obj_file = open(clus_cata,'r')  #The file contains the objects of interest
    pnames = obj_file.readline().split() #The names of the parameters given 
                                         #in the first line in the clus_cata
    for FailedParam in pnames:
        f_failed.writelines([str(FailedParam), ' '])
    f_failed.writelines(['flag \n'])
    pdb = {}                        #The parameter dictionary
    psfcounter = 0                  #For getting psf in the case of unknown ra
    for line_j in obj_file:
	RaDecInfo = 0
        try:
            values = line_j.split()
            k = 0
            for pname in pnames:
                pdb[pname] = values[k]
                k += 1
            try:
                gal_id = pdb["gal_id"]
            except:
                try:
                    gal_id = pdb["gimg"][:-5] #id will be filename without .fits
                except:
                    print "No image or gal_id found in the object catalogue." \
                          "Exiting"
                    os._exit(0)
            try:
                alpha1 = float(pdb["ra1"])
		RaDecInfo = 1
            except:
                alpha1 = -9999
            try:
                alpha2 = float(pdb["ra2"])
            except:
                alpha2 = 0
            try:
                alpha3 = float(pdb["ra3"])
            except:
                alpha3 = 0
            try:
                delta1 = float(pdb["dec1"])
		RaDecInfo = 1
            except:
                delta1 = -9999
            try:
                delta2 = float(pdb["dec2"])
            except:
                delta2 = 0
            try:
                delta3 = float(pdb["dec3"])
            except:
                delta3 = 0
            try:
                z = float(pdb["z"])
            except:
                z = 9999
            try:
                gimg = pdb["gimg"]    #Galaxy cutout
            except:
                if exists('I' + str(c.rootname) + '_' + str(gal_id) + '.fits'):
                    gimg = 'I' + str(c.rootname) + '_' + str(gal_id) + '.fits'
                elif exists(str(gal_id) + '.fits'): 
                    gimg = str(gal_id) + '.fits'
                else:
                    gimg = 'None'
            try:
                wimg = pdb["wimg"]   #Weight cut
            except:
                if exists('W' + str(c.rootname) + '_' + str(gal_id) + '.fits'):
                    wimg = 'W' + str(c.rootname) + '_' + str(gal_id) + '.fits'
                else:
                    wimg = 'None'
            try:
                cfile = pdb["cfile"]  #GALFIT configuration file
            except:
                if(c.repeat == True and c.galcut == False):
                    cfile = 'G_I' + str(c.rootname) + '_' + \
                             str(gal_id) + '.in'
                elif(c.repeat == True and c.galcut == True):
                    if ReSize:
                        cfile = 'G_I' + str(gimg)[:-5] + '.in' 
                    else:
                        cfile = 'G_' + str(gimg)[:-5] + '.in'
                else:
                    cfile = 'None'
            if exists(cfile):
                for line_c in open(cfile,'r'): #Reading config file if it exists
                    try:
                        valuec = line_c.split()
                        if(str(valuec[0]) == 'A)'):
                            gimg = str(valuec[1])
                        if(str(valuec[0]) == 'B)'):
                            oimg = (valuec[1])
                        if(str(valuec[0]) == 'C)'):
                            wimg = (valuec[1])
			    if exists(wimg):
				pass
			    else:
				wimg = 'None'
                        if(str(valuec[0]) == 'D)'):
                            pfile = (valuec[1])
                        if(str(valuec[0]) == 'F)'):
                            mimg = (valuec[1])
                        if(str(valuec[0]) == 'G)'):
                            confile = (valuec[1])
                    except:
                        pass
            else:
                cfile = 'None'
            if c.galcut:
                print 'Image is >>> ', gimg
            if(c.galcut == True):
                    ggimg = pyfits.open(gimg)
                    ggimage = ggimg[0].data
                    header0 = ggimg[0].header
                    if (header0.has_key('EXPTIME')):
                        EXPTIME = header0['EXPTIME']
                    else:
                        EXPTIME = -9999
                    if (header0.has_key('RDNOISE')):
                        RDNOISE= header0['RDNOISE']
                    else:
                        RDNOISE = -9999
                    if (header0.has_key('GAIN')):
                        GAIN = header0['GAIN']
                        c.SEx_GAIN = GAIN
                    else:
                        GAIN = -9999
                        c.SEx_GAIN = 1
                    if (header0.has_key('NCOMBINE')):
                        NCOMBINE= header0['NCOMBINE']
                    else:
                        NCOMBINE = -9999
                    if (header0.has_key('FILTER2') or \
                        header0.has_key('FILTER')):
                        try:
                            c.FILTER = header0['FILTER2']
                        except:
                            c.FILTER = header0['FILTER']
                        if c.Filter == 'UNKNOWN':
                            pass
                        else:
                            c.FILTER = c.Filter
                    ggimg.close()
                    SizeXX = ggimage.shape[1]
                    SizeYY = ggimage.shape[0]
            try:
                ximg = float(pdb["ximg"])
            except:
                if(c.galcut == True and RaDecInfo == 0):
                    ximg = SizeXX / 2.0
                else:
                    ximg = -9999
            try:
                yimg = float(pdb["yimg"])
            except:
                if(c.galcut == True and RaDecInfo == 0):
                    yimg = SizeYY / 2.0
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
            except:
                UserGivenPsf = 'None'
            if c.crashhandler and c.starthandle:
                CrashHandlerToRemove(gal_id)
                try:
                    CrashFlag = float(pdb["flag"])
                    CrashFlag = int(CrashFlag)
                except:
                    CrashFlag = 0
                if isset(CrashFlag, GetFlag("GALFIT_FAIL")) or\
                   isset(CrashFlag, GetFlag("BULGE_AT_LIMIT")) or \
                   isset(CrashFlag, GetFlag("DISK_AT_LIMIT")):
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
                if isset(CrashFlag, GetFlag("LARGE_CHISQ")):
                    if isset(CrashFlag, GetFlag("FIT_BULGE_CNTR")) and\
                       isset(CrashFlag, GetFlag("FIT_DISK_CNTR")):
                        pass
                    else:
                        c.fitting[0] = 1
                        c.fitting[1] = 1
                if isset(CrashFlag, GetFlag("FAKE_CNTR")):
                    c.center_deviated = 1
            if(c.galcut == True):   #Given galaxy cutouts
                if exists(sex_cata): #If the user provides sextractor catalogue
                                     #then it will not run SExtractor else do!
                    pass
                else: 
                    RunSex(gimg, wimg, 'None', 9999, 9999, 0)
		    SexShallow(gimg, wimg, 'None', 9999, 9999, 0)
            if(alpha1 == -9999 or delta1 == -9999):
                alpha_j = -9999
                delta_j = -9999
            else:
                alpha_j = (alpha1 + (alpha2 + alpha3 / 60.0) / 60.0) * 15.0
                if delta1 < 0.0:
                    delta_j = delta1 - (delta2 + delta3 / 60.0) / 60.0
                else:
                    delta_j = delta1 + (delta2 + delta3 / 60.0) / 60.0
            for line_s in open(sex_cata,'r'):
                try:
                    values = line_s.split()
                    if(c.galcut == False or RaDecInfo == 1):
                        alpha_s = float(values[3]) #- (c.shiftra) #This is the difference between the observed and the published coordinate for an object. It is used to correct the sextractor cordinate to compare with the published one.
                        delta_s = float(values[4]) # - (c.shiftdec) 
                    elif c.galcut and RaDecInfo == 0:
                        try:
                            alpha_j = float(values[3])
                            delta_j = float(values[4])
                            alpha_s = 9999
                            delta_s = 9999 
                        except:
                            alpha_s = 9999
                            delta_s = 9999 
                    sex_id = values[0]
#                    if(c.galcut == True):
                    xcntr  = float(values[1])
                    ycntr  = float(values[2])
#                    else:
#                        xcntr  = 9999
#                        ycntr  = 9999
                    try:
                        SearchRad = c.searchrad
                    except:
                        SearchRad = '1arc'
                    if SearchRad.endswith('arc'):
                        SeaDeg = float(SearchRad[:-3]) / (60.0 * 60.0)
                        SeaPix = 10.0
                    elif SearchRad.endswith('pix'):
                        SeaPix = float(SearchRad[:-3])
                        if alpha_s == 0.0 and delta_s == 0.0:
                            SeaDeg =  -0.00027
                        else:
                            SeaDeg = 0.00027
                    if(abs(alpha_j - alpha_s) < SeaDeg and \
                       abs(delta_s - delta_j) < SeaDeg or \
                       abs(xcntr - ximg) < SeaPix and \
                       abs(ycntr - yimg) < SeaPix):
                        print "SExtractor ID >>> ", values[0]
                        mag    = float(values[7]) #Magnitude
                        if c.UMag == -500.0:
                            c.UMag = mag - 7.0
                        if c.LMag == 500.0:
                            c.LMag = mag + 7.0
                        halfradius = float(values[9]) #Half light radius
			c.SexHalfRad = float(values[9]) #Sex halfrad to write
                        if c.URe == 500.0:
                            c.URe = halfradius * 10.0
                        if c.URd == 500.0:
                            c.URd = halfradius * 10.0
                        mag_zero = c.mag_zero #magnitude zero point
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
			    sky  = float(values[10]) #sky
                            c.SexSky = float(values[10])
			else:
			    sky  = ShallowSky
			    c.SexSky = ShallowSky
			print sky, float(values[10])
			c.SexMagAuto = float(values[17])
			c.SexMagAutoErr = float(values[18])
                        pos_ang = pa(float(values[11]))
                        axis_rat = 1.0 / float(values[12]) #axis ration b/a
                        eg = 1 - axis_rat
                        ArcR = float(values[11]) * (3.14 / 180.0)
                        if(eg<=0.05):
                            eg = 0.07
                        major_axis = float(values[14])
                        if(alpha_j == -9999 and delta_j == -9999):
                            alpha_j = alpha_s
                            delta_j = delta_s
                        f_err = open('error.log', 'a') 
                        if(c.galcut == True):
                            if ReSize:
                                cutimage = 'I' + gimg
                                whtimage = 'I' + wimg
                            else:
                                cutimage = gimg
                                whtimage = wimg                
                        else:
                            cutimage = 'I' + str(c.rootname) + '_' + \
                                        str(gal_id) + '.fits'
                            whtimage = 'W' + str(c.rootname) + '_' + \
                                        str(gal_id) + '.fits'
                        SizeX = halfradius * FracRad * abs(n.cos(ArcR)) + \
                             axis_rat * halfradius * FracRad * abs(n.sin(ArcR)) 
                        SizeY = halfradius * FracRad * abs(n.sin(ArcR)) + \
                             axis_rat * halfradius * FracRad * abs(n.cos(ArcR))
                        SizeX = int(SizeX)
                        SizeY = int(SizeY)
                        if Square:
                            SizeX = max(SizeX, SizeY)
                            SizeY = max(SizeX, SizeY)
                        if c.galcut:
                            if ReSize:
                                if VarSize:
                                    pass
                                else:
                                    SizeX = FixSize
                                    SizeY = FixSize
                            else:
                                SizeX = SizeXX
                                SizeY = SizeYY
                        else:
                            if VarSize:
                                pass
                            else:
                                SizeX = FixSize
                                SizeY = FixSize        
                        SizeXB = SizeX         #Bookkeeping the size
                        SizeYB = SizeY         #Bookkeeping the size
                        xcntr  = float(values[1])
                        ycntr  = float(values[2])
#                        print 'xcntr, ycntr, SizeX, SizeY', xcntr, ycntr, SizeX, SizeY
                        xmin = int(xcntr) - SizeX 
                        ymin = int(ycntr) - SizeY 
                        xmax = int(xcntr) + SizeX 
                        ymax = int(ycntr) + SizeY
                        xcntrFrac = xcntr - int(xcntr)
                        ycntrFrac = ycntr - int(ycntr) 
                        xminOut = 0
                        yminOut = 0
                        xmaxOut = 0
                        ymaxOut = 0
                        f_err.writelines(['\n\n###########   ', str(gal_id), \
                                          '   ###########\n'])
                        run = 1 #run =1 when pipeline runs sucessfuly
                        c.Flag = 0
                        if c.repeat:
                            c.Flag += 1
                        else:
                            pass
                        if c.fitting[0]:
                            c.Flag += 2
                        else:
                            pass
                        if c.fitting[1]:
                            c.Flag += 4
                        else:
                            pass
                        if c.fitting[2]:
                            c.Flag += 8
                        else:
                            pass
                        try:
                            if(c.repeat == False and c.galcut == False):
                                if(xmin < 0):
                                    xminOut = xmin
                                    xmin = 0
                                if(ymin < 0):
                                    yminOut = ymin
                                    ymin = 0
                                if(xmax > TX):
                                    xmaxOut = xmax
                                    xmax = TX
                                if(ymax > TY):
                                    ymaxOut = ymax
                                    ymax = TY
                                z1 = image[ymin:ymax,xmin:xmax]
                                hdu = pyfits.PrimaryHDU(z1.astype(n.float32))
                                try:
                                    hdu.header.update('RA_TARG', alpha_j)
                                    hdu.header.update('DEC_TARG', delta_j)
                                except:
                                    pass
                                if EXPTIME != -9999:
                                    hdu.header.update('EXPTIME', EXPTIME)
                                else:
                                    pass
                                if RDNOISE != -9999:
                                    hdu.header.update('RDNOISE', RDNOISE)
                                else:
                                    pass
                                if GAIN != -9999:
                                    hdu.header.update('GAIN', GAIN)
                                else:
                                    pass
                                if NCOMBINE != -9999:
                                    hdu.header.update('NCOMBINE', NCOMBINE)
                                else:
                                    pass
                                hdu.writeto(cutimage)
                                
                            if(c.repeat == False and c.galcut and ReSize):
                                fZcuT = pyfits.open(gimg)
                                ZcuT = fZcuT[0].data
                                fZcuT.close()
                                TX = ZcuT.shape[1]
                                TY = ZcuT.shape[0]
                                if(xmin < 0):
                                    xminOut = xmin
                                    xmin = 0
                                if(ymin < 0):
                                    yminOut = ymin
                                    ymin = 0
                                if(xmax > TX):
                                    xmaxOut = xmax
                                    xmax = TX
                                if(ymax > TY):
                                    ymaxOut = ymax
                                    ymax = TY
                                ZcuT1 = ZcuT[ymin:ymax,xmin:xmax]
                                hdu = pyfits.PrimaryHDU(ZcuT1.astype(n.float32))
                                try:
                                    hdu.header.update('RA_TARG', alpha_j)
                                    hdu.header.update('DEC_TARG', delta_j)
                                except:
                                    pass
                                if EXPTIME != -9999:
                                    hdu.header.update('EXPTIME', EXPTIME)
                                else:
                                    pass
                                if RDNOISE != -9999:
                                    hdu.header.update('RDNOISE', RDNOISE)
                                else:
                                    pass
                                if GAIN != -9999:
                                    hdu.header.update('GAIN', GAIN)
                                else:
                                    pass
                                if NCOMBINE != -9999:
                                    hdu.header.update('NCOMBINE', NCOMBINE)
                                else:
                                    pass
                                hdu.writeto(cutimage)
                            try:
                                if(c.repeat == False and c.galcut == False):
                                    if exists(whtfile): 
                                        z2 = weight[ymin:ymax,xmin:xmax]
                                        hdu = pyfits.PrimaryHDU(z2.astype\
                                              (n.float32))
                                        hdu.writeto(whtimage)
                                if(c.repeat == False and c.galcut and ReSize):
                                    if exists(wimg):
                                        fWZcuT = pyfits.open(wimg)
                                        WZcuT = fWZcuT[0].data
                                        fWZcuT.close()
                                        WZcuT1 = WZcuT[ymin:ymax,xmin:xmax]
                                        hdu = pyfits.PrimaryHDU(WZcuT1.astype\
                                                                (n.float32))
                                        hdu.writeto(whtimage)
                                Gal = pyfits.open(cutimage)
                                GalaxyCuT = Gal[0].data
                                Gal.close()
                                GalaxyCuT = n.swapaxes(GalaxyCuT, 0, 1) 
                                SizeX = GalaxyCuT.shape[0]
                                SizeY = GalaxyCuT.shape[1]
                                if c.galcut and ReSize == 0:
                                    pass
                                elif xminOut != 0 or yminOut !=0 or xmaxOut !=0\
                                     or ymaxOut != 0:
                                    c.Flag += 16
                                    if xminOut != 0:
                                        xcntr = SizeXB + xminOut + xcntrFrac
                                    else:
                                        xcntr = SizeX / 2 + xcntrFrac
                                    if yminOut != 0:
                                        ycntr = SizeYB + yminOut 
                                    else:
                                        ycntr = SizeY / 2 + ycntrFrac
                                    if xmaxOut != 0:
                                        xcntr = SizeXB + ycntrFrac 
                                    else:
                                        xcntr = SizeX / 2 + xcntrFrac
                                    if ymaxOut != 0:
                                        ycntr = SizeXB + ycntrFrac
                                    else:
                                        ycntr = SizeY / 2 + ycntrFrac
                                else:
                                    xcntr = SizeX / 2 + xcntrFrac
                                    ycntr = SizeY / 2 + ycntrFrac
#                                print cutimage,xcntr, ycntr, SizeX, SizeY, xminOut, yminOut, xmaxOut, ymaxOut
                                try:
                                #The following function provide the center of blank sky region and the sky sigma    
                                    ElliMaskFunc(cutimage, xcntr, ycntr, \
                                                 SizeX, SizeY, line_s, 0)
                                    try:
                                        Bkgd_Params = BkgdFunc(cutimage, \
                                                xcntr, ycntr, bxcntr, bycntr, \
                                                eg, pos_ang, sky)
                                        bxcntr = Bkgd_Params.bkgd[0]
                                        bycntr = Bkgd_Params.bkgd[1]
                                        skysig = Bkgd_Params.bkgd[2]
                                        print 'Sky Sigma >>> ', skysig
                                    except:
                                        f_err.writelines(['Could not',\
                                                  ' find the sky'\
                                                  ' sigma and mean\n'])

                                except:
                                    f_err.writelines(['Could not create mask ',\
                                                  'for casgm to find the sky'\
                                                  ' sigma and mean. Remove '\
                                                   'if BMask.fits exists\n'])
                                #major axis of the object
                                if(c.decompose):
                                    try:
                                        if(c.repeat == False):
                                            ElliMaskFunc(cutimage, xcntr, \
                                                         ycntr, SizeX, \
                                                         SizeY, line_s, 1)
                                        elif exists('TmpElliMask.fits'):
                                            pass
                                        else:
                                            ElliMaskFunc(cutimage, xcntr, \
                                                         ycntr, SizeX, \
                                                         SizeY, line_s, 1)
                                        ell_mask_file = 'EM_' + \
                                                         str(cutimage)[:-5] + \
                                                        '.fits'
                                        plfile = str(cutimage) + '.pl'
                                        if os.access(plfile, os.F_OK):
                                            os.remove(plfile)
                                        try:
                                            try:
                                                iraf.imcopy(ell_mask_file, \
                                                       plfile, verbose='no')
                                                iraf.flpr()
                                            except:
                                                pass
                                            try:
                                                ell_out = 'E_' + \
                                                     str(cutimage)[:-4] + 'txt'
                                                if os.access(ell_out, os.F_OK):
                                                    os.remove(ell_out)
                                                if os.access('GalEllFit.fits',\
                                                             os.F_OK):
                                                    os.remove('GalEllFit.fits')
                                                run_elli(cutimage, ell_out,\
                                                         xcntr, ycntr, eg, \
                                                      pos_ang, major_axis, sky)
                                                if os.access(plfile, os.F_OK):
                                                    os.remove(plfile)
                                                try:
                                                    iraf.flpr()
                                                except:
                                                    pass
                                            except:
                                                f_err.writelines(['Error '\
                                                           'in ellipse ',\
                                                           'task. Check ',\
                                                           'whether E_',\
                                                           str(cutimage)[:-4],\
                                                           'txt or ellip ',\
                                                           'or err  or ',\
                                                          'test.tab exists\n'])
                                                run = 0
                                                c.Flag += 32
                                        except:
                                            f_err.writelines(['Exists ',\
                                                       str(cutimage),'.pl or ',\
                                                       str(ell_mask_file),\
                                                       ' does not exist\n'])  
                                            run = 0
                                    except:
                                        f_err.writelines(['Error in making '\
                                                     'mask for ellipse task\n'])
                                        run = 0
                            except:
                                f_err.writelines(['The file ', str(whtimage), \
                                                  ' exists\n'])	
                                run = 0
                        except:
                            f_err.writelines(['The file ', str(cutimage),\
                                              ' exists\n'])
                            run = 0
                        if(c.cas):
                            try:
                                ell_mask_file = 'EM_' + \
                                                  str(cutimage)[:-5] + \
                                                  '.fits'
                                if(c.decompose == False):
                                    if c.repeat:
                                        if exists(ell_mask_file):
                                            pass
                                        else:
                                            ElliMaskFunc(cutimage, xcntr, \
                                                         ycntr, SizeX, \
                                                         SizeY, line_s, 1)
                                    else:
                                        ElliMaskFunc(cutimage, xcntr, ycntr,\
                                                     SizeX, SizeY, line_s, 1)
                                try:
                                    caSgm = casgm(cutimage, 'TmpElliMask.fits',\
                                                xcntr, ycntr, bxcntr, bycntr, \
                                                eg, pos_ang, sky, skysig)
                                    C = caSgm[0]
                                    C_err = caSgm[1]
                                    A = caSgm[2]
                                    A_err = caSgm[3]
                                    S = caSgm[4]
                                    S_err = caSgm[5]
                                    G = caSgm[6]
                                    M = caSgm[7]
                                    print 'C, C_err, A, A_err, S, S_err, G,'\
                                    ' M >>> ', str(C)[:5], str(C_err)[:5], \
                                    str(A)[:5], str(A_err)[:5], str(S)[:5], \
                                    str(S_err)[:5], str(G)[:5], str(M)[:5]
                                    if(c.decompose == False):
                                        f_res = open("result.csv", "ab")
                                        writer = csv.writer(f_res)
                                        GalId = str(cutimage)[:-5]
                                        writer.writerow([GalId, alpha_j, \
                                            delta_j, z, c.SexMagAuto, \
					    c.SexMagAutoErr, \
					    C, C_err, A, A_err, S, \
                                            S_err, G, M, c.Flag, c.SexHalfRad])
                                        f_res.close()
                                    f_err.writelines(['(((((CASGM '\
                                                      'Successful)))))'])
                                except:
                                    f_err.writelines(['The CASGM module',\
                                                          ' failed\n'])   
                                    c.Flag += 64
                            except:
                                f_err.writelines(['Could not make mask ',\
                                                      'image for casgm\n'])
                        f_err.close()
                        os.system('rm -f BMask.fits MRotated.fits \
                                  MaskedGalaxy.fits Rotated.fits')
                        if(c.decompose == False):
                            if os.access(ell_mask_file, os.F_OK):
                                os.remove(ell_mask_file)
                        f_err = open('error.log', 'a') 
                        if(c.decompose):
                            for psfelement in psflist:
                                psfradec(psfelement)
                            try:
                                if(c.repeat == False and cfile == 'None'):
                                    if(alpha_s == 9999 or delta_s == 9999):
                                        if UserGivenPsf == 'None':
                                            psffile = c.psflist[psfcounter]
                                        else:
                                            psffile = UserGivenPsf
                                            psfradec(psffile)
                                        try:
                                            p=pyfits.open(psffile)
                                            header = p[0].header
                                            if(header.has_key('RA_TARG')):
                                                ra_p = header['RA_TARG']
                                            if (header.has_key('DEC_TARG')):
                                                dec_p = header['DEC_TARG']
                                            p.close()
                                            r = 3.14159265 / 180.0
                                            distance = 3600.0*n.sqrt((delta_j\
                                                       - dec_p)**2.0 + \
                                                       ((alpha_j - ra_p) * \
                                                       n.cos(delta_j * r))**2.0)
                                        except:
                                            distance = 9999
                                        psfcounter += 1
                                    else:
                                        psffile, distance = \
                                               psf_select(alpha_j, delta_j)
                                        distance = distance * 60.0 * 60.0
                                else:
                                    if(alpha_s == 9999 or delta_s == 9999):
                                        distance = 9999
                                    else:
                                        psfradec(pfile)
                                        p=pyfits.open(pfile)
                                        header = p[0].header
                                        if(header.has_key('RA_TARG')):
                                            ra_p = header['RA_TARG']
                                        else:
                                            ra_p = 9999
                                        if (header.has_key('DEC_TARG')):
                                            dec_p = header['DEC_TARG']
                                        else:
                                            dec_p = 9999
                                        p.close()
                                        r = 3.14159265 / 180.0
                                        if(ra_p == 9999 or dec_p == 9999):
                                            distance = 9999
                                        else:
#                                            distance = n.sqrt((delta_j - \
#                                            dec)**2.0 + ((alpha_j - ra) * \
#                                            n.sin((0.5) *\
#                                            (delta_j + dec)))**2.0)
                                            distance = 3600.0*n.sqrt((delta_j\
                                                       - dec_p)**2.0 + \
                                                       ((alpha_j - ra_p) * \
                                                       n.cos(delta_j * r))**2.0)
#                                            distance = n.arccos(n.cos((90.0 - \
 #                                               delta_j) \
 #                                            * r) * n.cos((90.0 - dec_p) * r) \
  #                                           + n.sin((90.0 - delta_j) * r) *  \
  #                                           n.sin((90.0 - dec_p) * r) * \
  #                                           n.cos((alpha_j - ra_p) * r))
                                            #print 'alp dec alpsf decpsf d', alpha_j, delta_j, ra_p, dec_p, distance
                                if(cfile == 'None'):
                                    if c.manual_mask:
                                        ManualMaskManager(cutimage)
                                    else:
                                        MaskFunc(cutimage, xcntr, ycntr, \
                                                 SizeX, SizeY, line_s)
                                    maskimage = 'M_' + str(cutimage)[:-5] +\
                                                '.fits'
                                else:
                                    maskimage = mimg
                                    if exists(mimg):
                                        pass
                                    else:
                                        if c.manual_mask:
                                            ManualMask(cutimage)
                                        else:
                                            MaskFunc(cutimage, xcntr, ycntr, \
                                                     SizeX, SizeY, line_s)
                                try:
                                    if(cfile == 'None'):
                                        ConfigFunc(cutimage, whtimage,  xcntr,\
                                                   ycntr, SizeX, \
                                                   SizeY, line_s, psffile)
                                        config_file = 'G_' + \
                                                       str(cutimage)[:-5]+ '.in'
                                        outimage = 'O_' + str(cutimage)
                                    else:
                                        config_file = cfile
                                        outimage = str(oimg)
                                    f_fit = open('fit2.log','a')
                                    if exists('fit.log'):
                                        os.system('rm fit.log')
                                #Here the user should tell the location of the GALFIT excutable                     
                                    try:
                                        DetailFit = c.detail
                                    except:
                                        DetailFit = 0
                                    if c.galfit and DetailFit:
                                        ConfigIter(cutimage, whtimage,  xcntr,\
                                                   ycntr, SizeX, \
                                                   SizeY, line_s, psffile)
                                    elif c.galfit:
                                        cmd = str(c.GALFIT_PATH) + ' ' + \
                                                  config_file
                                        os.system(cmd)

#                                        os.system('/Vstr/vstr/vvinuv/galfit/modified/galfit "' + config_file + '"')
                                    if exists('fit.log'):
                                        for line in open('fit.log','r'):
                                            f_fit.writelines([str(line)])
                                    f_fit.close()
                                    try:
                                        if(c.repeat == False):
                                            OutMaskFunc(outimage, xcntr, \
                                                        ycntr,  SizeX, \
                                                        SizeY, line_s)
                                        out_mask_file = 'OEM_' + \
                                                         str(outimage)[:-5] + \
                                                        '.fits'
                                        outplfile = 'S' + str(outimage) + '.pl'
                                        if os.access(outplfile, os.F_OK):
                                            os.remove(outplfile)
                                        try:
                                            ell_output = 'OE_' + \
                                                   str(cutimage)[:-4] + 'txt'
                                            outmodel = 'S' + outimage
                                            try:
                                                iraf.imcopy(out_mask_file, \
                                                       outplfile, verbose='no')
                                                iraf.flpr()
                                                iraf.imcopy(outimage + '[2]', \
                                                       outmodel, verbose='no')
                                                iraf.flpr()
                                            except:
                                                pass
                                            if os.access(ell_output, \
                                                         os.F_OK):
                                                os.remove(ell_output)  
                                            try:
                                                FMo=pyfits.open(outimage)
                                                MoDel = f[2].data
                                                FMo.close()
                                                MoDel = n.swapaxes(MoDel, \
                                                        0, 1)
                                                MoShapX = MoDel.shape[0] /2
                                                MoShapY = MoDel.shape[1] /2
                                                MoCen = center_of_mass( \
                                                MoDel[MoShapX-5:MoShapX+5, \
                                                      MoShapY-5:MoShapY+5])
                                                MoX = MoShapX + MoCen[0] -5
                                                MoY = MoShapY + MoCen[1] -5
                                            except:
                                                MoX = xcntr
                                                MoY = ycntr
                                            try:
                                                if os.access('GalEllFit.fits',\
                                                          os.F_OK):
                                                    os.remove('GalEllFit.fits')
                                                run_elli(outmodel, ell_output,\
                                                     MoX, MoY, eg, \
                                                    pos_ang, major_axis, sky)
                                                try:
                                                    iraf.flpr()
                                                except:
                                                    pass
                                            except:
                                                f_err.writelines(['Error in '\
                                                          'ellipse '\
                                                          'task. Check ', \
                                                          'whether ' ,\
                                                           str(ell_output) ,\
                                                      ' or ellip or err  or',\
                                                      ' test.tab exists\n'])                                               
                                            for myfile in [outplfile, \
                                                           outmodel]:
                                                if os.access(myfile, \
                                                             os.F_OK):
                                                    os.remove(myfile) 
                                        except:
                                            f_err.writelines(['Exists ',\
                                                       str(outimage),'.pl or ',\
                                                       str(out_mask_file),\
                                                       ' does not exist.\n'])  
                                            f_err.writelines(['GALFIT '\
					              'MIGHT BE CRASHED\n'])
                                            c.Flag += 128
                                            failedgalfit(cutimage)
                                            run = 0
                                    except:
                                        f_err.writelines(['Error in making '\
                                                'out mask for ellipse task\n'])
                                        run = 0 
                                except:
                                    f_err.writelines(['Error in writing',\
                                                      ' configuration file\n'])	
                                    run = 0
                            except:
                                f_err.writelines(['Error in making mask for '\
                                                  'galfit\n'])
                                run = 0
#                        if exists('plot_' + str(cutimage)[6:-4] + 'png'):	
#                            os.system('rm ''plot_' + str(cutimage)[6:-4] + 'png''')
                            if(run == 1 or run == 0):
                                try:
                                    if exists('P_' + str(cutimage)[6:-4] \
                                              + 'png'):	
                                        os.system('rm ''P_' + str(cutimage)\
                                                   [6:-4] + 'png''')
                                    GoodNess = PlotFunc(cutimage, outimage, \
                                          maskimage, xcntr, ycntr, sky, skysig)
                                    Goodness = GoodNess.plot_profile
                                except:
                                    f_err.writelines(['Error in plotting. '])
                                    if(maskimage == 'None'):
                                        f_err.writelines(['Could not find '\
                                                          'Mask image\n'])
                                    run = 0	
                                    Goodness = 9999
                                    c.Flag += 256
                                try:
                                    EXPTIME = EXPTIME * 1.0
                                except:
                                    EXPTIME = 9999
                                try:
                                    write_params(cutimage, xcntr, ycntr, \
                                                 distance, alpha_j, \
                                                 delta_j, z, Goodness, \
                                                 C, C_err, A, A_err, S, S_err, \
                                                 G, M, EXPTIME)
#                                f_err.writelines(['(((((((((( Successful', \
 #                                                     ' ))))))))))\n'])
                                except:
                                    try:
                                        write_params(cutimage, xcntr, ycntr, \
                                                     distance, alpha_j,\
                                                     delta_j, z, Goodness, \
                                                     9999, 9999, 9999,\
                                                     9999, 9999, 9999, 9999, \
                                                     9999, EXPTIME)
                                    except:
                                        f_err.writelines(['Error in writing '\
                                                          'html\n'])
                                        run = 0
                            if(run == 1):
                                f_err.writelines(['((((( Decomposition '\
                                                  'Successful )))))\n'])
						
#iraf.imcopy(str(imagefile) + '[' + str(xmin) + ':' + str(xmax) + ',' + str(ymin) + ':' + str(ymax) + ']', cutimage)	
#iraf.imcopy(str(whtfile) + '[' + str(xmin) + ':' + str(xmax) + ',' + str(ymin) + ':' + str(ymax) + ']', whtimage)	
#					fitellifunc(gal_id, line_s)
                            if isset(c.Flag, GetFlag("GALFIT_FAIL")) or \
                               isset(c.Flag, GetFlag("LARGE_CHISQ")) or \
                               isset(c.Flag, GetFlag("FAKE_CNTR")) or \
                               isset(c.Flag, GetFlag("BULGE_AT_LIMIT")) or \
                               isset(c.Flag, GetFlag("DISK_AT_LIMIT")):
                                FailedValues = line_j.split()
                                for FailedValue in FailedValues:
                                    f_failed.writelines([str(FailedValue), ' '])
                                f_failed.writelines([str(c.Flag), '\n'])
                            f_err.close()
                            f_cat.writelines([str(gal_id), ' '])
                            f_cat.write(line_s)
                            for myfile in ['ellip','err','test.tab']:
                                if os.access(myfile,os.F_OK):
                                    os.remove(myfile)
                except:
                    pass
            if(c.galcut == True):
                if os.access(sex_cata, os.F_OK):
                    os.remove(sex_cata)
        except:
            pass
    f_cat.close()
    f_failed.close()
def selectpsf(ImG, CaT):
    c.psff = []
    im = pyfits.open(ImG)
    image = im[0].data
    im.close()
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
                    PsfSize = n.floor(float(values[14])) * c.starsize 
                    x1 = int(xcntr) + (PsfSize/2)
                    x2 = int(xcntr) - (PsfSize/2)
                    y1 = int(ycntr) + (PsfSize/2)
                    y2 = int(ycntr) - (PsfSize/2)
                    ra1 = int(float(values[3]) / 15.0)
                    ra2 = int((float(values[3]) / 15.0 - int(float(values[3]) / \
                          15.0))*60.0)
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
                        ra33 = '0' + (str(n.round(ra3, 1))[:3]).split('.')[0] + \
                                     (str(n.round(ra3, 1))[:3]).split('.')[1] 
                    else:
                        ra33 = (str(n.round(ra3, 1))[:4]).split('.')[0] + \
                               (str(n.round(ra3, 1))[:4]).split('.')[1]
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
                        dec33 = '0' + (str(n.round(dec3, 1))[:3]).split('.')[0]\
                                 + (str(n.round(dec3, 1))[:3]).split('.')[1]
                    else:
                        dec33 = (str(n.round(dec3, 1))[:4]).split('.')[0] +\
                                (str(n.round(dec3, 1))[:4]).split('.')[1]
                    psffile = 'psf_' + str(ra11) + str(ra22) + str(ra33) + str(dec11) +str(dec22) + str(dec33) + '.fits'
                    if psffile in c.psff:
                        pass
                    else:
                        c.psff.append(psffile)
                        psf = image[y2:y1, x2:x1]
                        psf = psf - BaKgR
                        if os.access(psffile, os.F_OK):
                            os.remove(psffile)
                        hdu = pyfits.PrimaryHDU(psf.astype(n.float32))
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
        c.SEx_SEEING_FWHM = 0.11
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
def UsageOfPyMorph():
    print "Usage: pymorph [--edit-conf[-e]] [--with-psf] [--force[-f]] "\
          "[--help[-h]] [--lmag] [--umag] [--lu] [--un] [--lre] [--ure] "\
	  "[--lrd] [--urd] [--with-in] [--with-filter] [--with-db] "\
	  "[--with-area]  [--no-mask] [--norm-mask] [--with-sg]"
    sys.exit(0)

if __name__ == '__main__':
    c.FirstCreateDB = 1 #Won't create table when c.FirstCreateDB=0
    c.VERSION = 1.7
    c.FILTER = 'UNKNOWN'
    c.Filter = 'UNKNOWN'
    try:
        if(c.repeat == False and c.galcut == False):
            img = pyfits.open(c.imagefile)
            c.ImAgE = img[0].data
            c.HeAdEr0 = img[0].header
            if (c.HeAdEr0.has_key('GAIN')):
                c.SEx_GAIN = c.HeAdEr0['GAIN']
            else:
                c.SEx_GAIN = 1
            if (c.HeAdEr0.has_key('FILTER2') or c.HeAdEr0.has_key('FILTER')):
                try:
                    c.FILTER = c.HeAdEr0['FILTER2']
                except:
                    c.FILTER = c.HeAdEr0['FILTER']
            img.close()
    except IOError, (errno, strerror):
        print imagefile, "I/O error(%s): %s" % (errno, strerror)
        os._exit(0)
    c.SEx_DETECT_MINAREA = 6
    c.SEx_DETECT_THRESH = 1.5
    c.SEx_ANALYSIS_THRESH = 1.5
    c.SEx_FILTER = 'Y'
    c.SEx_FILTER_NAME = 'default.conv'
    c.SEx_DEBLEND_NTHRESH = 32
    c.SEx_DEBLEND_MINCONT = 0.005
    c.SEx_PHOT_FLUXFRAC = 0.5
    c.SEx_PIXEL_SCALE = c.pixelscale
    c.SEx_SEEING_FWHM = 0.11
    c.SEx_BACK_SIZE = 64
    c.SEx_BACK_FILTERSIZE = 3
    c.SEx_BACKPHOTO_TYPE = 'GLOBAL'
    c.SEx_BACKPHOTO_THICK = 24
    c.SEx_WEIGHT_TYPE = 'DECIDE'
    c.WhichPsf = 0
    c.LMag = 500.0
    c.UMag = -500.0
    c.LN = 0.1
    c.UN = 20.0
    c.LRe = 0.0
    c.URe = 500.0
    c.LRd = 0.0
    c.URd = 500.0
    c.avoideme = 150.0
    c.AreaOfObj = 40.0 #Area of psf for selection
    c.StarGalProb = 0.8
    c.bdbox = 0
    c.bbox = 0
    c.dbox = 0
    c.NoMask = 0
    c.NormMask = 0
    sex_cata = c.sex_cata
    if len(sys.argv[1:]) > 0:
        try:
            options, args = getopt(sys.argv[1:], "efhti", ['edit-conf', \
                        'with-psf=', 'force', 'help', 'test', 'initial',\
                        'lmag=', 'umag=', 'ln=', 'un=', 'lre=', 'ure=', \
                        'lrd=', 'urd=', 'with-in=', 'with-filter=', \
                        'with-db=', 'with-area=', 'no-mask', 'norm-mask', \
			'with-sg='])
        except GetoptError, err:
            print str(err) 
            UsageOfPyMorph()
        for opt, arg in options:
            if opt in ('-c', '--edit-conf'):
                SExtractorConf()
            if opt in ('-i', '--initial'):
                DecideInitialParamers
            if opt in ('-f', '--force'):
                if exists(sex_cata):
                    os.remove(sex_cata)
                else:
                    pass
            if opt in ('-p', '--with-psf'):
                c.WhichPsf = int(arg)               #Nearest/farthest psf
            if opt in ('-t', '--test'):
                TestingOption
            if opt in ['--lmag']:
                c.LMag = float(arg)
            if opt in ['--umag']:
                c.UMag = float(arg)
            if opt in ['--ln']:
                c.LN = float(arg)
            if opt in ['--un']:
                c.UN = float(arg)
            if opt in ['--lre']:
                c.LRe = float(arg)
            if opt in ['--ure']:
                c.URe = float(arg)
                print c.URe
            if opt in ['--lrd']:
                c.LRd = float(arg)
            if opt in ['--ure']:
                c.URd = float(arg)
            if opt in ['--with-in']:
                c.avoideme = float(arg)
            if opt in ['--with-filter']:
                c.Filter = arg
            if opt in ['--with-db']:
                c.database = arg
            if opt in ['--with-area']:
                c.AreaOfObj = float(arg)
            if opt in ['--no-mask']:
		c.NoMask = 1
            if opt in ['--norm-mask']:
		c.NormMask = 1
            if opt in ['--with-sg']:
                c.StarGalProb = float(arg)
            if opt in ['--bdbox']:
                c.bdbox = 1
            if opt in ['--bbox']:
                c.bbox = 1
            if opt in ['--dbox']:
                c.dbox = 1
            if opt in ('-h', '--help'):
                UsageOfPyMorph()
    if c.Filter == 'UNKNOWN':
        pass
    else:
        c.FILTER = c.Filter
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
    if exists(sex_cata):
        pass
    elif(c.galcut == False):
        print 'The SExtractor catalogue for your frame is NOT found. ' \
              'One is making using the default values. It is always '\
              'recommended to make SExtractor catalogue by YOURSELF as '\
              'the pipeline keep the sky value at the SExtractor value '\
              'during the decomposition.'
        if exists(c.whtfile):
            RunSex(c.imagefile, c.whtfile, 'None', 9999, 9999, 0)
	    SexShallow(c.imagefile, c.whtfile, 'None', 9999, 9999, 0)
        else:
            RunSex(c.imagefile, 'None', 'None', 9999, 9999, 0)
	    SexShallow(c.imagefile, 'None', 'None', 9999, 9999, 0)
    def runpsfselect():
        if(c.galcut):   #Given galaxy cutouts
            obj_file = open(c.clus_cata,'r') 
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
                        if exists('I' + str(c.rootname) + '_' \
                                   + str(gal_id) + '.fits'):
                            gimg = 'I' + str(c.rootname) + '_' \
                                    + str(gal_id) + '.fits'
                        elif exists(str(gal_id) + '.fits'):
                            gimg = str(gal_id) + '.fits'
                        else:
                            print "No image found. Exiting"
                            os._exit(0)
                    try:
                        wimg = pdb["wimg"]   #Weight cut
                    except:
                        if exists('W' + str(c.rootname) + '_' + \
                                  str(gal_id) + '.fits'):
                            wimg = 'W' + str(c.rootname) + '_' + \
                                   str(gal_id) + '.fits'
                        else:
                            wimg = 'None'
                    GiMg = pyfit.open(gimg)
                    headerGiMg = GiMg[0].header
                    if (headerGiMg.has_key('GAIN')):
                        c.SEx_GAIN = headerGiMg['GAIN']
                    else:
                        c.SEx_GAIN = 1
                    GiMg.close()
                    if exists(sex_cata): 
                        pass
                    else:
                        RunSex(gimg, wimg, 'None', 9999, 9999, 0)
                    try:
                        selectpsf(gimg, sex_cata)
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
                cmd = 'mv CatWithPsf.cat ' + str(c.clus_cata)
                os.system(cmd)
            else:
                pass
        else:
            selectpsf(c.imagefile, sex_cata)
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
        c.psflist = '@psflist.list'
        FindAndFit()
        main()
        if c.crashhandler:
            c.starthandle = 1
            os.system('mv restart.cat CRASH.CAT')
            c.clus_cata = 'CRASH.CAT' 
            main()

    elif c.psfselect == 1:
        c.Interactive = 1
        os.system('ds9 &')
        time.sleep(2)
        runpsfselect()
        os.system('xpaset -p ds9 quit')
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
