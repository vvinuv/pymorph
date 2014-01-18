import numpy as np
import numpy.ma as ma
import pyfits
import types
import os
from os.path import exists
import config as c
import csv

#from ellimaskfunc_easy import ElliMaskFunc
from ellimaskfunc import ElliMaskFunc
from bkgdfunc import BkgdFunc
from casgm import casgm
from flagfunc import GetFlag, isset, SetFlag
#from outmaskfunc_easy import OutMaskFunc
from outmaskfunc import OutMaskFunc


def PymorphError(inst, msg):
    '''Some way to handle error'''
    print msg
    if c.noerror:
        pass
    else:
        raise inst

def GetSizeInfo():
    ReSize = c.size[0]
    VarSize = c.size[1]
    Square = c.size[3]
    FracRad = c.size[2]
    if VarSize == 0:
        FixSize = c.size[4]
    else:
        FixSize = 120
    return ReSize, VarSize, Square, FracRad, FixSize

def GetWhtImage():
    print "Using large image. c.imagefile >>> ", imagefile
    TY, TX = c.imagedata.shape
    if exists(os.path.join(c.datadir, c.whtfile)):
        print "Weight file >>> ", c.whtfile
        wht = pyfits.open(os.path.join(c.datadir, c.whtfile))
        if re.search("rms", c.whtfile.lower()):
            c.weightdata = wht[0].data
            c.weightexists = 1
        elif re.search("weight", c.whtfile.lower()):
            c.weightdata = 1 / np.sqrt(wht[0].data)
            c.weightexists = 1
        else:
            print 'Weight file is not understood. Please include ' + \
                  'the word weight/rms to the weight file name. ' + \
                  'If it is weight the rms will be found by 1/sqrt(w)'
        wht.close()
    else:
       print 'No weight image found\n'
       c.weightexists = 0
    if c.weightexist:
        WY, WX = c.weightdata.shape
        if WY != TY or WX != TX:
            raise ValueError('Dimensions of image and weight images \
                              does not match')
           
def GetModel():
    try:
        ComP = c.components
    except:
        print "c.components undefined. Asuming bulge+disk model"
        ComP = ['bulge', 'disk']
    if len(ComP) == 0:
        print "No model specified. Asuming bulge+disk model"
        ComP = ['bulge', 'disk']
    return ComP

def GetInputParams(pnames, line):
    pdb = {} 
    values = line.split()
    for pname, value in zip(pnames, values):
        pdb[pname] = value
     
    if "gal_id" in pnames:
        gal_id = pdb["gal_id"]
    else:
        print "no gal_id using gimg"
        if "gimg" in pnames:
            gal_id = pdb["gimg"].split('.')[0] #gal_id will be 
                                               #filename without .fits
        else:
            raise ValueError("No image or gal_id found in the input\
                              catalogue. Exiting")
    c.fstring = '%s_%s'%(c.rootname, gal_id)
    if 'ra1' in pnames:
        alpha1 = float(pdb["ra1"])
    else:
        print "No ra1 (hour) is given"
        alpha1 = -9999
    if 'ra2' in pnames:
        alpha2 = float(pdb["ra2"])
    else:
        print "No ra2 (minuite) is given. Asuming ra is in deg"
        alpha2 = 0
    if 'ra3' in pnames:
        alpha3 = float(pdb["ra3"])
    else:
        print "No ra3 (second) is given"
        alpha3 = 0
    if 'dec1' in pnames:
        delta1 = float(pdb["dec1"])
    else:
        print "No dec1 (deg) is given"
        delta1 = -9999
    if 'dec2' in pnames:
        delta2 = float(pdb["dec2"])
    else:
        print "No dec2 (min) is given"
        delta2 = 0
    if 'dec3' in pnames:
        delta3 = float(pdb["dec3"])
    else:
        print "No dec3 (sec) is given"
        delta3 = 0
    if 'ra1' in pnames or 'dec1' in pnames:
        RaDecInfo = 1 #Understood position is given
    else:
        RaDecInfo = 0 #Understood position is given
    if 'z' in pnames:
        z = float(pdb["z"])
    else:
        print "No z is given"
        z = 9999
    if 'bxcntr' in pnames:
        bxcntr = float(pdb["bxcntr"])
    else:
        bxcntr = 9999
    if 'bycntr' in pnames:
        bycntr = float(pdb["bycntr"])
    else:
        bycntr = 9999
    if 'mzero' in pnames:
        c.mag_zero = float(pdb["mzero"])
    else:
        c.mag_zero = c.mag_zero
    if 'star' in pnames:
        UserGivenPsf = pdb["star"]
        print 'Psf is assigned individually to galaxies'
    else:
        UserGivenPsf = 'None'
    if 'sky' in pnames:
        UserGivenSky = pdf['sky']
        print 'Sky is assigned individually to galaxies'
    else:
        UserGivenSky = -9999

    if 'gimg' in pnames:
        gimg = pdb["gimg"]    #Galaxy cutout
    else:
        gimg = None
    if 'wimg' in pnames:
        wimg = pdb["wimg"]   #Weight cut
    else:
        wimg = None

    if 'ximg' in pnames:
        ximg = float(pdb["ximg"])
        if ReSize and c.galcut and c.repeat:
            ximg = TX / 2.0
    else:
        print 'No ximg is given. Trying to find from the cutout if '\
              'no ra dec in the image header'
        if(c.galcut == True and RaDecInfo == 0):
            ximg = TX / 2.0
        else:
            ximg = -9999
    if 'yimg' in pnames:
        yimg = float(pdb["yimg"])
        if ReSize and c.galcut and c.repeat:
            yimg = TY / 2.0
    else:
        print 'No yimg is given. Trying to find from the cutout if '\
              'no ra dec in the image header'
        if(c.galcut == True and RaDecInfo == 0):
            yimg = TY / 2.0
        else:
            yimg = -9999

    return alpha1, alpha2, alpha3, delta1, delta2, delta3, z, bxcntr, bycntr, UserGivenPsf, UserGivenSky, gimg, wimg, ximg, yimg, RaDecInfo


def DecisionMaker(gimg, wimg): 
    if c.repeat:
        cfile = 'G_%s.in'%c.fstring
        gimg, oimg, wimg, pfile, mimg, confile = ReadGalfitConfig(cfile)
        cutimage = gimg
    else:
	cfile = 'G_%s.in'%c.fstring
	confile = '%s.con'%c.fstring
	mimg = 'M_%s.fits'%c.fstring
	oimg = 'O_%s.fits'%c.fstring
        pfile = None
        if c.galcut:
            if gimg is not None:
		if not exist(gimg) and not noerror:
		    raise ValueError('c.galcut=True, but no cutout image is found')
		elif not exist(gimg) and noerror:
		    print 'Not found %s. Skipping this galaxy'%gimg
		else exist(gimg):
		    if ReSize:
			cutimage = 'I%s.fits'%c.fstring
			whtimage = 'W%s.fits'%c.fstring
		    else:
			cutimage = gimg
			whtimage = wimg
            else:
		if not exist(gimg) and not noerror:
		    raise ValueError('c.galcut=True, but no cutout image is found')
		elif not exist(gimg) and noerror:
		    print 'Not found %s. Skipping this galaxy'%gimg
        else:
	    cutimage = 'I%s.fits'%c.fstring
	    whtimage = 'W%s.fits'%c.fstring
    return  gimg, oimg, wimg, cfile, pfile, mimg, confile

def SetSkies(id, deep_sky):
    # The shallow sky from the shallow run. If no shallow
    # sky it uses the deep sky
    ShallowSky = 9999
    if exists('%s.Shallow'%c.sex_cata):
	f_sex_shallow = open('%s.Shallow'%c.sex_cata, 'r')
	for line_shallow in f_sex_shallow:
	    v_shallow = line_shallow.split()
	    try:
		if str(v_shallow[19]) == str(id):
		    ShallowSky = float(v_shallow[10])
	    except:
		pass
	f_sex_shallow.close()
    if ShallowSky == 9999:
	c.SexSky = float(deep_sky)
	c.GalSky = 9999
    else:
	c.SexSky = ShallowSky
	c.GalSky = 9999

def SetCrashHandler():
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


def FindSexObj(sex_cata, RaDecInfo, SeaDeg, SeaPix, alpha_j, delta_j, ximg, yimg):
    # first count the number of "potential" targets in the search radius
    c.SexTargets = 0
    good_object = ''

    os.system('cp %s sex_%s.txt' %(sex_cata, c.fstring))
    new_distance = 999.0 #the distance from the center to the best target
    for line_s in open(sex_cata, 'r'):
        try:
            values = line_s.split()
            sex_id = values[0]
            alpha_s = float(values[3])
            delta_s = float(values[4])
            xcntr  = float(values[1])
            ycntr  = float(values[2])

            if RaDecInfo:
                curr_distance = np.sqrt((alpha_j - alpha_s)**2+(delta_s - delta_j)**2)
                if curr_distance < SeaDeg:
                    print "Candidate distance: %.3f" %curr_distance
                    c.SexTargets +=1
                    if curr_distance < new_distance:
                        new_distance = curr_distance
                        print "New Preferred target!!"
                        good_object = line_s
     
            else:
                curr_distance = np.sqrt((xcntr - ximg)**2+(ycntr - yimg)**2)

                if curr_distance < SeaPix:
                    print "Candidate distance: %.3f" %curr_distance    
                    c.SexTargets +=1
                    if curr_distance < new_distance:
                        new_distance = curr_distance
                        print "New Preferred target!!"
                        good_object = line_s
        except Exception as inst:
            if values[0].strip().isdigit():
                PymorphError(inst, 'Something happend in the pipeline. Check error.log')
            else:
                pass
    print "Target distance: %.3f" %new_distance
    return good_object

def Gamma(z):
    """This is the Lanczos approximation for Gamma function"""
    lanczosG = 7
    lanczos_coef = [0.99999999999980993, 676.5203681218851,\
                    -1259.1392167224028, 771.32342877765313,\
                    -176.61502916214059, 12.507343278686905, \
                    -0.13857109526572012, 9.9843695780195716e-6,\
                    1.5056327351493116e-7]
    if z < 0.5:
        answer =  n.pi / (n.sin(n.pi*z)*Gamma(1-z))
    else:
        z -= 1
        x = lanczos_coef[0]
        for i in range(1, lanczosG + 2):
            x += lanczos_coef[i]/(z + i)
        t = z + lanczosG + 0.5
        answer =  n.sqrt(2*n.pi) * t**(z + 0.5) * n.exp(-t) * x
    return answer

def bfunc(x):
        """ This function gives value of b_n given the Sersic index"""
        return 0.868242*x -0.142058 # Khosroshahi et al. 2000 approximation
    
def Get_R():
    """Return 3.14159265 / 180.0"""
    return np.pi / 180.0

def RaDegToHMS(ra):
    """Convert to Ra deg to H:M:S"""
    ra = ra / 15.0
    ra1 = int(ra)
    ra22 = (ra - ra1) * 60.0
    ra2 = int(ra22)
    ra3 = (ra22 - ra2) * 60.0
    return ra1, ra2, ra3

def DecDegToDMS(dec):
    """Convert to Dec deg to D:M:S"""
    dec1 = int(dec)
    dec22 = abs(dec - dec1) * 60.0
    dec2 = int(dec22)
    dec3 = (dec22 - dec2) * 60.0
    return dec1, dec2, dec3

def HMSToDeg(h, m, s):
    """Convert H:M:S to deg"""
    ra = (h + (m + s / 60.0) / 60.0) * 15.0
    return ra

def DMSToDeg(d, m, s):
    """Convert D:M:S to deg"""
    if d < 0.0:
        dec = (d - (m + s / 60.0) / 60.0)
    else:
        dec = (d + (m + s / 60.0) / 60.0)
    return dec

def PsfArr():
    """Return psf list if the given input is a file"""
    if type(c.psflist) is types.StringType:
        psffi = open(os.path.join(c.datadir, c.psflist.split('@')[1]), 'r')
        c.psflist = []
        for pline in psffi:
            c.psflist.append(pline.split()[0])
    elif type(c.psflist) is types.ListType:
        pass
    else:
        print "The psf list is not understood. Please use either \
               @filename or a list of psfs"


def UpdatePsfRaDec(element):
    """The function which will update the psf header if the psf files
       are in the specified format"""
    try:
        ra1 = float(str(element)[4:6])
        ra2 = float(str(element)[6:8])
        ra3 = float(str(element)[8:10]) + float(str(element)[10]) / 10.0
        dec1 = float(str(element)[11:-10])
        dec2 = float(str(element)[-10:-8])
        dec3 = float(str(element)[-8:-6]) + float(str(element)[-6]) / 10.0
        ra = HMSToDeg(ra1, ra2, ra3)
        dec = DMSToDeg(dec1, dec2, dec3)
        data, header = pyfits.getdata(element, header=True)
        header.update('RA_TARG', ra, "RA")
        header.update('DEC_TARG', dec, "DEC")
        pyfits.update(element, data, header)
    except Exception as inst:
        PymorphError(inst, 'Check PSF name is not in the format')


def FailedGalfit(WhichGalaxy):
    """This function creates fake fit.log if galfit fails"""
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
    if(float(x) >= 0 and float(x) <= 180.0):
        pos_ang = float(x) - 90.0 #position angle
    if(float(x) < 0 and float(x) >= -180.0):
        pos_ang = 90.0 - abs(float(x))  #position angle
    if(float(x) > 180 and float(x) <= 360.0):
        pos_ang = float(x) - 360.0 + 90.0 #position angle
    if(float(x) >= -360 and float(x) < -180.0):
        pos_ang = float(x) + 360.0 - 90.0 #position angle   
    return pos_ang


def SelectPsf(alpha_j, delta_j):                                   
    """This function will select the nearest psf from the psflist.
       The distance is calculated by using the following equation
       d = Sqrt((dec_a - dec_b) ^ 2 + ((ra_a - ra_b) * sin(0.5) * 
       (dec_a - dec_b)) ^ 2.0 )"""
    PsfDistanceDict = {}
    psffile = 'test.fits'
    r = Get_R()
    for element in c.psflist:
        p = pyfits.open(os.path.join(c.datadir, element))
        header = p[0].header
        if (header.has_key('RA_TARG')):
            ra = header['RA_TARG']
        elif (header.has_key('RA')):
            ra = header['RA']
        else:
            ra = 9999
        if (header.has_key('DEC_TARG')):
            dec= header['DEC_TARG']
        elif (header.has_key('DEC')):
            dec= header['DEC']
        else:
            dec= 9999
        p.close()
#       d = sqrt((ra - alpha_j) ** 2.0 + (dec - delta_j) ** 2.0)
#       d = np.arccos(np.cos((90.0 - delta_j) * r) * np.cos((90.0 - dec) *\
#           r) + np.sin((90.0 - delta_j) * r) *  np.sin((90.0 - dec) * r) * \
#           np.cos((alpha_j - ra) * r))
#       d = np.sqrt((delta_j - dec)**2.0 + ((alpha_j-ra)*np.sin((0.5) *\
#           (delta_j+dec)))**2.0)
        d = np.sqrt((delta_j - dec)**2.0 + ((alpha_j - ra) * \
            np.cos(delta_j * r))**2.0)
        PsfDistanceDict[element] = d
    ItemS = PsfDistanceDict.items()
    ItemS = [(v, k) for (k, v) in ItemS]
    ItemS.sort()
    ItemS = [(k, v) for (v, k) in ItemS]
    psffile = ItemS[c.WhichPsf][0]
    distance = ItemS[c.WhichPsf][1]
    return psffile, distance 

def PyMorphFiles(gal_id):
    """Returns all the temporary files produced by pymorph"""
    files = ['I' + c.fstring + '.fits', 'I' + c.fstring + '.con', \
             'W' + c.fstring + '.fits', 'E_I' + c.fstring + '.txt', \
             'EM_I' + c.fstring + '.fits', 'G_I' + c.fstring + '.in', \
             'M_I' + c.fstring + '.fits', 'OE_I' + c.fstring + '.fits', \
             'OEM_O_I' + c.fstring + '.fits', 'O_I' + c.fstring + '.fits', \
             'SO_I' + c.fstring + '.fits', 'R_I' + c.fstring + '.html']
    return files

def CrashHandlerToRemove(gal_id):
    """Removing all the temp files"""
    f_R_crash =  'R_I' + c.fstring + '_1.html'
    f_R_cra = open(f_R_crash, 'w')
    P_cra = 'P_I' + c.fstring + '.png'
    P_new = 'P_I' + c.fstring + '_1.png'
    for line_crash in open('R_I' + c.fstring + '.html', 'r'):
        line_crash2wri = line_crash.replace(P_cra, P_new)
        f_R_cra.write(line_crash2wri)
    f_R_cra.close()
    os.rename(P_cra, P_new)
    for f in PyMorphFiles(gal_id):
        if os.path.exists(f):
            os.remove(c.outdir + f)
    for f in glob.glob(c.outdir + 'Tmp*'):
        if os.path.exists(f):
            os.remove(c.outdir + f)

def PyMorphOutputParams(dbparams, decompose = 0):
    """Returns the output parameters from pymorph in two shapes. One shape is for writing and the other is for passing to writehtml"""
    params = {1:['Name','varchar(500)'],2:['ra_gal','float'],3:['dec_gal','float'],
              4:['z','float'], 5:[ 'MorphType','int'],6:[ 'mag_auto','float'],
              7:['magerr_auto','float'],8:['SexHalfRad','float'],
              9:[ 'num_targets','float'],10:[ 'C','float'],
              11:['C_err','float'], 12:[ 'A','float'],13:[ 'A_err','float'],
              14:['S','float'], 15:[ 'S_err','float'],16:[ 'G','float'],
              17:[ 'M','float'], 18:['magzp', 'float']}
    
    if decompose:
        extra_params = [['bulge_xctr','float'],['bulge_xctr_err','float'],
                         ['bulge_yctr','float'],['bulge_yctr_err','float'],
                         [ 'Ie','float'],['Ie_err','float'],
                         [ 'AvgIe','float'],[ 'AvgIe_err','float'],['re_pix','float'],
                         [ 're_pix_err','float'],['re_kpc','float'],
                         [ 're_kpc_err','float' ],['n','float'],['n_err','float'],
                         ['eb','float'],[ 'eb_err','float'],['bpa','float'],
                         ['bpa_err','float'], [ 'bboxy','float'],[ 'bboxy_err','float'],
                         ['disk_xctr','float'],['disk_xctr_err','float'],
                         ['disk_yctr','float'],['disk_yctr_err','float'],
                         ['Id','float'],[ 'Id_err','float'],
                         [ 'rd_pix','float'],['rd_pix_err','float'],[ 'rd_kpc','float'],
                         ['rd_kpc_err','float'],[ 'ed','float'],[ 'ed_err','float'],
                         ['dpa','float'],[ 'dpa_err','float'],['dboxy','float'],
                         [ 'dboxy_err','float'],['BD','float'],[ 'BT','float'],
                         ['p_xctr','float'],['p_xctr_err','float'],
                         ['p_yctr','float'],['p_yctr_err','float'],
                         [ 'Ip','float'],['Ip_err','float'],[ 'Pfwhm','float'],
                         [ 'Pfwhm_kpc','float'],['bar_xctr','float'],
                         ['bar_xctr_err','float'],['bar_yctr','float'],
                         ['bar_yctr_err','float'],['Ibar','float'],
                         [ 'Ibar_err','float'],['rbar_pix','float'],
                         [ 'rbar_pix_err','float'],['rbar_kpc','float'],
                         [ 'rbar_kpc_err','float'],['n_bar','float'],
                         [ 'n_bar_err','float'],['ebar','float'],['ebar_err','float'],
                         [ 'barpa','float'],[ 'barpa_err','float'],[ 'barboxy','float'],
                         [ 'barboxy_err','float'],[ 'chi2nu','float'],['Goodness','float'],
                         [ 'run','int'],['SexSky','float'],[ 'YetSky','float'],
                         ['GalSky','float'],['GalSky_err','float'],['dis_modu','float'],
                         ['distance','float'],[ 'fit','int'],['FitFlag','bigint']]
    else:
        extra_params = []

    extra_params +=[['flag','bigint'],['Manual_flag','bigint'],['Comments','varchar(1000)'], ['Date', 'varchar(50)'], 
                    ['Version','float'],['Filter','varchar(500)'],['Total_Run','int'], ['rootname', 'varchar(500)']]

    # add additional keys to params
    fill_key = max(params.keys()) 
    for add_param in extra_params:
        fill_key +=1
        params[fill_key] = add_param

    fill_key = max(params.keys()) 
    for dbparam in dbparams:
        fill_key +=1
        DBparam = dbparam.split(':')
        try:
            params[fill_key] = [DBparam[0],DBparam[2]]
        except:
            params[fill_key] = [DBparam[0],'varchar(500)']
        
    return params

def ReadGalfitConfig(cfile):
    """Return gimg, oimg, wimg, pfile, mimg, confile from galfit config file"""
    for line_c in open(cfile, 'r'): #Reading config file if it exists
        try:
            valuec = line_c.split()
            if str(valuec[0]) == 'A)':
                gimg = valuec[1]
            if str(valuec[0]) == 'B)':
                oimg = valuec[1]
            if str(valuec[0]) == 'C)':
                wimg = valuec[1]
                if exists(wimg):
                    pass
                else:
                    wimg = 'None'
            if str(valuec[0]) == 'D)':
                pfile = valuec[1]
            if str(valuec[0]) == 'F)':
                mimg = valuec[1]
            if str(valuec[0]) == 'G)':
                cofile = valuec[1]
        except:
            pass
    gimg = os.path.split(gimg)[1]
    wimg = os.path.split(wimg)[1]
    oimg = os.path.split(oimg)[1]
    mimg = os.path.split(mimg)[1]
    pfile = os.path.split(pfile)[1]
    cofile = os.path.split(cofile)[1]
    return gimg, oimg, wimg, pfile, mimg, cofile 

def CheckHeader(header0):
    """Set the global keywords from the fits header"""
    if 'EXPTIME' in header0.keys(): # Old was 'if header0.has_key('EXPTIME'):
        c.EXPTIME = header0['EXPTIME']
    else:
        c.EXPTIME = 1.0
    if 'RDNOISE' in header0.keys():
        c.RDNOISE= header0['RDNOISE']
    else:
        c.RDNOISE = 0.0
    if 'GAIN' in header0.keys():
        c.GAIN = header0['GAIN']
        c.SEx_GAIN = header0['GAIN']
    else:
        c.GAIN = 1.0
        c.SEx_GAIN = 1.0
    if 'NCOMBINE' in header0.keys():
        c.NCOMBINE= header0['NCOMBINE']
    else:
        c.NCOMBINE = 1
    if 'FILTER2' in header0.keys() or 'FILTER' in header0.keys():
        try:
            c.FILTER = header0['FILTER2']
        except:
            c.FILTER = header0['FILTER']

def FindCutSize(ReSize, VarSize, Square, FracRad, FixSize, SizeXX, SizeYY):
    """Return the size of the cutout. SizeXX, SizeYY are the size of 
      cutout if c.galcut is True"""
    ArcR = c.pos_ang * Get_R() #pa in radian
    SizeX = c.SexHalfRad * FracRad * np.abs(np.cos(ArcR)) + \
            c.axis_rat * c.SexHalfRad * FracRad * np.abs(np.sin(ArcR))
    SizeY = c.SexHalfRad * FracRad * np.abs(np.sin(ArcR)) + \
            c.axis_rat * c.SexHalfRad * FracRad * np.abs(np.cos(ArcR))
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
            # FIX
            SizeX = SizeXX
            SizeY = SizeYY
            # END
    else:
        if VarSize:
            pass
        else:
            SizeX = FixSize
            SizeY = FixSize
    return SizeX, SizeY

def MakeCutOut(xcntr, ycntr, alpha_j, delta_j, SizeX, SizeY, TX, TY, cutimage, whtimage, ReSize):
    """Make cutout image. The xcntr, ycntr are like iraf.  SizeX, SizeY are half size"""
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
    # FIX check c.imagedata or c.ggimage
    CutImDa = c.imagedata[ymin:ymax, xmin:xmax].copy()
    # END
    hdu = pyfits.PrimaryHDU(CutImDa.astype(np.float32))
    SizeY, SizeX = CutImDa.shape
    try:
        hdu.header.update('RA_TARG', alpha_j)
        hdu.header.update('DEC_TARG', delta_j)
    except:
        print 'Problem updating the Ra and Dec in cutout image'
    if c.EXPTIME != -9999:
        hdu.header.update('EXPTIME', c.EXPTIME)
    else:
        print 'c.EXPTIME have value -9999. Something wrong?'
    if c.RDNOISE != -9999:
        hdu.header.update('RDNOISE', c.RDNOISE)
    else:
        print 'c.RDNOISE have value -9999. Something wrong?'
    if c.GAIN != -9999:
        hdu.header.update('GAIN', c.GAIN)
    else:
        print 'c.GAIN have value -9999. Something wrong?'
    if c.NCOMBINE != -9999:
        hdu.header.update('NCOMBINE', c.NCOMBINE)
    else:
        print 'c.NCOMBINE have value -9999. Something wrong?'
    hdu.writeto(os.path.join(c.datadir, cutimage))
    # FIX
    #Making weight image cut
    if c.weightexists:
        z2 = c.weightdata[ymin:ymax,xmin:xmax].copy()
        hdu = pyfits.PrimaryHDU(z2.astype(np.float32))
        hdu.writeto(os.path.join(c.datadir, whtimage))
    else:
        print 'Cannot creat weight image. If you supply weight image please ',\
              'check whether it exists or report a bug'
    #END
    return cut_xcntr, cut_ycntr, SizeX, SizeY, ExceedSize 

def WriteError(err):
    """Write error.log file"""
    f_err = open('error.log', 'a')
    f_err.write(err) 
    f_err.close()

def FitEllipseManual(cutimage, xcntr, ycntr, SizeX, SizeY, sky, out):
    """Find 1-d profile of image"""
    print 'Working on FitEllipseManual '
    #print cutimage, xcntr, ycntr, SizeX, SizeY, sky, out
    if out:
        ell_mask_file = 'OEM_' + c.fstring + '.fits'
        ell_out = 'OE_' + c.fstring + '.txt'
        cutimage = c.outdir + cutimage 
    else:
        ell_mask_file = 'EM_' + c.fstring + '.fits'
        ell_out = 'E_' + c.fstring + '.txt'
        cutimage = os.path.join(c.datadir, cutimage)
    if exists(cutimage) and exists(ell_mask_file):
        f = pyfits.open(cutimage)
        galaxy = f[0].data
        galaxy = galaxy - sky
        f.close()
        f = pyfits.open(ell_mask_file)
        mask = f[0].data
        f.close()
        x = np.reshape(np.arange(SizeX * SizeY), (SizeY, SizeX)) % SizeX
        x = x.astype(np.float32)
        y = np.reshape(np.arange(SizeX * SizeY), (SizeY, SizeX)) / SizeX
        y = y.astype(np.float32)
        # r is the radius parameter
        co = np.cos(c.SexPosAng * Get_R())
        si = np.sin(c.SexPosAng * Get_R())
        xsq = ((x - xcntr)* co + (y - ycntr) * si)**2.0
        ysq = ((xcntr - x) * si + (y - ycntr) * co)**2.0
        one_minus_eb_sq = (1.0 - c.axis_rat)**2.0
        r = np.sqrt(xsq + ysq / one_minus_eb_sq)
        maskedgalaxy = ma.masked_array(galaxy, mask)
        R = []
        IntR = []
        IntRE = []
        MaxRad = np.min([np.log10(8 * c.SexHalfRad), \
                 np.log10(np.min(galaxy.shape))])
        NoOfPoints = int(30 * 10**MaxRad / 50.)
        #print MaxRad, NoOfPoints
        # FIX The EXPTIME factor in the error and intensity. Otherwise the 
        # S/N will be different
        for i in np.logspace(0, MaxRad, NoOfPoints, endpoint=True):
            Isub = maskedgalaxy[np.where(np.abs(r - i) <= 1.0)]
            NonMaskNo = len(ma.compressed(Isub))
            if NonMaskNo > 0 and ma.sum(Isub) > 0.0:
                R.append(i)
                IntRE.append(np.sqrt(ma.sum(Isub)) / (1.0 * NonMaskNo))
                IntR.append(np.mean(Isub))
            #print i, NonMaskNo
            # If you want to see the ellipse anulus, uncoment the following
            # START
            #if i > 10 and out:
            # mafilled = maskedgalaxy.copy()
            # mafilled = ma.filled(mafilled, 0)
            # mafilled[np.where(np.abs(r - i) <= 1.0)] = 10000
            # hdu = pyfits.PrimaryHDU(mafilled)
            # if exists('try.fits'):
            #  os.remove('try.fits')
            # hdu.writeto('try.fits')
            # raw_input()
            # STOP
        IntR = np.array(IntR)
        IntRE = np.array(IntRE)
        # END
        mag = -2.5 * np.log10(IntR)
        mag_l = np.abs(-2.5 * (np.log10(IntR) - np.log10(IntR - IntRE))) 
        mag_u = np.abs(-2.5 * (np.log10(IntR + IntRE) - np.log10(IntR)))
        f = open(ell_out, 'w')
        writer = csv.writer(f, delimiter=' ')
        writer.writerow(['sma', 'inte', 'intee', 'mag', 'magl', 'magu'])
        for i in range(len(IntR)):
            p = [R[i], IntR[i], IntRE[i], mag[i], mag_l[i], mag_u[i]]
            writer.writerow(p)
        f.close()
        print 'Done' 
def CleanEllipse(ell_out, after):
    """Cleaning temp files from Ellipse task. after=1 means cleaning after \
       the task. before = 0"""
    if after:
        ftorlist = ['ellip', 'err', 'test.tab', 'GalEllFit.fits', \
                'GalEllFit.fits.pl']
    else:
        ftorlist = [ell_out, 'ellip', 'err', 'test.tab', 'GalEllFit.fits', \
                'GalEllFit.fits.pl']
    for ftor in ftorlist:
        if os.access(ftor, os.F_OK):
            os.remove(ftor)
 
def HandleEllipseTask(cutimage, xcntr, ycntr, SizeX, SizeY, sky, out):
    """Running the ellipse task. SizeX, SizeY are the total size"""
    manual_profile = 0
    try:
        import unknown
        from pyraf import iraf
        from fitellifunc import run_elli
        use_pyraf = 1
    except ImportError:
        use_pyraf = 0
        print 'No pyraf installed!'
        WriteError('Cannot find pyraf installation! Trying manual 1d ' + \
                   'profile finder\n')
    if use_pyraf:
        if out:
            ell_mask_file = 'OEM_' + c.fstring + '.fits'
            ell_out = 'OE_' + c.fstring + '.txt'
        else:
            ell_mask_file = 'EM_' + c.fstring + '.fits'
            ell_out = 'E_' + c.fstring + '.txt'
        plfile = 'GalEllFit.fits.pl'
        CleanEllipse(ell_out, 0)
        try:
            iraf.imcopy(ell_mask_file, plfile, verbose='no')
            iraf.flpr()
        except:
            pass
        try:
            run_elli(cutimage, ell_out, xcntr, ycntr, c.eg, \
                     c.pos_ang, c.major_axis, sky)
            CleanEllipse(ell_out, 1)
            try:
                iraf.flpr()
            except:
                pass
            if exists(ell_out):
                pass
            else:
                manual_profile = 1
        except:
            manual_profile = 1
            WriteError('Error in ellipse task. Trying manual profile finder\n')
            c.Flag = SetFlag(c.Flag, GetFlag('ELLIPSE_FAIL'))
    if use_pyraf == 0 or manual_profile:
        FitEllipseManual(cutimage, xcntr, ycntr, SizeX, SizeY, sky, out)

def HandleGalfitOutput(cutimage, outimage, xcntr, ycntr,  SizeX, SizeY, line_s):
    """Handling the output images of galfit inclding ellipse fit"""
    if exists(outimage):
        if(c.repeat == False):
            OutMaskFunc(outimage, xcntr, ycntr,  SizeX, SizeY, line_s)
        outmodel = 'S' + outimage
        if os.access(outmodel, os.F_OK):
            os.remove(outmodel)
        # FIX. Add the header also in the file
        f = pyfits.open(outimage)
        z  = f[2].data
        f.close()
        hdu = pyfits.PrimaryHDU(z)
        hdu.writeto(outmodel)
        # END
        # FIX . Now sky = c.SexSky is passing to the ellipse task. This \
        # means that 
        # we are comparing images with same sky. may be fixed later
        HandleEllipseTask(outmodel, xcntr, ycntr, SizeX, SizeY, c.SexSky, 1)
        # END 
    else:
        WriteError('The output image is not found ' + \
                      'GALFIT MIGHT BE CRASHED\n')
        c.Flag = SetFlag(c.Flag, GetFlag('GALFIT_FAIL'))
        FailedGalfit(cutimage)
        c.run = 0

def Distance(psffile, ra, dec):
    """Find the distance between psf and object in arcsec. Ra and dec is the \
       position of the object. Psf coordinates will be read from the header"""
    try:
        p = pyfits.open(psffile)
        header = p[0].header
        if(header.has_key('RA_TARG')):
            ra_p = header['RA_TARG']
        if (header.has_key('DEC_TARG')):
            dec_p = header['DEC_TARG']
        p.close()
        distance = 3600.0 * np.sqrt((dec - dec_p)**2.0 + \
                   ((ra - ra_p) * np.cos(dec * Get_R()))**2.0)
    except:
        distance = 9999
    if distance == 9999:
        print 'Distance between psf and image is 9999'
    return distance
    
def HandlePsf(cfile, UserGivenPsf, ra, dec):
    """Determine the psf used for fitting""" 
    if not c.repeat:
        if UserGivenPsf != 'None':
            psffile = UserGivenPsf
            UpdatePsfRaDec(psffile)
            distance = Distance(psffile, ra, dec)
        elif np.abs(ra) == 9999 or np.abs(dec) == 9999:
            psffile = c.psflist[c.psfcounter]
            UpdatePsfRaDec(psffile)
            distance = 9999
            c.psfcounter += 1
        else:
            psffile, distance = SelectPsf(ra, dec)
            distance = distance * 60.0 * 60.0
    else:
        if np.abs(ra) == 9999 or np.abs(dec) == 9999:
            distance = 9999
            psffile = c.pfile
        else:
            psffile = c.pfile
            UpdatePsfRaDec(c.pfile)
            distance = Distance(c.pfile, ra, dec)
    return psffile, distance 

def HandleCASGMBack(cutimage, cut_xcntr, cut_ycntr, SizeX, SizeY, \
                    line_s, bxcntr, bycntr):
    """Find the blank sky region for casgm"""
    try:
        ElliMaskFunc(cutimage, cut_xcntr, cut_ycntr, \
                     SizeX, SizeY, line_s, 0)
        try:
            Bkgd_Params = BkgdFunc(cutimage, cut_xcntr, cut_ycntr, \
                                   bxcntr, bycntr, c.eg, c.pos_ang, \
                                   c.SexSky)
            bxcntr = Bkgd_Params.bkgd[0]
            bycntr = Bkgd_Params.bkgd[1]
#            c.skysig = Bkgd_Params.bkgd[2]
#            print 'Sky Sigma >>> ', c.skysig
            return bxcntr, bycntr         
        except:
            WriteError('Could not find blank sky region for casgm\n')
            return 9999, 9999
    except:
        WriteError('Could not create mask for casgm to find the sky' +\
                   ' sigma and mean \n')
        return 9999, 9999


def HandleCasgm(cutimage, xcntr, ycntr, alpha_j, delta_j, redshift, SizeX, SizeY, line_s, bxcntr, bycntr):
    """Run casgm module and its associated functions"""
    #Finding blank sky
    mask_file = 'EM_' + c.fstring + '.fits'
    if bxcntr == 9999 or bycntr == 9999:
        bxcntr, bycntr = HandleCASGMBack(cutimage, xcntr, ycntr, SizeX, \
                                         SizeY, line_s, bxcntr, bycntr)
    else:
        pass
    try:
        caSgm = casgm(cutimage, mask_file, xcntr, ycntr, bxcntr, \
                      bycntr, c.eg, c.pos_ang, c.SexSky, c.SkySig)
        C = caSgm[0]
        C_err = caSgm[1]
        A = caSgm[2]
        A_err = caSgm[3]
        S = caSgm[4]
        S_err = caSgm[5]
        G = caSgm[6]
        M = caSgm[7]
        print 'C, C_err, A, A_err, S, S_err, G, M >>> ', \
              str(C)[:5], str(C_err)[:5], str(A)[:5], str(A_err)[:5], \
              str(S)[:5], str(S_err)[:5], str(G)[:5], str(M)[:5]
        
        #if(c.decompose == False):
        #    f_res = open("result.csv", "ab")
        #    writer = csv.writer(f_res)
        #    GalId = c.fstring
        #    writer.writerow([GalId, alpha_j, delta_j, redshift, c.SexMagAuto, \
        #                     c.SexMagAutoErr, C, C_err, A, A_err, S, \
        #                     S_err, G, M, c.Flag, c.SexHalfRad])
        #    f_res.close()
        WriteError('(((((CASGM Successful)))))\n')
        return C, C_err, A, A_err, S, S_err, G, M
    except:
        WriteError('The CASGM module failed\n')
        c.Flag = SetFlag(c.Flag, GetFlag('CASGM_FAIL'))
        return 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999

def returngimg():
    print "No gimg given."
    if exists(os.path.join(c.datadir, 'I' + c.fstring + '.fits')):
        gimg = 'I' + c.fstring + '.fits'
    elif exists(os.path.join(c.datadir, str(gal_id) + '.fits')):
        gimg = str(gal_id) + '.fits'
    else:
        print "No possible gimg found"
        gimg = 'None'
    return gimg

def returnthenames(gal_id, gimg, wimg):
    """This function returns the names of cutout, weight \
       based on the parameters set in the config.py"""
    if c.galcut:
        if c.ReSize:
            cutimage = 'I' + gimg
            whtimage = 'W' + wimg 
        else:
            cutimage = gimg
            whtimage = wimg
    else:
        cutimage = 'I' + gimg
        whtimage = 'W' + wimg
    return gimg, wimg

