import os
import csv
import numpy as np
import fitsio
from .pymorphutils import HMSToDeg, DMSToDeg

verbose = False

def getpsf(datadir, psflist, which_psf, alpha_j, delta_j):

    """

    This function will select the nearest psf from the psflist.
    The distance is calculated by using the following equation
    d = sqrt((dec_a - dec_b)^2 + ((ra_a - ra_b) * sin(0.5) * 
    (dec_a - dec_b))^2.0 )

    """

    psf_distance_dict = {}
    psffile = 'test.fits'
    r = np.pi / 180.

    for element in psflist:
        p = fitsio.FITS(os.path.join(datadir, element))
        header = p[0].read_header()
        p.close()
        #print(header)
        if 'RA_TARG' in header:
            ra = header['RA_TARG']
        elif 'RA' in header:
            ra = header['RA']
        elif 'CRVAL1' in header:
            ra = header['CRVAL1']
        else:
            ra = 9999
        if 'DEC_TARG' in header:
            dec= header['DEC_TARG']
        elif ('DEC' in header):
            dec= header['DEC']
        elif 'CRVAL2' in header:
            dec = header['CRVAL2']
        else:
            dec= 9999
        #print(ra, dec)
        #print(alpha_j, delta_j)
#       d = sqrt((ra - alpha_j) ** 2.0 + (dec - delta_j) ** 2.0)
#       d = np.arccos(np.cos((90.0 - delta_j) * r) * np.cos((90.0 - dec) *\
#           r) + np.sin((90.0 - delta_j) * r) *  np.sin((90.0 - dec) * r) * \
#           np.cos((alpha_j - ra) * r))
#       d = np.sqrt((delta_j - dec)**2.0 + ((alpha_j-ra)*np.sin((0.5) *\
#           (delta_j+dec)))**2.0)
        d = np.sqrt((delta_j - dec)**2.0 + ((alpha_j - ra) * \
            np.cos(delta_j * r))**2.0)
        if verbose:
            print('d', d)
        psf_distance_dict[element] = d

    items = sorted(psf_distance_dict.items(), key=lambda x: x[1])
    if verbose:
        print(items)
    psffile = items[which_psf][0]
    distance = items[which_psf][1]
    if verbose:
        print(psffile, distance)
    return psffile, distance


def PSFArr(psflist):
    """Return psf list if the given input is a file"""
    #print('P1', psflist)
    for pf in psflist:
        ra1 = float(str(pf)[4:6]) 
        ra2 = float(str(pf)[6:8])
        ra3 = float(str(pf)[8:10]) + float(str(pf)[10]) / 10.0
        dec1 = float(str(pf)[11:-10])
        dec2 = float(str(pf)[-10:-8])
        dec3 = float(str(pf)[-8:-6]) + float(str(pf)[-6]) / 10.0
        ra = HMSToDeg(ra1, ra2, ra3)
        dec = DMSToDeg(dec1, dec2, dec3)

        pfits = fitsio.FITS(pf, 'rw')
        pfits[0].write_key('RA', ra)
        pfits[0].write_key('DEC', dec)
        pfits.close()
 
    return psflist



def UpdatePsfRaDec(datadir, element):
    """The function which will update the psf header if the psf files
       are in the specified format"""
    print('U1', os.path.join(datadir, element))
    print(element, str(element)[4:6])
    if 1:
        ra1 = float(str(element)[4:6])
        ra2 = float(str(element)[6:8])
        ra3 = float(str(element)[8:10]) + float(str(element)[10]) / 10.0
        dec1 = float(str(element)[11:-10])
        dec2 = float(str(element)[-10:-8])
        dec3 = float(str(element)[-8:-6]) + float(str(element)[-6]) / 10.0
        ra = HMSToDeg(ra1, ra2, ra3)
        dec = DMSToDeg(dec1, dec2, dec3)
        pfits = fitsio.FITS(os.path.join(datadir, element), 'rw')
        pfits[0].write_key('RA', ra)
        pfits[0].write_key('DEC', dec)
        pfits.close()
    #except:
    #    print('Problem updating PSF coordinates')
    #    pass



def FindPsf(catalog, stargal_prob, area_obj, star_size):
    '''
    Find psf file list
    '''

    psffile = []
    for line in open(catalog, 'r'):
        values = line.split()
        try:
            bkg  = float(values[10])
            if float(values[16]) >= stargal_prob & \
               float(values[13]) > area_obj &\
               float(values[14]) < 50.0:
                xcntr = float(values[1]) - 1
                ycntr = float(values[2]) - 1
                #size of the psf is 8 times the sigma assuming the star has 
                #Gaussian profile
                PsfSize = np.floor(float(values[14])) * star_size

                x1 = int(xcntr) + (PsfSize / 2)
                x2 = int(xcntr) - (PsfSize / 2)
                y1 = int(ycntr) + (PsfSize / 2)
                y2 = int(ycntr) - (PsfSize / 2)

                ra1 = int(float(values[3]) / 15.0)
                ra2 = int((float(values[3]) / 15.0 - int(float(values[3])\
                      / 15.0)) * 60.0)
                ra3 = (((float(values[3]) / 15.0 - int(float(values[3]) / \
                      15.0)) * 60.0) - ra2) * 60.0

                dec1 = int(float(values[4]))
                dec2 = abs(int((float(values[4]) - dec1) * 60.0))
                dec3 = (abs(((float(values[4]) - dec1) * 60.0)) - dec2) * 60.0

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

                pf = 'psf_{}{}{}{}{}{}.fits'.format(ra11, ra22, ra33, dec11, dec22, dec33)
                psffile.append(pf)
                psf = image[y2:y1, x2:x1]
                psf = psf - bkg
                if os.access(pf, os.F_OK):
                    os.remove(pf)
                hdu = pyfits.PrimaryHDU(psf.astype(np.float32))
                hdu.header.update('XCNTR', int(xcntr))
                hdu.header.update('YCNTR', int(ycntr))
                hdu.writeto(pf)
        except:
            pass

    return psffile



def selectpsf(imgfile, catalog):
    '''
    Select PSF
    '''

    psffile = []
    im = pyfits.open(imgfile)
    image = im[0].data
    im.close()

    while len(psffile) < 5 & area_obj > 10:
        psffile = FindPsf(catalog, stargal_prob, area_obj, star_size)
        area_obj -= 5

    if len(psffile) == 0:
        manualpsf = raw_input("No psf is found. Please enter "\
                              "the psf name >>> ")
        psffile.append(manualpsf)

    if interactive:
        PsfList = []
        TmPLST = []
        for element in psffile:
            TmPLST.append(element)
        print('Checking Started. You can just visually check the psf. You can'\
               ' do the thorough checking later')

        for element in psffile:
            os.system('xpaset -p ds9 frame clear')
            if exists(element):
                os.system('cat ' + str(element) + ' | xpaset ds9 fits')
                os.system('xpaset -p ds9 scale mode zscale')
                os.system('xpaset -p ds9 zoom to fit')
            else:
                print('The psf you have given does NOT exists!!!')

            write = raw_input("Do you need this psf? ('y' if yes, 'c'"\
                              " to cancel psf checking) " )
            if write == 'y':
                PsfList.append(element)
                TmPLST.remove(element)
#                psffile.remove(element)
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
        print('\nFinal Checking Started. If you are using psfselect = 2, be '\
              'carefull!!! This is your last chance for selecting psf. '\
              'Select ONLY GOOD psf. Press "y" to accept the previous psf. '\
             'ENTER to go to the next one. This will continue until you press '\
             '"1" when it asked Finished? ALL THE BEST! \n')
        UrPsfChk = raw_input("Do you want to use your own psf? Enter 'y' or " \
                             "'n' >>> ")
        if UrPsfChk == 'y':
            UrPsf = raw_input("Enter your psf >>> ")

            f = open('psflist.list', 'ab')
            f.writelines([str(UrPsf), '\n'])
            f.close()

            if(galcut):
                ValueS.append(UrPsf)
                fwithpsf = open('CatWithPsf.cat', 'ab')
                for v in ValueS:
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
                    print('The psf you have given is NOT exists!!!')
                    pass
                write = raw_input("Do you REALLY need this psf? ('y' or," \
                                  "'n' or press any key to continue) ")
                if write == 'y':
                    TmpPsfList.append(element)
                    f = open('psflist.list', 'ab')
                    f.writelines([str(element), '\n'])
                    if(galcut):
                        if UpdateCounter:
                            ValueS.append(element)
                            fwithpsf = open('CatWithPsf.cat', 'ab')
                            for v in ValueS:
                                fwithpsf.writelines([str(v), ' '])
                            fwithpsf.writelines(['\n'])
                            fwithpsf.close()
                            UpdateCounter = 0
                    f.close()
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
        for element in psffile:
            f = open('psflist.list', 'ab')
            f.writelines([str(element), '\n'])
        f.close()




def run_psfselect(imagefile, datadir, clus_cata, galcut):

    if galcut:   #Given galaxy cutouts
        obj_file = open(os.path.join(datadir, clus_cata), 'r')
        pnames = obj_file.readline().split()
        pnames.append('star')

        fwithpsf = open('CatWithPsf.cat', 'ab')
        for p in pnames:
            fwithpsf.writelines(['{} '.format(p)])
        fwithpsf.writelines(['\n'])
        fwithpsf.close()

        pdb = {}                        #The parameter dictionary
        for line_j in obj_file:
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
                        gal_id = pdb["gimg"][:-5]
                    except:
                        print("No image or gal_id found in the object" \
                              "catalogue. Exiting")
                        os._exit(0)
                try:
                    #gimg = pdb["gimg"]    #Galaxy cutout
                    cutimage = pdb["gimg"]    #Galaxy cutout
                except:
                    if os.path.exists(os.path.join(datadir,
                                     'I{}_{}.fits'.format(rootname, gal_id))):
                        #gimg = 'I{}_{}.fits'.format(rootname, gal_id)
                        cutimage = 'I{}_{}.fits'.format(rootname, gal_id)
                    elif os.path.exists(os.path.join(datadir,
                                        '{}.fits'.format(gal_id))):
                        #gimg = '{}.fits'.format(gal_id)
                        cutimage = '{}.fits'.format(gal_id)
                    else:
                        print("No image found. Exiting")
                        os._exit(0)
                try:
                    wimg = pdb["wimg"]   #Weight cut
                except:
                    if os.path.exists(os.path.join(datadir,
                                     'W{}_{}.fits'.format(rootname, gal_id))):
                        wimg = 'W{}_{}.fits'.format(rootname, gal_id)
                    else:
                        wimg = 'None'

                #gfits = pyfits.open(os.path.join(datadir, gimg))
                gfits = pyfits.open(os.path.join(datadir, cutimage))
                header = gfits[0].header
                if 'GAIN' in header:
                    SEx_GAIN = header['GAIN']
                else:
                    SEx_GAIN = 1
                gfits.close()

                if os.path.exists(sex_cata):
                    pass
                else:
                    #RS = RunSex(os.path.join(datadir, gimg), os.path.join(datadir, wimg), 'None', 9999, 9999, 0)
                    RS = RunSex(os.path.join(datadir, cutimage), os.path.join(datadir, wimg), 'None', 9999, 9999, 0)
                    RS.sex()

                try:
                    #selectpsf(os.path.join(datadir, gimg), sex_cata)
                    selectpsf(os.path.join(datadir, cutimage), sex_cata)
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
            cmd = 'mv CatWithPsf.cat {}'.format(clus_cata)
            os.system(cmd)
        else:
            pass
    else:
        selectpsf(os.path.join(datadir, imagefile), sex_cata)

