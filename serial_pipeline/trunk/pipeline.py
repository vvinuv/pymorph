#!/usr/bin/env python
"""This is the python pipeline for GALFIT """

import os
from os.path import exists
import sys
import csv
import pyfits
import numarray as n
from pyraf import iraf
import config as c
from maskfunc import *
from configfunc import *
from ellimaskfunc import *
from outmaskfunc import *
from fitellifunc import *
from plotfunc import *
from writehtmlfunc import *
from runsexfunc import *
from casgm import *

def main():
    imagefile = c.imagefile
    whtfile = c.whtfile
    sex_cata = c.sex_cata
    clus_cata = c.clus_cata
    out_cata = c.out_cata
    size = c.size/2
    threshold = c.threshold
    thresh_area = c.thresh_area
    if exists('index.html'):
        pass
    else:
        indexfile = open('index.html', 'w')
        indexfile.writelines(['<HTML>\n<BODY>\n'])
        indexfile.writelines(['</BODY></HTML>'])
        indexfile.close()
    try:
        if(c.repeat == False and c.galcut == False):
            img = pyfits.open(imagefile)
            image = img[0].data
            img.close()
            print imagefile
    except IOError, (errno, strerror):
        print imagefile, "I/O error(%s): %s" % (errno, strerror)
        os._exit(0)
    try:
        if exists(whtfile):
            if(c.repeat == False and c.galcut == False):
                wht = pyfits.open(whtfile)
                weight = wht[0].data
                wht.close()
                print whtfile
        else:
           print 'No weight image found\n'
    except IOError, (errno, strerror):
        print whtfile, "I/O error(%s): %s" % (errno, strerror)
        pass
    psflist = c.psflist
    try:        #The function which will update the psf header if the psf files
                #are the specified format
        for element in psflist:
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
            print element
            iraf.hedit(element, 'RA_TARG', ra, add= 'yes', verify= 'no', \
                       update='yes')
            iraf.hedit(element, 'DEC_TARG', dec, add= 'yes', verify= 'no', \
                       update='yes')
    except:
        pass
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
        """This function will select the nearest psf from the psflist"""
        distance = 9999.0
        psffile = 'test.fits'
        psflist = c.psflist
        r = 3.14159265 / 180.0
        for element in psflist:
            p=pyfits.open(element)
            header = p[0].header
            if (header.has_key('RA_TARG')):
                ra = header['RA_TARG']
            if (header.has_key('DEC_TARG')):
                dec= header['DEC_TARG']
            p.close()
#		d = sqrt((ra - alpha_j) ** 2.0 + (dec - delta_j) ** 2.0)
            d = n.arccos(n.cos((90.0 - delta_j) * r) * n.cos((90.0 - dec) *\
                r) + n.sin((90.0 - delta_j) * r) *  n.sin((90.0 - dec) * r) * \
                n.cos((alpha_j - ra) * r))
            print 'alp dec alpsf decpsf d', alpha_j, delta_j, ra, dec, d
            if(d < distance):
                psffile = element
                distance = d
        return psffile, distance

#weight = where(weight1 > 0, 1.0 / sqrt(weight1), 0.0)
    if exists('result.csv'):
        pass
    else:
        f_res = open("result.csv", "ab")
        writer = csv.writer(f_res)
        if(c.decompose):
            writer.writerow(['Name','ra','dec','z','Ie','Ie_err','re(pixels)',\
                         're_err(pixels)', 're(kpc)', 're_err(kpc)' ,'n',\
                         'n_err','Id','Id_err','rd(pixels)','rd_err(pixels)',\
                         'rd(kpc)', 'rd_err(kpc)', 'BD', 'BT', 'chi2nu', \
                         'run', 'C', 'C_err', 'A', 'A_err', 'S', 'S_err', \
                         'G', 'M', 'Comments'])
        else:
            writer.writerow(['Name','ra','dec','z', 'C', \
                         'C_err', 'A', 'A_err', 'S', 'S_err', 'G', 'M', \
                         'Comments'])
        f_res.close()
    f_cat = open(out_cata,'w')
    obj_file = open(clus_cata,'r')  #The file contains the objects of interest
    pnames = obj_file.readline().split() #The names of the parameters given 
                                         #in the first line in the clus_cata
    pdb = {}                        #The parameter dictionary
    psfcounter = 0                  #For getting psf in the case of unknown ra
                                    #and dec
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
                gal_id = pdb["gimg"][:-5] #id will be filename without .fits
            try:
                alpha1 = float(pdb["ra1"])
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
                gimg = 'None'
            try:
                wimg = pdb["wimg"]   #Weight cut
            except:
                wimg = 'None'
            try:
                cfile = pdb["cfile"]  #GALFIT configuration file
            except:
                if(c.repeat == True and c.galcut == False):
                    cfile = 'G_I' + str(c.rootname) + '_' + \
                             str(gal_id) + '.in'
                elif(c.repeat == True and c.galcut == True):
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
            try:
                ximg = pdb["ximg"]
            except:
                if(c.galcut == True):
                    ggimg = pyfits.open(gimg)
                    ggimage = ggimg[0].data
                    ggimg.close()
                    c.size = ggimage.shape[0] 
                    ximg = c.size / 2.0
                else:
                    ximg = -9999
            try:
                yimg = pdb["yimg"]
            except:
                if(c.galcut == True):
                    yimg = c.size / 2.0
                else:
                    yimg = -9999
            if(c.galcut == True):   #Given galaxy cutouts
                if exists(sex_cata): #If the user provides sextractor catalogue
                                     #then it will not run SExtractor else do!
                    pass
                else: 
                    RunSex(gimg, wimg)
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
                    if(c.galcut == False):
                        alpha_s = float(values[3]) - (c.shiftra) #This is the difference between the observed and the published coordinate for an object. It is used to correct the sextractor cordinate to compare with the published one.
                        delta_s = float(values[4]) - (c.shiftdec) 
                    else:
                        alpha_s = 9999
                        delta_s = 9999
                    sex_id = values[0]
                    if(c.galcut == True):
                        xcntr  = float(values[1])
                        ycntr  = float(values[2])
                    else:
                        xcntr  = 9999
                        ycntr  = 9999
                    if(abs(alpha_j - alpha_s) < 0.00027/1.0 and \
                       abs(delta_s - delta_j) < 0.00027/1.0 or \
                       abs(xcntr - ximg) < 3.5 and abs(ycntr - yimg) < 3.5):
                        f_err = open('error.log', 'a') 
                        if(c.galcut == True):
                            cutimage = gimg
                            whtimage = wimg
                        else:
                            cutimage = 'I' + str(c.rootname) + '_' + \
                                        str(gal_id) + '.fits'
                            whtimage = 'W' + str(c.rootname) + '_' + \
                                        str(gal_id) + '.fits'
                            xcntr  = float(values[1])
                            ycntr  = float(values[2])
                            xmin = int(xcntr) - (size + 1)
                            ymin = int(ycntr) - (size + 1)
                            xmax = int(xcntr) + (size - 1)
                            ymax = int(ycntr) + (size - 1)
                        f_err.writelines(['\n\n###########   ', str(gal_id), \
                                          '   ###########\n'])
                        run = 1 #run =1 when pipeline runs sucessfuly
                        try:
                            if(c.repeat == False and c.galcut == False):
                                z1 = image[ymin:ymax,xmin:xmax]
                                hdu = pyfits.PrimaryHDU(z1.astype(Float32))
                                try:
                                    hdu.header.update('RA_TARG', alpha_j)
                                    hdu.header.update('DEC_TARG', delta_j)
                                except:
                                    pass
                                hdu.writeto(cutimage)
                            try:
                                if(c.repeat == False and c.galcut == False):
                                    if exists(whtfile): 
                                        z2 = weight[ymin:ymax,xmin:xmax]
                                        hdu = pyfits.PrimaryHDU(z2.astype\
                                              (Float32))
                                        hdu.writeto(whtimage)
                                if(c.galcut == False):
                                    xcntr_o  = xcntr #x center of the object
                                    ycntr_o  = ycntr #y center of the object
                                    xcntr = size + 1.0 + xcntr_o - int(xcntr_o)
                                    ycntr = size + 1.0 + ycntr_o - int(ycntr_o)
                                mag    = float(values[7]) #Magnitude
                                radius = float(values[9]) #Half light radius
                                mag_zero = 25.256 #magnitude zero point
                                sky  = float(values[10]) #sky
                                pos_ang = pa(values[11])
                                axis_rat = 1.0 / float(values[12]) #axis ration b/a
                                eg = 1 - axis_rat
                                if(eg<=0.05):
                                    eg = 0.07
                                major_axis = float(values[14])
                                #major axis of the object
                                if(c.decompose):
                                    try:
                                        if(c.repeat == False):
                                            ElliMaskFunc(cutimage, c.size, \
                                                         line_s, 1)
                                        ell_mask_file = 'EM_' + \
                                                         str(cutimage)[:-5] + \
                                                        '.fits'
                                        plfile = str(cutimage) + '.pl'
                                        if os.access(plfile, os.F_OK):
                                            os.remove(plfile)
                                        try:
                                            iraf.imcopy(ell_mask_file, \
                                                       plfile)
                                            try:
                                                ell_out = 'E_' + \
                                                     str(cutimage)[:-4] + 'txt'
                                                if os.access(ell_out, os.F_OK):
                                                    os.remove(ell_out)
                                                run_elli(cutimage, ell_out,\
                                                             xcntr, ycntr, eg, \
                                                            pos_ang, major_axis)
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
                                ElliMaskFunc(cutimage, c.size, line_s, 0)
                                try:
                                    if(c.decompose == False):
                                        ElliMaskFunc(cutimage, c.size, \
                                                         line_s, 1)
                                    ell_mask_file = 'EM_' + \
                                                      str(cutimage)[:-5] + \
                                                      '.fits'
                                    try:
                                        caSgm = casgm(cutimage, ell_mask_file,\
                                              xcntr, ycntr, eg, pos_ang, sky)
                                        C = caSgm[0]
                                        C_err = caSgm[1]
                                        A = caSgm[2]
                                        A_err = caSgm[3]
                                        S = caSgm[4]
                                        S_err = caSgm[5]
                                        G = caSgm[6]
                                        M = caSgm[7]
                                        print C, C_err, A, A_err, S, S_err, G,M
                                        if(c.decompose == False):
                                            f_res = open("result.csv", "ab")
                                            writer = csv.writer(f_res)
                                            GalId = str(cutimage)[:-5]
                                            writer.writerow([GalId, alpha_j, \
                                            delta_j, z, C, C_err, A, A_err, S, \
                                            S_err, G, M])
                                            f_res.close()
                                        f_err.writelines(['(((((CASGM '\
                                                          'Successful)))))'])
                                    except:
                                        f_err.writelines(['The CASGM module',\
                                                          ' failed\n'])   

                                except:
                                    f_err.writelines(['Could not make mask ',\
                                                      'image for casgm\n'])
                            except:
                                f_err.writelines(['Could not create mask ',\
                                                  'for casgm to find the sky'\
                                                  ' sigma and mean. Remove '\
                                                   'if BMask.fits exists\n'])
                        f_err.close()
                        os.system('rm -f BMask.fits MRotated.fits \
                                  MaskedGalaxy.fits Rotated.fits')
                        if(c.decompose == False):
                            if os.access(ell_mask_file, os.F_OK):
                                os.remove(ell_mask_file)
                        f_err = open('error.log', 'a') 
                        if(c.decompose):
                            try:
                                if(c.repeat == False and cfile == 'None'):
                                    if(alpha_j == -9999 or delta_j == -9999):
                                        psffile = c.psflist[psfcounter]
                                        distance = 9999
                                        psfcounter += 1
                                    else:
                                        psffile = psf_select(alpha_j, delta_j)\
                                                  [0]
                                        distance = psf_select(alpha_j, delta_j)\
                                                  [1] * 60.0 * 60.0
                                else:
                                    if(alpha_j == -9999 or delta_j == -9999):
                                        distance = 9999
                                    else:
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
                                            distance = n.arccos(n.cos((90.0 - \
                                                delta_j) \
                                             * r) * n.cos((90.0 - dec_p) * r) \
                                             + n.sin((90.0 - delta_j) * r) *  \
                                             n.sin((90.0 - dec_p) * r) * \
                                             n.cos((alpha_j - ra_p) * r))
                                            print 'alp dec alpsf decpsf d', alpha_j, delta_j, ra_p, dec_p, distance
                                if(cfile == 'None'):
                                    MaskFunc(cutimage, c.size, line_s)
                                    maskimage = 'M_' + str(cutimage)[:-5] +\
                                                '.fits'
                                else:
                                    maskimage = mimg
                                try:
                                    if(cfile == 'None'):
                                        ConfigFunc(cutimage, whtimage, c.size,\
                                                   line_s, psffile)
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
                                    if(c.galfit):
                                        cmd = str(c.GALFIT_PATH) + ' ' + config_file
                                        os.system(cmd)

#                                        os.system('/Vstr/vstr/vvinuv/galfit/modified/galfit "' + config_file + '"')
                                    if exists('fit.log'):
                                        for line in open('fit.log','r'):
                                            f_fit.writelines([str(line)])
                                    f_fit.close()
                                    try:
                                        if(c.repeat == False):
                                            OutMaskFunc(outimage, c.size, \
                                                       line_s)
                                        out_mask_file = 'OEM_' + \
                                                         str(outimage)[:-5] + \
                                                        '.fits'
                                        outplfile = str(outimage) + '.pl'
                                        if os.access(outplfile, os.F_OK):
                                            os.remove(outplfile)
                                        try:
                                            iraf.imcopy(out_mask_file, \
                                                       outplfile)

                                            try:
                                                ell_output = 'OE_' + \
                                                     str(cutimage)[:-4] + 'txt'
                                                if os.access(ell_output, \
                                                             os.F_OK):
                                                    os.remove(ell_output)  
                                                outmodel = outimage + '[2]'  
                                                run_elli(outmodel, ell_output,\
                                                         xcntr, ycntr, eg, \
                                                         pos_ang, major_axis)
                                            except:
                                                f_err.writelines(['Error in \
                                                           ellipse '\
                                                          'task. Check ', \
                                                          'whether ' ,\
                                                           str(ell_output) ,\
                                                      ' or ellip or err  or',\
                                                      ' test.tab exists OR ',\
                                                      'GALFIT MIGHT BE '\
                                                      'CRASHED\n'])
                                                run = 0
                                        except:
                                            f_err.writelines(['Exists ',\
                                                       str(outimage),'.pl or ',\
                                                       str(out_mask_file),\
                                                       ' does not exist\n'])  
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
                            if(run == 1):
                                try:
                                    if exists('P_' + str(cutimage)[6:-4] \
                                              + 'png'):	
                                        os.system('rm ''P_' + str(cutimage)\
                                                   [6:-4] + 'png''')
                                    PlotFunc(cutimage, outimage, maskimage)
                                except:
                                    f_err.writelines(['Error in plotting. '])
                                    if(maskimage == 'None'):
                                        f_err.writelines(['Could not find '\
                                                          'Mask image\n'])
                                    run = 0	
                                try:	
                                    write_params(cutimage, distance, alpha1, alpha2, alpha3, delta1, delta2, delta3, z, C, C_err, A, A_err, S, S_err, G, M)
#                                f_err.writelines(['(((((((((( Successful', \
 #                                                     ' ))))))))))\n'])
                                except:
                                    try:
                                        write_params(cutimage, distance, alpha1, alpha2, alpha3, delta1, delta2, delta3, z, 9999, 9999, 9999, 9999, 9999, 9999, 9999, 9999)
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

main()
#Filename conventions
#Suppose imagefile = j8f645-1-1_drz_sci.fits
#and gal_id = 9999, then
#cutimage = image_j8f645_9999.fits, the cut out of the galaxy
#whtimage = wht_j8f645_9999.fits, correspoding weight image for the cuts
#maskfile = mask_j8f645_9999.fits, galfit mask
#ell_mask_file = ell_mask_j8f645_9999.fits, mask for ellipse task
#ell_mask_j8f645_9999.fits will be converted to ell_mask_j8f645_9999.fits.pl for ellipse task
#config_file = gal_j8f645_9999.in , configuration file for GALFIT
#The ouput image from galfit out_j8f645_9999.fits
#The parameters from galfit is image_j8f645_9999.txt, if use the modified GALFIT by Vinu Vikram
#The output parametrs will be append to the file fit2.log
#The process status of the pipeline can be seen in the file error.log
#The ellipse task output of input image elli_j8f645_9999.txt
#The ellipse task output of output image out_elli_j8f645_9999.txt
#The plot of input, output, residue images and the 1-D profile comparison can be seen in plot_j8f645_9999.png
#The html output including the figures and parameters will be in the file result_j8f645_9999.html
#The index file of all the fit will be in index.html
