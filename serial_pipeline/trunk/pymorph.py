#!/usr/bin/env python
"""PyMorph [Py MOrphological Parameters' Hunter], is a pipeline to find the Morphological parameters of galaxy. Written by Vinu Vikram in collaboration with 
Yogesh Wadadekar and Ajit K. Kembhavi. 2008 Feb"""

import os
import time
from os.path import exists
import sys
import csv
import pyfits
import numpy as n
from pyraf import iraf
from ndimage import center_of_mass
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
from bkgdfunc import *


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
            header0 = img[0].header
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
            img.close()
            print "imagefile >>> ", imagefile
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
    if(c.decompose):
        try:   #The function which will update the psf header if the psf files
                #are the specified format
            for element in psflist:
                ra1 = float(str(element)[4:6])
                ra2 = float(str(element)[6:8])
                ra3 = float(str(element)[8:10]) + float(str(element)[10]) / 10.0
                dec1 = float(str(element)[11:-10])
                dec2 = float(str(element)[-10:-8])
                dec3 = float(str(element)[-8:-6]) + \
                       float(str(element)[-6]) / 10.0
                ra = (ra1 + (ra2 + ra3 / 60.0) / 60.0) * 15.0
                if dec1 < 0.0:
                    dec = (dec1 - (dec2 + dec3 / 60.0) / 60.0)
                else:
                    dec = (dec1 + (dec2 + dec3 / 60.0) / 60.0)
                #print element
                iraf.hedit(element, 'RA_TARG', ra, add= 'yes', verify= 'no', \
                           show='no', update='yes')
                iraf.flpr() 
                iraf.hedit(element, 'DEC_TARG', dec, add= 'yes', verify= 'no', \
                           show='no', update='yes')
                iraf.flpr()
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
        """This function will select the nearest psf from the psflist.
           The distance is calculated by using the following equation
           d = Sqrt((dec_a - dec_b) ^ 2 + ((ra_a - ra_b) * sin(0.5) * 
           (dec_a - dec_b)) ^ 2.0 )"""
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
#            d = n.arccos(n.cos((90.0 - delta_j) * r) * n.cos((90.0 - dec) *\
#                r) + n.sin((90.0 - delta_j) * r) *  n.sin((90.0 - dec) * r) * \
#                n.cos((alpha_j - ra) * r))
#            d = n.sqrt((delta_j - dec)**2.0 + ((alpha_j-ra)*n.sin((0.5) *\
#                (delta_j+dec)))**2.0)
            d = n.sqrt((delta_j - dec)**2.0 + ((alpha_j - ra) * \
                n.cos(delta_j * r))**2.0)
            #print 'alp dec alpsf decpsf d', alpha_j, delta_j, ra, dec, d
            if(d < distance):
                psffile = element
                distance = d
        return psffile, distance 
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
            ParamToWrite = ['Name','ra','dec','z', 'Ie','Ie_err','re(pixels)',\
                            're_err(pixels)', 're(kpc)', 're_err(kpc)' ,'n', \
                            'n_err', 'Id','Id_err','rd(pixels)',\
                            'rd_err(pixels)', 'rd(kpc)', 'rd_err(kpc)', 'BD', \
                            'BT', 'Point', 'Point_err', 'Pfwhm', 'Pfwhm(kpc)', \
                            'chi2nu', 'Goodness', 'run', 'C', 'C_err', 'A', \
                            'A_err', 'S', 'S_err', 'G', 'M', 'distance', \
                            'fit', 'flag', 'Comments']
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
            writer.writerow(['Name','ra','dec','z', 'C', \
                         'C_err', 'A', 'A_err', 'S', 'S_err', 'G', 'M', \
                         'flag', 'Comments'])
        f_res.close()
    f_cat = open(out_cata,'w')
    obj_file = open(clus_cata,'r')  #The file contains the objects of interest
    pnames = obj_file.readline().split() #The names of the parameters given 
                                         #in the first line in the clus_cata
    pdb = {}                        #The parameter dictionary
    psfcounter = 0                  #For getting psf in the case of unknown ra
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
                    gal_id = pdb["gimg"][:-5] #id will be filename without .fits
                except:
                    print "No image or gal_id found in the object catalogue." \
                          "Exiting"
                    os._exit(0)
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
                    gimg = 'W' + str(c.rootname) + '_' + str(gal_id) + '.fits'
                else:
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
                    else:
                        GAIN = -9999
                    if (header0.has_key('NCOMBINE')):
                        NCOMBINE= header0['NCOMBINE']
                    else:
                        NCOMBINE = -9999
                    ggimg.close()
                    SizeXX = ggimage.shape[0]
                    SizeYY = ggimage.shape[1]
            try:
                ximg = float(pdb["ximg"])
            except:
                if(c.galcut == True):
                    ximg = SizeXX / 2.0
                else:
                    ximg = -9999
            try:
                yimg = float(pdb["yimg"])
            except:
                if(c.galcut == True):
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
                        SeaDeg =  0.00027
                    if(abs(alpha_j - alpha_s) < SeaDeg and \
                       abs(delta_s - delta_j) < SeaDeg or \
                       abs(xcntr - ximg) < SeaPix and \
                       abs(ycntr - yimg) < SeaPix):
                        print "SExtractor ID >>> ", values[0]
                        mag    = float(values[7]) #Magnitude
                        halfradius = float(values[9]) #Half light radius
                        mag_zero = c.mag_zero #magnitude zero point
                        sky  = float(values[10]) #sky
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
                        try:
                            if(c.repeat == False and c.galcut == False):
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
                                    c.Flag = 16384 
                                    xminOut = xmin
                                    xmin = 0
                                if(ymin < 0):
                                    c.Flag = c.Flag + 32768
                                    yminOut = ymin
                                    ymin = 0
                                if(xmax > TX):
                                    c.Flag = c.Flag + 65536
                                    xmaxOut = xmax
                                    xmax = TX
                                if(ymax > TY):
                                    c.Flag = c.Flag + 131072
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
                                    if xminOut != 0:
                                        xcntr = SizeXB + xminOut
                                    if yminOut != 0:
                                        ycntr = SizeYB + xminOut 
                                    if xmaxOut != 0:
                                        xcntr = xcntr
                                    if ymaxOut != 0:
                                        ycntr = ycntr
                                    xcntr = xcntr + xcntrFrac
                                    ycntr = ycntr + ycntrFrac
                                else:
                                    xcntr = SizeX / 2 + xcntrFrac
                                    ycntr = SizeY / 2 + ycntrFrac
#                                print cutimage,xcntr, ycntr, SizeX, SizeY
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
                                        ell_mask_file = 'EM_' + \
                                                         str(cutimage)[:-5] + \
                                                        '.fits'
                                        plfile = str(cutimage) + '.pl'
                                        if os.access(plfile, os.F_OK):
                                            os.remove(plfile)
                                        try:
                                            iraf.imcopy(ell_mask_file, \
                                                       plfile, verbose='no')
                                            iraf.flpr()
                                            try:
                                                ell_out = 'E_' + \
                                                     str(cutimage)[:-4] + 'txt'
                                                if os.access(ell_out, os.F_OK):
                                                    os.remove(ell_out)
                                                run_elli(cutimage, ell_out,\
                                                         xcntr, ycntr, eg, \
                                                         pos_ang, major_axis)
                                                if os.access(plfile, os.F_OK):
                                                    os.remove(plfile)
                                                iraf.flpr()
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
                                                c.Flag = c.Flag + 4
                                        except:
                                            f_err.writelines(['Exists ',\
                                                       str(cutimage),'.pl or ',\
                                                       str(ell_mask_file),\
                                                       ' does not exist\n'])  
                                            run = 0
                                            c.Flag = c.Flag + 2
                                    except:
                                        f_err.writelines(['Error in making '\
                                                     'mask for ellipse task\n'])
                                        run = 0
                                        c.Flag = c.Flag + 1
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
                                            delta_j, z, C, C_err, A, A_err, S, \
                                            S_err, G, M])
                                        f_res.close()
                                    f_err.writelines(['(((((CASGM '\
                                                      'Successful)))))'])
                                except:
                                    f_err.writelines(['The CASGM module',\
                                                          ' failed\n'])   
                                    c.Flag = c.Flag + 16
                            except:
                                f_err.writelines(['Could not make mask ',\
                                                      'image for casgm\n'])
                                c.Flag = c.Flag + 8
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
                                    if(alpha_s == 9999 or delta_s == 9999):
                                        psffile = c.psflist[psfcounter]
                                        try:
                                            p=pyfits.open(pfile)
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
                                        psffile = psf_select(alpha_j, delta_j)\
                                                  [0]
                                        distance = psf_select(alpha_j, delta_j)\
                                                  [1] * 60.0 * 60.0
                                else:
                                    if(alpha_s == 9999 or delta_s == 9999):
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
                                    if(c.galfit):
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
                                            outmodel = 'S' + outimage 
                                            iraf.imcopy(outimage + '[2]', \
                                                       outmodel, verbose='no')
                                            iraf.flpr()
                                            iraf.imcopy(out_mask_file, \
                                                       outplfile, verbose='no')
                                            iraf.flpr()
                                            try:
                                                ell_output = 'OE_' + \
                                                     str(cutimage)[:-4] + 'txt'
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
                                                run_elli(outmodel, ell_output,\
                                                         MoX, MoY, eg, \
                                                         pos_ang, major_axis)
                                                iraf.flpr()
                                                for myfile in [outplfile, \
                                                               outmodel]:
                                                    if os.access(myfile, \
                                                                 os.F_OK):
                                                        os.remove(myfile) 
                                            except:
                                                f_err.writelines(['Error in '\
                                                          'ellipse '\
                                                          'task. Check ', \
                                                          'whether ' ,\
                                                           str(ell_output) ,\
                                                      ' or ellip or err  or',\
                                                      ' test.tab exists OR ',\
                                                      'GALFIT MIGHT BE '\
                                                      'CRASHED\n'])
                                                run = 0
                                                c.Flag = c.Flag + 512
                                        except:
                                            f_err.writelines(['Exists ',\
                                                       str(outimage),'.pl or ',\
                                                       str(out_mask_file),\
                                                       ' does not exist\n'])  
                                            run = 0
                                            c.Flag = c.Flag + 256
                                    except:
                                        f_err.writelines(['Error in making '\
                                                'out mask for ellipse task\n'])
                                        run = 0 
                                        c.Flag = c.Flag + 128
                                except:
                                    f_err.writelines(['Error in writing',\
                                                      ' configuration file\n'])	
                                    run = 0
                                    c.Flag = c.Flag + 64
                            except:
                                f_err.writelines(['Error in making mask for '\
                                                  'galfit\n'])
                                run = 0
                                c.Flag = c.Flag + 32
#                        if exists('plot_' + str(cutimage)[6:-4] + 'png'):	
#                            os.system('rm ''plot_' + str(cutimage)[6:-4] + 'png''')
                            if(run == 1 or run == 0):
                                try:
                                    if exists('P_' + str(cutimage)[6:-4] \
                                              + 'png'):	
                                        os.system('rm ''P_' + str(cutimage)\
                                                   [6:-4] + 'png''')
                                    GoodNess = PlotFunc(cutimage, outimage, \
                                               maskimage, xcntr, ycntr, skysig)
                                    Goodness = GoodNess.plot_profile
                                except:
                                    f_err.writelines(['Error in plotting. '])
                                    if(maskimage == 'None'):
                                        f_err.writelines(['Could not find '\
                                                          'Mask image\n'])
                                    run = 0	
                                    Goodness = 9999
                                    c.Flag = c.Flag + 1024
                                try:
                                    write_params(cutimage, xcntr, ycntr, \
                                                 distance, alpha_j, \
                                                 delta_j, z, Goodness, \
                                                 C, C_err, A, A_err, S, S_err, \
                                                 G, M)
#                                f_err.writelines(['(((((((((( Successful', \
 #                                                     ' ))))))))))\n'])
                                except:
                                    try:
                                        write_params(cutimage, xcntr, ycntr, \
                                                     distance, alpha_j,\
                                                     delta_j, z, Goodness, \
                                                     9999, 9999, 9999,\
                                                     9999, 9999, 9999, 9999, \
                                                     9999)
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
def selectpsf(ImG, CaT):
    c.psff = []
    im = pyfits.open(ImG)
    image = im[0].data
    im.close()
    AreaOfObj = 40
    def FindPsf(AreaOfObj, CaT):
        for line in open(CaT,'r'):
            values = line.split()
            try:
                if float(values[16]) >= 0.8 and float(values[13]) > AreaOfObj \
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
                    c.psff.append(psffile)
                    psf = image[y2:y1, x2:x1]
                    if os.access(psffile, os.F_OK):
                        os.remove(psffile)
                    hdu = pyfits.PrimaryHDU(psf.astype(n.float32))
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
    PsfList = []
    print 'Checking Started. You can just visually check the psf. You can do', \
          'the thorough checking later'
    for element in c.psff:
        os.system('xpaset -p ds9 frame clear')
        if exists(element):
            os.system('cat ' + str(element) + ' | xpaset ds9 fits')
            os.system('xpaset -p ds9 scale mode zscale')
            os.system('xpaset -p ds9 zoom to fit')
        else:
            print 'The psf you have give is NOT exists!!!'
            pass
        write = raw_input("Do you need this psf? ('y' if yes, 'c'"\
                          " to cancel psf checking) " )
        if write == 'y':
            PsfList.append(element)
#            c.psff.remove(element)
        elif write == 'c':
            for element1 in c.psff:
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
    finish = 0
    TmpPsfList = []
    while finish == 0:
        for element in PsfList:
            if exists(element):
                os.system('xpaset -p ds9 frame clear')
                os.system('cat ' + str(element) + ' | xpaset ds9 fits')
                os.system('xpaset -p ds9 scale mode zscale')
                os.system('xpaset -p ds9 zoom to fit')
            else:
                print 'The psf you have give is NOT exists!!!'
                pass
            write = raw_input("Do you REALLY need this psf? ('y' or," \
                              "'n' or press any key to continue) ") 
            if write == 'y':
                TmpPsfList.append(element)
                fff = open('psflist.list', 'ab')
                fff.writelines([str(element), '\n'])
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

if __name__ == '__main__':
    sex_cata = c.sex_cata
    if exists(sex_cata):
        pass
    elif(c.galcut == False):
        print 'The SExtractor catalogue for your frame is NOT found. ' \
              'One is making using the default values. It is always '\
              'recommended to make SExtractor catalogue by YOURSELF as '\
              'the pipeline keep the sky value at the SExtractor value '\
              'during the decomposition.'
        if exists(c.whtfile):
            RunSex(c.imagefile, c.whtfile)
        else:
            RunSex(c.imagefile, 'None')
    def runpsfselect():
        if(c.galcut):   #Given galaxy cutouts
            obj_file = open(c.clus_cata,'r') 
            pnames = obj_file.readline().split() 
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
                    if exists(sex_cata): 
                        pass
                    else:
                        RunSex(gimg, wimg)
                    try:
                        selectpsf(gimg, sex_cata)
                    except:
                        pass
                    if os.access(sex_cata, os.F_OK):
                        os.remove(sex_cata)
                except:
                    pass
            obj_file.close()  
        else:
            selectpsf(c.imagefile, sex_cata)
    if c.psfselect == 2:
        os.system('ds9 &')
        time.sleep(2)
        runpsfselect()
        os.system('xpaset -p ds9 quit')
        c.psflist = '@psflist.list'
        main()
    elif c.psfselect == 1:
        os.system('ds9 &')
        time.sleep(2)
        runpsfselect()
        os.system('xpaset -p ds9 quit')
    elif c.psfselect == 0:
        main()
