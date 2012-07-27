#!/usr/bin/env python
"""This is the python pipeline for GALFIT """

import os
from os.path import exists
import sys
import pyfits
import numarray as n
from pyraf import iraf
import config as c
from maskfunc import *
from configfunc import *
from ellimaskfunc import *
from fitellifunc import *
from plotfunc import *
from writehtmlfunc import *


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
        img = pyfits.open(imagefile)
        image = img[0].data
        img.close()
        print imagefile
    except IOError, (errno, strerror):
        print imagefile, "I/O error(%s): %s" % (errno, strerror)
        os._exit(0)
    try:
        if exists(whtfile):
            wht = pyfits.open(whtfile)
            weight = wht[0].data
            wht.close()
            print whtfile
        else:
           print 'No weight image found\n'
    except IOError, (errno, strerror):
        print whtfile, "I/O error(%s): %s" % (errno, strerror)
        os._exit(0)

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
            if(d < distance):
                psffile = element
                distance = d
        return psffile, distance

#weight = where(weight1 > 0, 1.0 / sqrt(weight1), 0.0)

    f_cat = open(out_cata,'w')
    for line_s in open(sex_cata,'r'):
        try:
            values = line_s.split()
            alpha_s = float(values[3]) - (c.shiftra) #This is the difference between the observed and the published coordinate for an object. It is used to correct the sextractor cordinate to compare with the published one.
            delta_s = float(values[4]) - (c.shiftdec) 
            sex_id = values[0]
            xcntr  = float(values[1])
            ycntr  = float(values[2])
    	    #print sex_id,xcntr,ycntr 
            for line_j in open(clus_cata,'r'):
                try:
                    values = line_j.split()
                    clus_id = values[0]
                    alpha1 = float(values[1])
                    alpha2 = float(values[2])
                    alpha3 = float(values[3])
                    delta1 = float(values[4])
                    delta2 = float(values[5])
                    delta3 = float(values[6])
                    alpha_j = (alpha1 + (alpha2 + alpha3 / 60.0) / 60.0) * 15.0
                    delta_j = delta1 - (delta2 + delta3 / 60.0) / 60.0
                    if(abs(alpha_j - alpha_s) < 0.00027/1.0 and \
                       abs(delta_s - delta_j) < 0.00027/1.0):
                        f_err = open('error.log', 'a') 
                        xmin = int(xcntr) - (size + 1)
                        ymin = int(ycntr) - (size + 1)
                        xmax = int(xcntr) + (size - 1)
                        ymax = int(ycntr) + (size - 1)
                        cutimage = 'image_' + str(imagefile)[:6] + '_' + \
                                    str(clus_id) + '.fits'
                        whtimage = 'wht_' + str(imagefile)[:6] + '_' + \
                                    str(clus_id) + '.fits'
#                       mask_file = 'mask_' + str(imagefile)[:6] + '_'  + str(clus_id) + '.fits'
                        run = 1 #run =1 when pipeline runs sucessfuly
                        f_err.writelines(['\n\n###########   ', str(clus_id), \
                                          '   ###########\n'])
                        try:
                            z1 = image[ymin:ymax,xmin:xmax]
                            hdu = pyfits.PrimaryHDU(z1.astype(Float32))
                            hdu.header.update('RA_TARG', alpha_j)
                            hdu.header.update('DEC_TARG', delta_j)
                            hdu.writeto(cutimage)
                            try:
                                if exists(whtfile): 
                                    z2 = weight[ymin:ymax,xmin:xmax]
                                    hdu = pyfits.PrimaryHDU(z2.astype(Float32))
                                    hdu.writeto(whtimage)
                                try:
                                    ElliMaskFunc(clus_id, line_s)
                                    values = line_s.split()
                                    ell_mask_file = 'ell_mask_' + \
                                                     str(imagefile)[:6] + \
                                                    '_'  + str(clus_id) + \
                                                    '.fits'
                                    xcntr_o  = xcntr #x center of the object
                                    ycntr_o  = ycntr #y center of the object
                                    xcntr = (size) + 1.0 + xcntr_o - \
                                            int(xcntr_o)
                                    ycntr = (size) + 1.0 + ycntr_o - \
                                            int(ycntr_o)
                                    mag    = float(values[7]) #Magnitude
                                    radius = float(values[9]) #Half light radius
                                    mag_zero = 25.256 #magnitude zero point
                                    sky	 = float(values[10]) #sky
                                    pos_ang = pa(values[11])
                                    axis_rat = 1.0 / float(values[12]) #axis ration b/a
                                    eg = 1 - axis_rat
                                    if(eg<=0.05):
                                        eg = 0.07
                                    major_axis = float(values[14])
                                    #major axis of the object
                                    try:
                                        iraf.imcopy(ell_mask_file, 'image' + \
                                                    str(ell_mask_file)[8:] + \
                                                    '.pl')
                                        try:
                                            run_elli(cutimage, xcntr, ycntr, \
                                                     eg, pos_ang, major_axis)
                                        except:
                                            f_err.writelines(['Error in',\
                                                              ' ellipse',\
                                                              ' task. Check ',\
                                                              'whether elli_',\
                                                              str(cutimage)[6:-4],\
                                                               'txt or ellip ',\
                                                               'or err  or ',\
                                                              'test.tab exists\n'])
                                            run = 0
                                    except:
                                        f_err.writelines(['Exists image',\
                                                           str(ell_mask_file)[8:],\
                                                          '.pl or ',\
                                                           str(ell_mask_file),\
                                                          ' does not exist\n'])  
                                        run = 0
                                except:
                                    f_err.writelines(['Error in making mask for ',\
                                                      'ellipse task\n'])
                                    run = 0
                            except:
                                f_err.writelines(['The file ', str(whtimage), \
                                                  ' exists\n'])	
                                run = 0
                        except:
                            f_err.writelines(['The file ', str(cutimage),\
                                              ' exists\n'])
                            run = 0
                        try:
                            MaskFunc(clus_id, line_s)
                            try:
                                psffile = psf_select(alpha_j, delta_j)[0]
                                distance = psf_select(alpha_j, delta_j)[1] * \
                                           60.0 * 60.0
                                ConfigFunc(clus_id, line_s, psffile)
                                outimage = 'out_' + str(imagefile)[:6] + '_'  + \
                                           str(clus_id) + '.fits[2]'
                                try:
                                    run_elli(outimage, xcntr, ycntr, eg, \
                                             pos_ang, major_axis)
                                except:
                                    f_err.writelines(['Error in ellipse task.'\
                                                      ' Check ', \
                                                      'whether out_elli_' ,\
                                                      str(cutimage)[6:-4] ,\
                                                      'txt or ellip or err  or',\
                                                      'test.tab exists OR',\
                                                      'GALFIT MIGHT BE CRASHED\n'])
                                    run = 0
                            except:
                                f_err.writelines(['Error in writing',\
                                                  ' configuration file\n'])	
                                run = 0
                        except:
                            f_err.writelines(['Error in making mask for galfit\n'])
                            run = 0
                        if exists('plot_' + str(cutimage)[6:-4] + 'png'):	
                            os.system('rm ''plot_' + str(cutimage)[6:-4] + 'png''')
                        if(run == 1):
                            try:	
                                PlotFunc(cutimage)
                                try:	
                                    write_params(cutimage, distance)
                                    f_err.writelines(['(((((((((( Successful', \
                                                      ' ))))))))))\n'])
                                except:
                                    f_err.writelines(['Error in writing html\n'])
                                    run = 0
                            except:
                                f_err.writelines(['Error in plotting\n'])
                                run = 0	
						
#iraf.imcopy(str(imagefile) + '[' + str(xmin) + ':' + str(xmax) + ',' + str(ymin) + ':' + str(ymax) + ']', cutimage)	
#iraf.imcopy(str(whtfile) + '[' + str(xmin) + ':' + str(xmax) + ',' + str(ymin) + ':' + str(ymax) + ']', whtimage)	
#					fitellifunc(clus_id, line_s)

                        f_err.close()
                        f_cat.writelines([str(clus_id), ' '])
                        f_cat.write(line_s)
                        for myfile in ['ellip','err','test.tab']:
                            if os.access(myfile,os.F_OK):
                                os.remove(myfile)
                except:
                    pass
        except:
            pass
    f_cat.close()

main()
#Filename conventions
#Suppose imagefile = j8f645-1-1_drz_sci.fits
#and clus_id = 9999, then
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
