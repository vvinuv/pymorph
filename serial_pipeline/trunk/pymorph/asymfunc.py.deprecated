import os
import numpy as n	
import numpy.ma as ma
import pyfits
from rotate_new import *
import config as c

class asymmetry:
    """Finding Asymmetry parameter. The algorithm is as follows
       1. Rotate the galaxy through 180 degrees about its center. Bilinear 
          interpolation was used to find out the rotated image.

       2. Extract a circular region of the image of size 1.5 times the 
          Petrosian radius of the galaxy.

       3. Find the residue of the two images and find the asymmetry value 
                  A = Sum(abs(I_0 - I_r) / Sum(I_0) 
          where I_0 is the galaxy pixel value and I_r is that of rotated image

        4. Centering correction: 
	  a. Asymmetry is computed for centers at the surrounding eight points 
             in a 3X3 grid
	  b. This procedure repeats until a minimum is found for the asymmetry.

        5. Noise correction:
	  a. The uncorrelated noise can be corrected by substracting the 
             asymmetry of the background.

        6. The final formula to compute asymmetry is 			
	  A = min(Sum(abs(I_0 - I_r) / Sum(abs(I_0)) 
              - min(Sum(abs(B_0 - B_r) / Sum(abs(I_0))
        where B_0 is the background pixel value and I_r is that of rotated 
        background"""		
    def __init__(self, cutimage, maskimage, ini_xcntr, ini_ycntr, pa, eg, r50, extraction_radius, background, angle, flag_image, ABS_ZSUM):
        self.cutimage		= cutimage
        self.maskimage          = maskimage
        self.ini_xcntr		= ini_xcntr 
        self.ini_ycntr		= ini_ycntr
        self.extraction_radius	= extraction_radius
        self.pa			= pa
        self.eg			= eg
        self.r50                = r50
        self.background		= background
        self.one_minus_eg_sq	= (1.0-eg)**2.0
        self.flag_image		= flag_image
        self.angle		= angle
        self.ABS_ZSUM		= ABS_ZSUM
#The flag_image=1 for image asymmetry and flag_image=0 for background
        self.image_asymm=ASYM(self.cutimage, self.maskimage, self.ini_xcntr, self.ini_ycntr, self.pa, self.one_minus_eg_sq, r50, self.background, self.extraction_radius, self.angle, self.flag_image, self.ABS_ZSUM)
        return 

def ASYM(cutimage, maskimage, ini_xcntr, ini_ycntr, pa, one_minus_eg_sq, r50, background, extraction_radius, angle, flag_image, ABS_ZSUM):
    if flag_image == 0:
	maskimage = 'BackMask.fits'
    co = n.cos(pa * n.pi / 180.0)
    si = n.sin(pa * n.pi / 180.0)
    Aabs = n.zeros([9])
    Aabs = Aabs.astype(n.float32)
    rot_sum = n.zeros([9])
    rot_sum = rot_sum.astype(n.float32)
    abs_zsum = n.zeros([9])
    abs_zsum = abs_zsum.astype(n.float32)
    sh_sum = n.zeros([9])
    sh_sum = sh_sum.astype(n.float32)
    absres_sum = n.zeros([9])
    absres_sum = absres_sum.astype(n.float32)
    f=pyfits.open(c.datadir +cutimage)
    z = f[0].data
    z = z - background
    f.close()
    z = n.swapaxes(z, 0, 1)
    ff = pyfits.open(maskimage)
    ZZZ = ff[0].data
    ff.close()
    RoTz = rotate_modi(ZZZ, 180.0, ini_xcntr, ini_ycntr, prefilter=False)
    os.system('rm -f MRotated.fits')
    hdu = pyfits.PrimaryHDU(RoTz)
    hdu.writeto('MRotated.fits')
#    f=pyfits.open("AMask.fits")
#    mask = f[0].data
#    f.close()
    f=pyfits.open("MRotated.fits")
    mask1 = f[0].data
    f.close()
    mask1 = n.swapaxes(mask1, 0, 1)
    f=pyfits.open(maskimage)
    mask2 = f[0].data
    f.close()
    mask2 = n.swapaxes(mask2, 0, 1)
    mask = mask1 + mask2
    maskedgalaxy = ma.masked_array(z, mask)
    z1 = ma.filled(maskedgalaxy, background)
    
#    hdu = pyfits.PrimaryHDU(z)
#    hdu.writeto('Z.fits')
#    hdu = pyfits.PrimaryHDU(mask)
#    hdu.writeto('MM.fits')
    hdu = pyfits.PrimaryHDU(n.swapaxes(z, 0, 1).astype(n.float32))
    hdu.writeto('MaskedGalaxy.fits')
    os.system('rm -f MaskedGalaxy1.fits')
    hdu = pyfits.PrimaryHDU(n.swapaxes(z1, 0, 1).astype(n.float32))
    hdu.writeto('MaskedGalaxy1.fits')
#    os._exit(0)
    NXPTS = z.shape[0]
    NYPTS = z.shape[1]
    #r50=c.r50
    center_rad = 0.01 * r50
#    print 'center_rad', center_rad 
#    center_rad = 0.5#0.01*r50 The centering correction has to be done by calculating asymmetric parameters at different points 0.1% of r50 away around the center
    flag_center = 0 #flag_center=1 if the program finds the exact center else center=0
	#extraction_radius=total_rad ie. 1.5 times the radius at eta=0.2
    x = n.reshape(n.arange(NXPTS * NYPTS), (NXPTS, NYPTS)) / NYPTS
    x = x.astype(n.float32)
    y = n.reshape(n.arange(NXPTS * NYPTS), (NXPTS, NYPTS)) % NYPTS
    y = y.astype(n.float32)
    nn = 0 #This nn is given here and the following loop run either it finds the minimum asymmetry or n=20
    while flag_center == 0:
        flag_out = 0
		#x and y coordinates of different center points in 1-d array
        xcntr = n.reshape(n.array([[ini_xcntr] * 3, [ini_xcntr - center_rad] *\
                3, [ini_xcntr + center_rad] * 3]), (9, )) 
        xcntr = xcntr.astype(n.float32)
        ycntr = n.array([ini_ycntr, ini_ycntr - center_rad, ini_ycntr + \
                center_rad]*3)
        ycntr = ycntr.astype(n.float32)
        for iter in range(9):#The below is done in order to keep the extraction radius inside the image. If the extraction radius goes outside the image then flag_out=1
            if(flag_image):#flag_image will tell whether the input is image (flag_image=1) or background (flag_image=0)
                if(xcntr[iter] - 2 - extraction_radius < 0):
                    extraction_radius = xcntr[iter] - 2
                    flag_out = 1
                if(ycntr[iter] - 2 - extraction_radius < 0):
                    extraction_radius = ycntr[iter] - 2
                    flag_out = 1
                if(xcntr[iter] + 2 + extraction_radius > NXPTS):
                    extraction_radius = NXPTS - xcntr[iter] - 2
                    flag_out = 1
                if(ycntr[iter] + 2 + extraction_radius > NYPTS):
                    extraction_radius = NYPTS - ycntr[iter] - 2
                    flag_out = 1
            tx = (x - xcntr[iter]) * co + (y - ycntr[iter]) * si
            ty = (xcntr[iter] - x) * si + (y - ycntr[iter]) * co
            R = n.sqrt(tx**2.0 + ty**2.0 / one_minus_eg_sq)
#	    print xcntr[iter] , ycntr[iter]
#	    os.system('ds9 MaskedGalaxy.fits MaskedGalaxy1.fits')
            ff = pyfits.open("MaskedGalaxy.fits")
	    ZZZ = ff[0].data
            ff.close()
	    RoTz = rotate_modi(ZZZ, 180.0, xcntr[iter], ycntr[iter], prefilter=False)
            os.system('rm -f Rotated.fits')
	    hdu = pyfits.PrimaryHDU(RoTz)
            hdu.writeto('Rotated.fits')
            f = pyfits.open("Rotated.fits")
            rz = f[0].data
            f.close()
#            f = pyfits.open("Rotated1.fits")
#            rz1 = f[0].data
#            f.close()
#	    hdu = pyfits.PrimaryHDU(rz)
#	    hdu.writeto('RR.fits')
            rz = n.swapaxes(rz, 0, 1)
#	    rz1 = n.swapaxes(rz1, 0, 1)
#	    hdu = pyfits.PrimaryHDU(rz)
#	    hdu.writeto('RRSwap.fits')
	    res = z - rz
#            res1 = z1 - rz1
#	    print ini_xcntr, ini_ycntr, flag_image
	    if flag_image < 0:
#	    if ini_xcntr>-145 or ini_xcntr<1135:
 	        do = raw_input('Image ? ')
	        if len(do) == 0:
		    pass
  	        else:
		    os.system('rm -f Z.fits RRSwap.fits BackMM.fits ReS.fits')
 		    hdu = pyfits.PrimaryHDU(z)
	 	    hdu.writeto('Z.fits')
		    hdu = pyfits.PrimaryHDU(rz)
		    hdu.writeto('RRSwap.fits')
                    mm = ma.masked_array(res, mask)
                    mm1 = ma.filled(mm, -1)
		    hdu = pyfits.PrimaryHDU(mm1)
		    hdu.writeto('BackMM.fits')
		    ReS = n.where(R <= extraction_radius, res, 1)
		    hdu = pyfits.PrimaryHDU(ReS)
		    hdu.writeto('ReS.fits')
		    os.system('ds9 Z.fits RRSwap.fits BackMM.fits ReS.fits') 

#	    os._exit(0)
            for myfile in ['Rotated.fits', 'Rotated1.fits', 'Residual.fits', 'AResidual.fits']:
                if os.access(myfile, os.F_OK):
                    os.remove(myfile)
#            sh_sum[iter] = z1[n.where(R <= extraction_radius)].sum()
	    sh_sum[iter] = ma.masked_array(z, mask)[n.where(R <= extraction_radius)].sum()
#            print sh_sum[iter], sh_sum1
#            rot_sum[iter] = rz1[n.where(R <= extraction_radius)].sum()
	    rot_sum[iter] = ma.masked_array(rz, mask)[n.where(R <= extraction_radius)].sum()
#	    print rot_sum[iter], rot_sum1
#            abs_zsum[iter] = abs(z1[n.where(R <= extraction_radius)]).sum()
            abs_zsum[iter] = abs(ma.masked_array(z, mask)[n.where(R <= extraction_radius)]).sum()
#	    print abs_zsum[iter], abs_zsum1
#            absres_sum[iter] = abs(res1[n.where(R <= extraction_radius)]).sum()
	    absres_sum[iter] = abs(ma.masked_array(res, mask)[n.where(R <= extraction_radius)]).sum()
#	    print absres_sum[iter], absres_sum1
            if(flag_image):
                Aabs[iter] = absres_sum[iter] / (abs_zsum[iter])
            else:
                Aabs[iter] = absres_sum[iter] / (ABS_ZSUM)
        Aabso = Aabs.min()
        index = Aabs.argmin()
        if(nn > 20):
            xcntr[index] = ini_xcntr
            ycntr[index] = ini_ycntr
        if(xcntr[index] == ini_xcntr and ycntr[index] == ini_ycntr):
            ff = pyfits.open("MaskedGalaxy.fits")
	    ZZZ = ff[0].data
	    ff.close()
	    RoTz = rotate_modi(ZZZ, 180.0, xcntr[iter], ycntr[iter], prefilter=False)
	    os.system('rm -f Rotated.fits')
	    hdu = pyfits.PrimaryHDU(RoTz)
	    hdu.writeto('Rotated.fits')
            f = pyfits.open("Rotated.fits")
            rz = f[0].data
            f.close()
            rz = n.swapaxes(rz, 0, 1)
            res = z - rz
            hdu = pyfits.PrimaryHDU(n.swapaxes(res, 0, 1).astype(n.float32))
            hdu.writeto('AResidual.fits')
            flag_center = 1
            rot_sum = rot_sum[index]
            abs_zsum = abs_zsum[index]
            if(flag_image == 0):
                abs_zsum = ABS_ZSUM
            sh_sum = sh_sum[index]
            absres_sum = absres_sum[index]
            if(flag_image):
                error_asym = n.sqrt( Aabso**2*( ((sh_sum + rot_sum + 4 * \
                             background * \
                             R[n.where(R <= n.floor(extraction_radius))].size) \
                             / absres_sum**2.0) + ((sh_sum + 2.0 * background *\
                             R[n.where(R <= n.floor(extraction_radius))].size) \
                             / abs_zsum**2) ) )
            else:
                error_asym = n.sqrt( Aabso**2 * (((sh_sum+rot_sum + 2 * \
                             background * \
                             R[n.where(R <= n.floor(extraction_radius))].size) \
                             / absres_sum**2.0) + ((sh_sum + 2.0 * background \
                           * R[n.where(R <= n.floor(extraction_radius))].size) \
                             / abs_zsum**2) ) )
        else:
            ini_xcntr = xcntr[index]
            ini_ycntr = ycntr[index]
            nn += 1
    for myfile in ['MaskedGalaxy.fits', 'MRotated.fits']:
        if os.access(myfile, os.F_OK):
            os.remove(myfile)
    return Aabso, error_asym, ini_xcntr, ini_ycntr, nn, flag_out,\
           abs_zsum, sh_sum, extraction_radius

#asymmetry('n5585_lR.fits','BMask.fits', 192.03,157.42,0.0,0.0, 100.0, 1390.377,180.0, 0, 1000000)
