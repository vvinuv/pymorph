import os
import numpy as np	
import numpy.ma as ma
import pyfits
from rotate import rotate, ImSec
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
        self.image_asymm = ASYM(self.cutimage, self.maskimage, \
                                self.ini_xcntr, self.ini_ycntr, self.pa, \
                                self.one_minus_eg_sq, r50, self.background, \
                                self.extraction_radius, self.angle, \
                                self.flag_image, self.ABS_ZSUM)
        return 

def ASYM(cutimage, maskimage, ini_xcntr, ini_ycntr, pa, one_minus_eg_sq, r50, background, extraction_radius, angle, flag_image, ABS_ZSUM):
    if flag_image == 0:
	maskimage = 'BackMask.fits'
    co = np.cos(pa * np.pi / 180.0)
    si = np.sin(pa * np.pi / 180.0)
    Aabs = np.zeros(9)
    Aabs = Aabs.astype(np.float32)
    rot_sum = np.zeros(9)
    rot_sum = rot_sum.astype(np.float32)
    abs_zsum = np.zeros(9)
    abs_zsum = abs_zsum.astype(np.float32)
    sh_sum = np.zeros(9)
    sh_sum = sh_sum.astype(np.float32)
    absres_sum = np.zeros(9)
    absres_sum = absres_sum.astype(np.float32)
    f = pyfits.open(c.datadir + cutimage)
    z = f[0].data
    z = z - background
    f.close()
    fm = pyfits.open(c.outdir + maskimage)
    zm = fm[0].data
    fm.close()
    # Rotate the mask
    rzm = rotate(zm, 180.0, ini_xcntr, ini_ycntr)
    mask = zm + rzm # Adding rotated and original mask
    center_rad = 0.01 * r50 # 0.01 * r50 The centering correction has to 
                            # be done by calculating asymmetric parameters
                            # at different points 0.1% of r50 away around 
                            #the center

    flag_center = 0 # flag_center=1 if the program finds the 
                    # exact center else center=0
    nn = 0 # This nn is given here and the following loop run either 
           # it finds the minimum asymmetry or n=20
    while flag_center == 0:
        flag_out = 0
        # x and y coordinates of different center points in 1-d array
        xcntr = np.reshape(np.array([[ini_xcntr] * 3, \
                          [ini_xcntr - center_rad] * 3, \
                          [ini_xcntr + center_rad] * 3]), (9, )) 
        xcntr = xcntr.astype(np.float32)
        ycntr = np.array([ini_ycntr, ini_ycntr - center_rad, \
                          ini_ycntr + center_rad]*3)
        ycntr = ycntr.astype(np.float32)
        for iter in range(9): # The below is done in order to keep the 
                              # extraction radius inside the image. If the 
                              # extraction radius goes outside the image 
                              # then flag_out=1
            CutImDa, cut_xcntr, cut_ycntr, SizeY, SizeX, ymin, \
           ymax, xmin, xmax, flag_out = \
                       ImSec(z, xcntr[iter], ycntr[iter], extraction_radius)
            if flag_image == 0: # flag_image will tell whether the input is \
                                # image (flag_image=1) or background 
                flag_out = 0
            x = np.reshape(np.arange(SizeX * SizeY), (SizeY, SizeX)) % SizeX
            x = x.astype(np.float32)
            y = np.reshape(np.arange(SizeX * SizeY), (SizeY, SizeX)) / SizeX
            y = y.astype(np.float32)
            tx = (x - cut_xcntr) * co + (y - cut_ycntr) * si
            ty = (cut_xcntr - x) * si + (y - cut_ycntr) * co
            R = np.sqrt(tx**2.0 + ty**2.0 / one_minus_eg_sq)
            # print xcntr[iter] , ycntr[iter]
	    rz = rotate(CutImDa, 180.0, cut_xcntr, cut_ycntr)
	    res = CutImDa - rz
            masksub = mask[ymin:ymax, xmin:xmax].copy()

            # If want to check the rotated image and cutout image
            # uncomment the following
            # for myfile in ['Rotated.fits', 'CutImDa.fits']:
            #    if os.access(myfile, os.F_OK):
            #        os.remove(myfile)
	    # hdu = pyfits.PrimaryHDU(rz)
            # hdu.writeto('Rotated.fits')
	    # hdu = pyfits.PrimaryHDU(CutImDa)
            # hdu.writeto('CutImDa.fits')
            # raw_input()
            # END

	    sh_sum[iter] = ma.masked_array(CutImDa, masksub) \
                           [np.where(R <= extraction_radius)].sum()
	    rot_sum[iter] = ma.masked_array(rz, masksub) \
                            [np.where(R <= extraction_radius)].sum()
            abs_zsum[iter] = abs(ma.masked_array(CutImDa, masksub) \
                             [np.where(R <= extraction_radius)]).sum()
	    absres_sum[iter] = abs(ma.masked_array(res, masksub) \
                               [np.where(R <= extraction_radius)]).sum()
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
            flag_center = 1
	    rz = rotate(CutImDa, 180.0, cut_xcntr, cut_ycntr)
            res = CutImDa - rz
            rot_sum = rot_sum[index]
            abs_zsum = abs_zsum[index]
            if flag_image == 0:
                abs_zsum = ABS_ZSUM
            sh_sum = sh_sum[index]
            absres_sum = absres_sum[index]
            if(flag_image):
                error_asym = np.sqrt( Aabso**2*( ((sh_sum + rot_sum + 4 * \
                          background * \
                          R[np.where(R <= np.floor(extraction_radius))].size) \
                          / absres_sum**2.0) + ((sh_sum + 2.0 * background * \
                          R[np.where(R <= np.floor(extraction_radius))].size) \
                          / abs_zsum**2) ) )
            else:
                error_asym = np.sqrt( Aabso**2 * (((sh_sum+rot_sum + 2 * \
                          background * \
                          R[np.where(R <= np.floor(extraction_radius))].size) \
                          / absres_sum**2.0) + ((sh_sum + 2.0 * background * \
                          R[np.where(R <= np.floor(extraction_radius))].size) \
                          / abs_zsum**2) ) )
        else:
            ini_xcntr = xcntr[index]
            ini_ycntr = ycntr[index]
            nn += 1
    # print Aabso, error_asym, ini_xcntr, ini_ycntr, nn, flag_out,\
    #       abs_zsum, sh_sum, extraction_radius
    return Aabso, error_asym, ini_xcntr, ini_ycntr, nn, flag_out,\
           abs_zsum, sh_sum, extraction_radius

#Aabso, error_asym, ini_xcntr, ini_ycntr, nn, flag_out, abs_zsum, sh_sum, extraction_radius = asymmetry('n5585_lR.fits','BMask.fits', 192.03,157.42,0.0,0.0, 50, 100.0, 1390.377,180.0, 1, 1000000)
#print Aabso, error_asym, ini_xcntr, ini_ycntr, nn, flag_out, abs_zsum, sh_sum, extraction_radius
