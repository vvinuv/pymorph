from numarray import *
import numarray.nd_image as image
import numarray.ma as ma
import pyfits
from pyraf import iraf

class asymmetry:
    """asymmetry parameter"""		
    def __init__(self, cutimage, maskimage, ini_xcntr, ini_ycntr, pa, eg, extraction_radius, background, angle, flag_image, ABS_ZSUM):
        self.cutimage		= cutimage
        self.maskimage          = maskimage
        self.ini_xcntr		= ini_xcntr 
        self.ini_ycntr		= ini_ycntr
        self.extraction_radius	= extraction_radius
        self.pa			= pa
        self.eg			= eg
        self.background		= background
        self.one_minus_eg_sq	= (1.0-eg)**2.0
        self.flag_image		= flag_image
        self.angle		= angle
        self.ABS_ZSUM		= ABS_ZSUM
#The flag_image=1 for image asymmetry and flag_image=0 for background
        self.image_asymm=ASYM(self.cutimage, self.maskimage, self.ini_xcntr, self.ini_ycntr, self.pa, self.one_minus_eg_sq, self.background, self.extraction_radius, self.angle, self.flag_image, self.ABS_ZSUM)
        return 

def ASYM(cutimage, maskimage, ini_xcntr, ini_ycntr, pa, one_minus_eg_sq, background, extraction_radius, angle, flag_image, ABS_ZSUM):
    iraf.images(_doprint=0)
    iraf.images.imgeom(_doprint=0)
    iraf.images.imutil(_doprint=0)
    co= math.cos(pa*pi/180.0)
    si= math.sin(pa*pi/180.0)
    Aabs=zeros([9],type=Float32)
    rot_sum=zeros([9],type=Float32)
    abs_zsum=zeros([9],type=Float32)
    sh_sum=zeros([9],type=Float32)
    absres_sum=zeros([9],type=Float32)
    f=pyfits.open(cutimage)
    z = f[0].data
    f.close()
    iraf.rotate("".join(maskimage), output="MRotated.fits",\
                        rotation=180, xin = ini_xcntr + 1.0,\
                        yin = ini_ycntr + 1.0, \
                        xout = ini_xcntr + 1.0, yout = ini_ycntr + 1.0, \
                        ncols= "INDEF",\
                        nlines = "INDEF", interpo= "linear", \
                        boundar="nearest", verbose="no")
    iraf.imarith("".join(maskimage), op="*", operand2="MRotated.fits", \
                    result="AMask.fits")
    f=pyfits.open("AMask.fits")
    mask = f[0].data
    f.close()
    maskedgalaxy = ma.masked_array(z, mask)
    z = ma.filled(maskedgalaxy, value = background)
    hdu = pyfits.PrimaryHDU(z.astype(Float32))
    hdu.writeto('MaskedGalaxy.fits')
    NXPTS = z.shape[0]
    NYPTS = z.shape[1]
    #r50=c.r50
    center_rad=0.2#0.01*r50 The centering correction has to be done by calculating asymmetric parameters at different points 0.1% of r50 away around the center
    flag_center=0 #flag_center=1 if the program finds the exact center else center=0
	#extraction_radius=total_rad ie. 1.5 times the radius at eta=0.2
    x = reshape(arange(NXPTS * NYPTS), (NXPTS, NYPTS)) % NYPTS
    x = x.astype(Float32)
    y = reshape(arange(NXPTS * NYPTS), (NXPTS, NYPTS)) / NYPTS
    y = y.astype(Float32)
    n = 0 #This n is given here and the following loop run either it finds the minimum asymmetry or n=20
    while flag_center == 0:
        flag_out = 0
		#x and y coordinates of different center points in 1-d array
        xcntr = reshape(array([[ini_xcntr]*3,[ini_xcntr-center_rad]*3,\
                [ini_xcntr+center_rad]*3],type=Float32),(9,)) 
        ycntr = array([ini_ycntr,ini_ycntr-center_rad,ini_ycntr+center_rad]*3,\
              type = Float32)
        for iter in range(9):#The below is done in order to keep the extraction radius inside the image. If the extraction radius goes outside the image then flag_out=1
            if(flag_image):#flag_image will tell whether the input is image (flag_image=1) or background (flag_image=0)
                if(xcntr[iter]-2-extraction_radius<0):
                    extraction_radius=xcntr[iter]-2
                    flag_out=1
                if(ycntr[iter]-2-extraction_radius<0):
                    extraction_radius=ycntr[iter]-2
                    flag_out=1
                if(xcntr[iter]+2+extraction_radius>NXPTS):
                    extraction_radius=NXPTS-xcntr[iter]-2
                    flag_out=1
                if(ycntr[iter]+2+extraction_radius>NYPTS):
                    extraction_radius=NYPTS-ycntr[iter]-2
                    flag_out=1
            tx = (x-xcntr[iter])*co + (y-ycntr[iter])*si
            ty = (xcntr[iter]-x)*si + (y-ycntr[iter])*co
            R=sqrt(tx**2.0+ty**2.0/one_minus_eg_sq)
            iraf.rotate("".join("MaskedGalaxy.fits"), output="Rotated.fits",\
                        rotation=180, xin = xcntr[iter] + 1.0, \
                        yin = ycntr[iter] + 1.0, \
                        xout = xcntr[iter] + 1.0, yout = ycntr[iter] + 1.0, \
                        ncols= "INDEF",\
                        nlines = "INDEF", interpo= "linear", \
                        boundar="nearest", verbose="no")
            f=pyfits.open("Rotated.fits")
            rz = f[0].data
            f.close()
            res = z - rz
            for myfile in ['Rotated.fits','Residual.fits']:
                if os.access(myfile,os.F_OK):
                    os.remove(myfile)
            sh_sum[iter]=z[where(R<=extraction_radius)].sum()
            rot_sum[iter]=rz[where(R<=extraction_radius)].sum()
            abs_zsum[iter]=abs(z[where(R<=extraction_radius)]).sum()	
            absres_sum[iter]=abs(res[where(R<=extraction_radius)]).sum()
            if(flag_image):
                Aabs[iter]=absres_sum[iter]/(abs_zsum[iter])
            else:
                Aabs[iter]=absres_sum[iter]/(ABS_ZSUM)
        Aabso=Aabs.min()
        index=Aabs.argmin()
        if(n>20):
            xcntr[index] = ini_xcntr
            ycntr[index] = ini_ycntr
        if(xcntr[index] == ini_xcntr and ycntr[index] == ini_ycntr):
            flag_center=1
            rot_sum=rot_sum[index]
            abs_zsum=abs_zsum[index]
            if(flag_image==0):
                abs_zsum=ABS_ZSUM
            sh_sum=sh_sum[index]
            absres_sum=absres_sum[index]
            if(flag_image):
                error_asym=sqrt( Aabso**2*( ((sh_sum+rot_sum+4*background*\
                           R[where(R<=floor(extraction_radius))].nelements())/\
                           absres_sum**2.0)+((sh_sum+2.0*background*\
                           R[where(R<=floor(extraction_radius))].nelements())/\
                           abs_zsum**2) ) )
            else:
                error_asym=sqrt( Aabso**2*( ((sh_sum+rot_sum+2*background*\
                           R[where(R<=floor(extraction_radius))].nelements())/\
                           absres_sum**2.0)+((sh_sum+2.0*background*\
                           R[where(R<=floor(extraction_radius))].nelements())/\
                           abs_zsum**2) ) )
        else:
            ini_xcntr=xcntr[index]
            ini_ycntr=ycntr[index]
            n+=1
    for myfile in ['AMask.fits', 'MaskedGalaxy.fits']:
        if os.access(myfile, os.F_OK):
            os.remove(myfile)
    return Aabso,error_asym,ini_xcntr, ini_ycntr,n, flag_out,\
           abs_zsum,sh_sum,extraction_radius

