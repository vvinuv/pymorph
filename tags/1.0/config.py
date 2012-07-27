"""Configure file for automation"""
imagefile = 'j8f645-1-1_drz_sci.fits'
whtfile = 'j8f645-1-1_drz_rms.fits' #The weight from multidrizzle. The program will calculate the weight to galfit as 1/ sqrt(whtimage)
sex_cata = 'j8f645_sex.cat'
clus_cata = 'cl1216-1201.cat' #catalogue of galaxies from online catalogu service
out_cata = 'cl1216-1201_out.cat' ##catalogue of galaxies in the field
psflist = ('psf_1216435-1203120.fits', 'psf_1216392-1202402.fits', 'psf_1216382-1200443.fits', 'psf_1216433-1202030.fits', 'psf_1216446-1203464.fits', 'psf_1216384-1201179.fits')
mag_zero = 25.256 #magnitude zero point
mask_reg = 2.0
thresh_area = 1650.0
threshold = 2.5 #Masking will start for neighbours whose distace from the object greater than threshold * semi-major axis of the object and area of the neighbour less than thresh_area sq.pixel. The masking will be for a circular region of radius mask_reg*semi-major axis of the nighbour with respect to the center of the neightbour.
size = 120 # size of the stamp image
shiftra = 0.000212500000004 
shiftdec =  6.66666666653e-05
