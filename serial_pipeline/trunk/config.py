"""Configure file for automation"""
###----Specify the input images and Catalogues----###
imagefile = 'j8f647-1-1_drz_sci.fits'
whtfile = 'j8f647-1-1_drz_rms.fits'   #The weight image. 
sex_cata = 'j8f647_sex.cat'           #The sextractor catalogue which has 
                                      #the format given in the file
clus_cata = 'cl1216-1201.cat'         #catalogue of galaxies from
                                      #online catalogu service
                                      #(name ra1 ra2 ra2 dec1 dec2 dec3)

###----Specify the output names of images and catalogues----###
out_cata = 'cl1216-1201_out.cat'      #catalogue of galaxies in the field
rootname = 'j8f647'

###----Psf list----###
psflist = ('psf_1216382-1200443.fits', 'psf_1216408-1200251.fits', 'psf_1216424-1202057.fits','psf_1216487-1201246.fits','psf_1216504-1202104.fits')   
                                      #List of psf containg their 
                                      #position information in the 
                                      #header (RA_TARG, DEC_TARG). 
                                      #Make psf with the names as here 
                                      #and use psf_header_update.py. 
                                      #It will update the header information.
mag_zero = 25.256                     #magnitude zero point

###----Conditions for Masking----###
manual_mask = 0
mask_reg = 2.0
thresh_area = 0.2
threshold = 3.0                       #Masking will be done for neighbours 
                                      #whose semimajor*threshold overlaps with 
                                      #threshold * semi-major axis of 
                                      #the object and area of the neighbour 
                                      #less than thresh_area * object area in
                                      #sq.pixel. 
                                      #The masking will be for a circular 
                                      #region of radius mask_reg*semi-major 
                                      #axis of the nighbour with respect to 
                                      #the center of the neightbour.

###Size of the cut out and search conditions---###
size = 120                            #size of the stamp image
shiftra = 0.0 
shiftdec =  0.0                       #If the image WCS is not same as the 
                                      #coordinate given in the clus_cata, 
                                      #the appropriateshiftra and shiftdec  
                                      #should be used. It will be better to 
                                      #correct WCS using iraf command ccmap 
                                      #so that the programcan identify the 
                                      #correct objects. Remember: Shift 
                                      #in the frame in not LINEAR! and it 
                                      #can leads to detect wrong objects

###----Parameters for calculating the physical parameters of galaxy----###
pixelscale = 0.045                    #Pixel scale (arcsec/pixel)
H0 = 71                               #Hubble parameter
WM = 0.27                             #Omega matter
WV = 0.73                             #Omega Lambda

###----Parameters to be set for calculating the CASGM----###
back_extraction_radius = 15.0
#back_ini_xcntr = 32.0 
#back_ini_ycntr = 22.0
angle = 180.0

###----Fitting modes----###
repeat = False                        #Repeat the pipeline manually
galcut = False                        #True if we provide cutouts
decompose = True
galfit = True #Always keep this True as it is not functional yet!
cas = True

###----Set the SExtractor and GALFIT path here----###
GALFIT_PATH = '/home/vinu/software/galfit/modified/galfit' 
SEX_PATH = '/home/vinu/software/sextractor-2.5.0/sex/bin/sex'
PYMORPH_PATH = '/home/vinu/factory/cl1216-1201/j8f647/pymorph'

###----The following conditions are used to classify fit goo/bad----###
chi2sq = 1.9
Goodness = 0.80
center_deviation = 5.0
