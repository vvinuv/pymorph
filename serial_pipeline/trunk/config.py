"""Configure file for automation"""
imagefile = 'j8f645-1-1_drz_sci.fits'
whtfile = 'j8f645-1-1_drz_rms.fits'   #The weight image. 
sex_cata = 'j8f645_sex.cat'           #The sextractor catalogue which has 
                                      #the format given in the file
clus_cata = 'cl1216-1201.cat'         #catalogue of galaxies from
                                      #online catalogu service
                                      #(name ra1 ra2 ra2 dec1 dec2 dec3)
out_cata = 'cl1216-1201_out.cat'      #catalogue of galaxies in the field

psflist = ('psf_1216435-1203120.fits', 'psf_1216392-1202402.fits', 'psf_1216382-1200443.fits', 'psf_1216433-1202030.fits', 'psf_1216446-1203464.fits', 'psf_1216384-1201179.fits')                    #List of psf containg their 
                                      #position information in the 
                                      #header (RA_TARG, DEC_TARG). 
                                      #Make psf with the names as here 
                                      #and use psf_header_update.py. 
                                      #It will update the header information.
mag_zero = 25.256                     #magnitude zero point
mask_reg = 2.0
thresh_area = 1650.0
threshold = 2.5                       #Masking will start for neighbours 
                                      #whose distace from the object greater 
                                      #than threshold * semi-major axis of 
                                      #the object and area of the neighbour 
                                      #less than thresh_area sq.pixel. 
                                      #The masking will be for a circular 
                                      #region of radius mask_reg*semi-major 
                                      #axis of the nighbour with respect to 
                                      #the center of the neightbour.
size = 120                            #size of the stamp image
shiftra = 0.000212500000004 
shiftdec =  6.66666666653e-05         #If the image WCS is not same as the 
                                      #coordinate given in the clus_cata, 
                                      #the appropriateshiftra and shiftdec  
                                      #should be used. It will be better to 
                                      #correct WCS using iraf command ccmap 
                                      #so that the programcan identify the 
                                      #correct objects. Remember: Shift 
                                      #in the frame in not LINEAR! and it 
                                      #can leads to detect wrong objects
pixelscale = 0.045                    #Pixel scale (arcsec/pixel)

H0 = 71                               #Hubble parameter
WM = 0.27                             #Omega matter
WV = 0.73                             #Omega Lambda
repeat = False                        #Repeat the pipeline manually
