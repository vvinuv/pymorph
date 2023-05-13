PyMorph
=======
This pipeline was written as a part of some projects on galaxy morphology by Vinu Vikraman and Alan Meert. Yogesh Wadadekar, Ajit Kembhavi, G V Vijayagovidan and  Mariangela Bernardi were working on those projects. It REQUIRES SExtractor and GALFIT.  

### INSTALLATION
1. Download pymorph either by zip or by

```
https://github.com/vvinuv/pymorph.git
```

2. In the pymorph folder 

```
pip install .
``` 

### BUGS
Current version is buggy and not recommended to use this version.

### If you want to really use this

```
>>> import pymorph
>>> p = pymorph.PyMorph()
>>> p.pymorph()
```

This assumes that there is a configuration file in your directory. It runs also with default run parameters which are mostly different from the parameters in the configuration file. You can of course change this. Keep on read this


PyMorph is a pipeline which integrates SExtractor and GALFIT for estimating parametric and nonparametric quantities to describe the galaxy morphology. It can be used on a large frames with or without object coordinates of target galaxies or individual target objects. It can use PSF which is nearest to the target galaxy based on SExtractor star-galaxy probability or list of PSFs or individual PSF. It also gives the abilities to select the PSF interactively from a given PSF list or using SExtractor. It can repeat the selection procedure until you are satisfied. It can use SExtractor sky, GALFIT sky or the sky value implemented in PyMorph. It can repeat the fitting if GALFIT fails. It uses different initial parameters for GALFIT during refitting. It calculates nonparameteric quantities such as concentration, asymmetry, clumpness, Gini coefficient and moment parameters. These are implimented in PyMorph. It creates plots and a simple html page. It can write database for the output parameters.  

PyMorph assumes a configuration file in the current directory which is shown in the below

```
[imagecata]
;Specify the input images and Catalogues
imagefile = withsky.fits

;The weight image. If it contains the string 'rms', this will treated as 
;RMS_MAP and if it contains weight, then that will be treated as WEIGHT_MAP. 
;If nothing found, then by default it is treated as MAP_RMS 
whtfile = rms.fits   

;The sextractor catalogue which has the format given in the file
;catalogue of galaxies from online catalogu service
;(name ra1 ra2 ra2 dec1 dec2 dec3)
sex_cata = sex.cat            
clus_cata = input.cat         


;Specify the output names of images and catalogues
;catalogue of galaxies in the field
out_cata = out.cat      
rootname = cl1358

;the directory containing input images if commented out, then program uses
;current directory
datadir = /Users/vinu/github/pymorph_refactoring/examples/small_image/data/ 

;the directory containing output data if commented out, then program uses
;current directory
outdir = /Users/vinu/github/pymorph_refactoring/examples/small_image/results/

[external]
;;;----Set the SExtractor and GALFIT path here----;;;
GALFIT_PATH = /usr/local/bin/galfit
;/home/vinu/software/galfit/modified/galfit 
SEX_PATH = /usr/local/bin/sex
PYMORPH_PATH = /Users/vinu/github/pymorph_refactoring/

[psf]
;Star-galaxy classification 
stargal_prob = 0.8 

;Psf list
;0 => No psfselection
;1 => Only Select psf 
;2 => Select psf and run pipeline
;Recommended: Run with '1' and then run pipeline
psfselect = 0                         
                                      
;psf image size will be startsize times the SMA given by SExtractor
star_size = 20                         

;List of psf containg their position information in the 
;header (RA_TARG, DEC_TARG). Make psf with the names as here 
;and use psf_header_update.py. It will update the header information.

psflist = psf_1447219+0828392.fits
;[psf_1216382-1200443.fits, psf_1216408-1200251.fits, psf_1216424-1202057.fits,psf_1216487-1201246.fits,psf_1216504-1202104.fits]   
;psflist = @psflist.list

;0 nearest, 1 second nearest, 2 third nearest etc
which_psf = 0

#area of object for psf selection
area_obj = 40

[mask]
;;;----Conditions for Masking----;;;
;Masking will be done for neighbours whose semimajor*threshold overlaps with 
;threshold * semi-major axis of the object and area of the neighbour 
;less than thresh_area * object area in ;sq.pixel. 
;The masking will be for a circular region of radius mask_reg*semi-major 
;axis of the nighbour with respect to the center of the neightbour.
manual_mask = 0
mask_reg = 2.0
thresh_area = 0.2
threshold = 3.0                       
                                      
[size]
;Size of the cut out and search conditions
;size = [resize?, varsize?, fracrad, square?, fixsize]
;size of the stamp image
size_list = 0, 1, 6, 1, 120              
;The search radius 
searchrad = 0.3arc                     

[cosmology]
;Parameters for calculating the physical parameters of galax
;Pixel scale (arcsec/pixel)
pixelscale = 0.045                    
;Hubble parameter
H0 = 71                               
;Omega matter
WM = 0.27                             
;Omega Lambda
WV = 0.73                             
;redshift
redshift = 9999

[nonparams]
;Parameters to be set for calculating the CASGM----;;;
back_extraction_radius = 15.0
;back_ini_xcntr = 32.0 
;back_ini_ycntr = 22.0
angle = 180.0

[modes]
;PyMorph fitting modes
;Repeat the pipeline manually
repeat = False                        
;True if we provide cutouts
galcut = False                        
decompose = True
;Detailed fitting
detail = False 
;Always keep this True as it is not functional yet!
galfit = True 
cas = False
;if findandfit= True, then maglim = [faint_mag, bright_mag]
;findandfit takes all the objects in clus_cata and fit it in one go.
;Need to be very careful when setting the sextractor parameters
;images start with I and weight start with W
findandfit = True
crashhandler = False

[mag]
;magnitude zero point
mag_zero = 25.256                     
maglim = 22, 15 

[galfit]
;Galfit Controls
;The components to be fitted to the object
components = bulge, disk, bar, point
;set to False to fit sersic bulge, set to true to fit devacouler's bulge (n = 4)
devauc = False 

;;;---fixing = [bulge_center, disk_center, sky, bar_center, point_center]
; = 0, Fix params at SExtractor value
fitting = 1, 1, 0, 1, 0

[version]
galfitv = 3.0.2
pymorphv = 0.1

[diagnosis]
;;;----The following conditions are used to classify fit goo/bad----;;;
chi2sq = 1.9                          
Goodness = 0.60                       
;< abs(center - fitted center)
center_deviation = 3.0                
;Keep center within +/- center_constrain
center_constrain = 2.0                

[db]
;;;----Database Informations----;;;
host = localhost
database = Cluster
table = cluster
usr = vinu
pword = cluster
dbparams = Cluster:cl1216-1201, ObsID:1:int
```

We can use the following 'command line' like arguments to change the default parameters

```
-e, --sex_conf: Edit SExtractor configuration (default = False)
-f, --force: Remove already existing SExtractor catalog (default = False)
-t, --test: Runs the test instance-OVERRIDES ALL OTHER INPUT-User must supply a directory for output (default = False)
--lmag: Lower magnitude cutoff (default = 500)
--umag: Upper magnitude cutoff" (default = -500)
--ln: Lower Sersic (default = 0.1)
--un: Upper Sersic (default = 500)
--lre: Lower Bulge Radius (default = 0)
--ure: Upper Bulge Radius (default = 500)
--lrd: Lower Disk Radius (default = 0)
--urd: Upper Disk Radius (default = 500)

--bdbox: Turns on bdbox (default = False)
--bbox: Turns on bbox (default = False)
--dbox: Turns on dbox (default = False)
--devauc: Turns on DeVacouleur's bulge fitting (dafault: False)


--with-in: Remove boundary of images (default = 50)
--with-filter: Filter used (default = 'UNKNOWN')
                  

-p, --with-psf: Nearest/farthest PSF (default = 0)
--with-sg: For psf identification (default = 0.9)
--with-area: Min area of psf for selection (default = 40)
--no_mask: Turn off masking (default=False)
--norm_mask: Turns on Normal masking (default=False)

--with-host: MySql host used (default = 'localhost')
--with-db: Database used (default = 'UNKNOWN')
```


