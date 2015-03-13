# This is a wrapper script for galfit.pl.  ACS images are huge, and it takes
# more time for galfit to read in an image than it takes to fit a galaxy.
# Thus, this script takes an ACS image and splits it up into little chunks.
# Then, it outputs a temporary "perl.in" file which is read in by the
# "galfit.pl" Perl script.  The Perl script reads in a Sextractor file
# for galaxy information, formats it into galfit input file, figures out
# which galaxy falls into which panel, then finally runs galfit.
#
# Once all the galaxies have been fit, we return to this IRAF script
# This script then reformats the headers of all the postage stamp images to
# reflect the correct information relative to the big ACS image.
#
# Parameter explanations:
#   skyval - The value of the sky background in the image being fitted that
#            SExtractor knows nothing about.  This value gets added to 
#            background determined by SExtractor when creating the galfit 
#            input file, and doesn't actually change the image being fitted
#            or the catalog itself.  This is needed because sometimes people
#            run Sextractor on images that have the background removed, only
#            to realize they want GALFIT to create a sigma map internally,
#            which needs the sky background intact.  If this is the case,
#            the user has to add the sky background back into the image
#            manually.

procedure galfit (image, sigimg, sexcat, startnum, stopnum)

string image {prompt = "Input image to be fitted"}
string sigimg {prompt = "Sigma image to use"}
string sexcat {prompt = "Sextractor catalog name"}
int    startnum {1, prompt = "Start from which object number in the catalog?"}
int    stopnum {999999, prompt = "Stop at which number in the catalog?"}
string psfloc {prompt = "PSF location (directory)"}
string constraint {prompt = "Parameter constraint file"}
string fitfunc {enum='sersic|BD', prompt = "Single sersic or B/D decomposition?"}
string region {"A", prompt = "(A) Whole image catalog or, (B) region of image?"}
string boundary {"0 0 0 0", prompt = "Boundary of region if (B) (xmin xmax ymin ymax)"}
real   skyval {"0.", prompt = "Sky background to add back to Sextractor value"}
string mrange {"0. 99.", prompt = "Acceptable magnitude range"}
string stargal {"0. 1.", prompt = "Range of good galaxies: 0.0 (gal) to 1.0 (star)"}
string fwhmpsf {"0.06", prompt = "FWHM of the PSF for cosmic ray rejection"}
string magzpt {"22.", prompt = "Photometric zeropoint"}
string plate {"0.03", prompt = "Image plate scale"}
bool   exist {no, prompt= "Rid all existing temp-*.fits images before run?"}
bool   nofat {yes, prompt= "Rid all temp images, galfit input, after run?"}
real nsplit {10., prompt = "Split image into how many parts along one axis?"}

struct *fileptr

begin

string img, sex, stampname, section, errstr, param, tempout, dum, sig, sigout
int nrow, ncol, irow, icol, nsubx, nsuby, pos, pos2, objn, snum, stnum
int BUFFER, xlo, xhi, ylo, yhi, offx, offy, i, last
real val, xcboxcent, ycboxcent
struct line

BUFFER = 200

img = image
sig = sigimg
sex = sexcat
snum = startnum
stnum =  stopnum

if (!access ("galfit-tempfits"))
    mkdir ("galfit-tempfits")

print ("")
print ("Putting all sub-panel images into the directory: galfit-tempfits.")
print ("")

################  Split image into nsplit panels along a side  ################

imgets (img, "naxis1")
ncol = int (imgets.value)
imgets (img, "naxis2")
nrow = int (imgets.value)

nsubx = int (ncol / nsplit)
nsuby = int (nrow / nsplit)

# Get rid of the last / in the directory name:

last = strlen (psfloc)
dum = substr (psfloc, last, last)
if (dum == "/")
    psfloc = substr(psfloc, 1, last-1)

# Check to see if we should first delete existing sub-panels

if (exist) {
    imdel ("galfit-tempfits/temp-*.fits")
    imdel ("galfit-tempfits/sig-*.fits")
}

for (j=1; j <= nsplit; j=j+1) {
    for (i=1; i <= nsplit; i=i+1) {
        if (nsubx * (i-1) - BUFFER <= 1) {
	        xlo = 1
	} else {
	    xlo = nsubx * (i-1) - BUFFER + 1
        }

	if (nsuby * (j-1) - BUFFER <= 1) {
	    ylo = 1
	} else {
	    ylo = nsuby * (j-1) - BUFFER + 1
        }

	if (nsubx * i + BUFFER >= ncol) {
	    xhi = ncol
	} else {
	    xhi = nsubx * i + BUFFER
        }

	if (nsuby * j + BUFFER >= nrow) {
	    yhi = nrow
	} else {
	    yhi = nsuby * j + BUFFER
        }

	tempout = "galfit-tempfits/temp-"//i//"-"//j//".fits"
	sigout = "galfit-tempfits/sig-"//i//"-"//j//".fits"
	if (access (img) && !access (tempout))
	    imcopy (in=img//"["//xlo//":"//xhi//","//ylo//":"//yhi//"]", 
								out=tempout)
	if (access (sig) && !access (sigout))
	    imcopy (in=sig//"["//xlo//":"//xhi//","//ylo//":"//yhi//"]", 
								out=sigout)

    }
}

#######################  Set up and run Perl script  ##########################

if ( access ( "perl.in" ) ) delete ( "perl.in" )

print (img, "	# The ACS image name to fit", >> "perl.in")
print (sex, "	# The sextractor catalog", >> "perl.in")
print (psfloc, "		# PSF directory", >> "perl.in")
print (snum, "	# Starting from object number", >> "perl.in")
print (stnum, "	# Stop at object number", >> "perl.in")
print (constraint, "		# Parameter constraint file", >> "perl.in")
print ("galfit-tempfits		# Temporary working directory", >> "perl.in")
print (fitfunc, "		# Fit a single sersic or do B/D decomposition?", >> "perl.in")
print (region, "		# (A) Whole image catalog or, (B) region of the image", >> "perl.in")
print (boundary, "		# Boundary of the object region if (B)", >> "perl.in")
print (skyval, "		# Sky value to add to Sextractor determined value", >> "perl.in")
print (mrange, "		# Acceptable magnitude range", >> "perl.in")
print (stargal, "		# Range of good galaxies: 0.0 (galaxy) to 1.0 (star)", >> "perl.in")
print (fwhmpsf, "		# FWHM of the PSF for cosmic ray rejection", >> "perl.in")
print (magzpt, "		# Photometric zeropoint", >> "perl.in")
print (plate, "		# Plate scale", >> "perl.in")
print (ncol, nrow, "	# ncol, nrow of original ACS image", >> "perl.in")
print ("temp		# Temporary sub-panel file names", >> "perl.in");
print (nsplit, "		# Number of subpanels along a side to split ACS image", >> "perl.in");
print (BUFFER, "		# Annulus width to pad around the split image", >> "perl.in");


if (access ("psf.temp"))
    delete ("psf.temp")

dir (psfloc//"/psf*.fits", ncols=1, maxch = 80, > "psf.temp")

######### Enter into Perl script, which runs GALFIT

! perl galfit.pl perl.in

######### After GALFIT has been run on ALL the galaxies in the database 
######### modify the header info of the image blocks.

fileptr = "galfit-tempfits/offsets"

while (fscan (fileptr, objn, stampname, offx, offy, section, xcboxcent, ycboxcent) != EOF) {
    if (access (stampname//".fits")) {
	hedit (ima=stampname//"[2]", fie="objnum", val=objn, add+, ver-, update+)
        hedit (ima=stampname//"[1]", fie="object", val=img//section, ver-, update+)
	hedit (ima=stampname//"[2]", fie="datain", val=img, ver-, update+)
	hedit (ima=stampname//"[2]", fie="noise", val=sig, ver-, update+)
        hedit (ima=stampname//"[2]", fie="fitsect", val=section, ver-, update+)
        hedit (ima=stampname//"[2]", fie="cboxcent", val=xcboxcent//", "//ycboxcent, ver-, update+)

        i = 1
        param = i//"_xcent"
        imgets (image=stampname//"[2]", par=param)
        while (imgets.value != "0") {
	    pos = stridx ("+", imgets.value)
	    if (pos == 0) {
	        pos = strlen (imgets.value)
	        errstr = ""
	    } else {
		pos2 = strlen (imgets.value)
	   	errstr = substr (imgets.value, pos, pos2)
	    }
	    val = real (substr (imgets.value, 1, pos)) + offx
            hedit (ima=stampname//"[2]", fie=param, val=val//" "//errstr, ver-, update+)

	    param = i//"_ycent"
            imgets (image=stampname//"[2]", par=param)
	    pos = stridx ("+", imgets.value)
	    if (pos == 0) {
		pos = strlen (imgets.value)
	        errstr = ""
	    } else {
	    	pos2 = strlen (imgets.value)
		errstr = substr (imgets.value, pos, pos2)
	    }
	    val = real (substr (imgets.value, 1, pos)) + offy
            hedit (ima=stampname//"[2]", fie=param, val=val//" "//errstr, ver-, update+)

	    i=i+1
	    param = i//"_xcent"
            imgets (image=stampname//"[2]", par=param)
        }
    }
}

#  Delete all temporary images if the user wants this

if (nofat) {
    imdelete ("galfit-tempfits/temp-*.fits")
    imdelete ("galfit-tempfits/sig-*.fits")
    delete ("galfit-tempfits/offsets")
    !\rm obj??
    !\rm galfit.[0-9]?
    delete("perl.in")
    delete ("psf.temp")
    !\rmdir galfit-tempfits
}
    

print ("Done galfitting!")

end
