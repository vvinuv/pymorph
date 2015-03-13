#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fitsio.h"
#include "nrutil.h"
#include "structs.h"
#include "debug.h"

void printerror (int status);

/************************************************************************/
/* Read a FITS image and determine the minimum and maximum pixel values */
/************************************************************************/

int readfits (struct image *img)

{
    fitsfile *fptr;
    int status,  nfound, anynull=0, errstat=0, ix, iy;
    float *vec;
    long fpixel, npix, i, j, k, naxes[3], *nptr;
    char comment[FLEN_COMMENT];

    float datamin, datamax, nullval;

    status = 0;

    if (fits_open_file(&fptr, img->name, READONLY, &status)) {
	img->err = 1;
	return (1);
    } else
	img->err = 0;

    /* Read the image keywords to determine size, gain, rdnoise, etc.. */

    /* Note I want to put the output values in naxes[1] and naxes[2] */
    /* instead of the default naxes[0] and naxes[1] locations        */

    if (fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes+1, &nfound, &status))
        printerror( status );
    img->naxes[1]=naxes[1];
    img->naxes[2]=naxes[2];

    if (fits_read_key(fptr, TFLOAT, "GAIN", &(img->gain), comment, &status)){
 	status=0;
        if (fits_read_key(fptr, TFLOAT, "CCDGAIN", &(img->gain), comment, &status)){
	    status=0;
            if (fits_read_key(fptr, TFLOAT, "ATODGAIN", &(img->gain), comment, &status)) {
		img->gain = 7.;
		status = 0;
	    };
	};
    };
    if (fits_read_key(fptr, TFLOAT, "RDNOISE", &(img->rdnoise), comment, &status)) {
	img->rdnoise = 5.24;
	status = 0;
    };
    if (fits_read_key(fptr, TFLOAT, "NCOMBINE", &(img->ncombine), comment, &status)) {
	img->ncombine = 1.;
	status = 0;
    };
    if (fits_read_key(fptr, TFLOAT, "EXPTIME", &(img->exptime), comment, &status)) {
	img->exptime = 1.;
	status = 0;
	errstat = 2;      /* Return the fact the exposure time is missing */
    } else {
	img->magzpt += 2.5 * log10 (img->exptime);    /* adjust the mag */
	img->muzpt +=  2.5 * log10 (img->exptime);    /* zeropoints     */
    };

    /* number of pixels in the image */
    npix = img->naxes[1] * img->naxes[2];  

    fpixel   = 1;
    nullval  = 0;       /* don't check for null values in the image */
    datamin  = 1.0E30;
    datamax  = -1.0E30;

    /* Create one-offset vector and matrix.  Notice the needed "&"  */
    /* in passing the pointer to the address which then points to   */ 
    /* a newly created address space.                               */

    img->z = matrix (1, img->naxes[2], 1, img->naxes[1]);
    vec = vector (0, img->naxes[2] * img->naxes[1]-1);

    /* Note that because fits_read_img needs zero offset vectors,  */
    /* one is added to img->flux pointer.                          */

    if ( fits_read_img(fptr, TFLOAT, fpixel, npix, &nullval,
    						vec, &anynull, &status) ) {
        printerror( status );
    }

    i = 0;
    for (iy = 1; iy <= img->naxes[2]; iy++) {
	for (ix = 1; ix <= img->naxes[1]; ix++) {
	    img->z[iy][ix] = vec[i++];
	};
    };

    free_vector (vec, 0, img->naxes[2] * img->naxes[1]-1);

    if ( fits_close_file(fptr, &status) )
        printerror( status );

    return (errstat);   /* Normal return */
}


