#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <curses.h>
#include "fitsio.h"
#include "nrutil.h"
#include "structs.h"
#include "debug.h"

void printerror(int status);

/*************************************************************************/
/* Write a FITS image and determine the minimum and maximum pixel values */
/*    add = 0     ==>  new image                                         */
/*    add = 1     ==>  append to image                                   */
/*************************************************************************/

void writefits(char *phead, struct image *img, char *object, int append) 
{
    extern int (*pfunc)(const char *, ...);
    fitsfile *fptr, *hptr;  /* pointer to the FITS file, defined in fitsio.h */
    int status, i, j, ix, iy, dum;
    long fpixel, npix, exposure, naxes[2];
    float *vec;
    char comment[FLEN_COMMENT];
    char newobj[FLEN_CARD];

    /* initialize FITS image parameters */

    int bitpix = FLOAT_IMG; /* Single precision floating point values   */
    long naxis = 2;    /* 2-dimensional image                           */ 

    status = 0;        /* initialize status before calling fitsio routines */

    /* open file with the primary header parameters to use as a template */

    naxes[0] = img->naxes[1];
    naxes[1] = img->naxes[2];

    if (append) {         /* Append to old image? */

        if (fits_open_file (&fptr, img->name, READWRITE, &status)) {
	    pfunc ("\n Can't append to image %s!\n", img->name);
            printerror( status );
	};
        
        fits_open_file (&hptr, phead, READONLY, &status);
        fits_create_hdu(fptr, &status);
        fits_update_key(fptr, TSTRING, "XTENSION", "IMAGE",
       					  "IMAGE extension", &status);
	dum = -32;
        fits_update_key(fptr, TINT, "BITPIX", &dum,
       					  "Bits per pixel", &status);
	dum = 2;
        fits_update_key(fptr, TINT, "NAXIS", &dum,
       					  "Number of axes", &status);

        fits_update_key(fptr, TLONG, "NAXIS1", &naxes[0],
       					  "length of data axis 1", &status);
        fits_update_key(fptr, TLONG, "NAXIS2", &naxes[1],
       					  "length of data axis 2", &status);
        fits_update_key(fptr, TSTRING, "OBJECT", object,
       					  "object name", &status);
        if (status != 0)
            printerror( status );       

    } else {
	remove(img->name);  /* Delete old file if it already exists */

        if (fits_create_file(&fptr, img->name, &status))
	    printerror (status);
	if (fits_create_img (fptr, bitpix, 2, naxes, &status))
            printerror( status );

        /* Copy the primary header from the input file to the output */

        if (!fits_open_file (&hptr, phead, READONLY, &status)) {
            fits_copy_header(hptr, fptr, &status);
	    if ( fits_read_key(hptr, TSTRING, "OBJECT", newobj, comment, 
								&status) ) {
		status = 0;
	        if ( fits_read_key(hptr, TSTRING, "TARGNAME", newobj, comment, 
								&status) )
	            status = 0;
            }
	    fits_close_file (hptr, &status );
            sprintf (newobj, "%s%s", newobj, img->imgsect);
            if ( fits_update_key(fptr, TSTRING, "OBJECT", newobj, 
						"object name", &status) )
                printerror( status );
            fits_update_key(fptr, TLONG, "NAXIS1", &naxes[0],
       					  "length of data axis 1", &status);
            fits_update_key(fptr, TLONG, "NAXIS2", &naxes[1],
       					  "length of data axis 2", &status);

	} else {

	    status = 0;
            if ( fits_update_key(fptr, TSTRING, "OBJECT", object, 
						"object name", &status) )
                printerror( status );
	};
    };

    fpixel = 1;                                 /* first pixel to write */
    npix = img->naxes[1] * img->naxes[2];  /* number of pixels to write */

    /* write the array of floating point values to the FITS file.  Note */
    /* the conversion to zero offset with "img->flux+1"                 */

    i = 0;
    vec = vector (0, img->naxes[1] * img->naxes[2]-1);
    for (iy = 1; iy <= img->naxes[2]; iy++) {
        for (ix = 1; ix <= img->naxes[1]; ix++) {
	    vec[i++] = img->z[iy][ix];
	};
    };

    if ( fits_write_img(fptr, TFLOAT, fpixel, npix, vec, &status) )
        printerror (status);

    if ( fits_close_file(fptr, &status) )       /* close the file */
        printerror(status);

    free_vector (vec, 0, img->naxes[1] * img->naxes[2]-1);

    return;
}
