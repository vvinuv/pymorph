/***********************************************************************\
*  Figure out how big of a convolution box we can have for each object  *
\***********************************************************************/

#include <math.h>
#include <stdio.h>
#include "structs.h"
#include "nrutil.h"
#include "mymath.h"

void convregion (struct derivs *df, float a[], struct image *d, 
    struct convpars *cpar)
{
    int nx, ny, xmin, xmax, ymin, ymax, xlo, ylo, xhi, yhi;
    extern struct inpars input;

    /****************************************************************\
    *                                                                *
    *  Make sure the user-specified convolution box can fit within   *
    *  the stamp sized data image.  If not, trim the convolution     *
    *  box to fit, i.e. keep track of the convolution box size that  *
    *  is within the stamp sized image.  The center of the           *
    *  convolution box is always centered on object centers.         *
    *                                                                *
    \****************************************************************/

    sscanf (d->imgsect, "[%d:%d,%d:%d]", &xlo, &xhi, &ylo, &yhi);

    if ( (xmax = NINT (input.convbox[1]/2.+a[1])-1) > d->naxes[1] )
	xmax = d->naxes[1];

    if ( (ymax = NINT (input.convbox[2]/2.+a[2])-1) > d->naxes[2] )
	ymax = d->naxes[2];

    if ( (xmin = NINT (a[1] - input.convbox[1]/2.)) < 1 )
	xmin = 1;

    if ( (ymin = NINT (a[2] - input.convbox[2]/2.)) < 1 )
	ymin = 1;

    /*  The entire convolution box is outside the fitting region...  */

    if ((xmin > d->naxes[1] && xmax == d->naxes[1]) ||
	(ymin > d->naxes[2] && ymax == d->naxes[2]) ||
	(xmax < 1 && xmin == 1) ||
	(ymax < 1 && ymin == 1)) {
        sprintf (df->imgsect, "[%d:%d,%d:%d]", 0, 0, 0, 0); 
	df->naxes[1] = 0;
	df->naxes[2] = 0;

    } else {

        sprintf (df->imgsect, "[%d:%d,%d:%d]", xmin, xmax, ymin, ymax); 

        /********************************************************************\
        *  The full image size includes the fitting region plus PSF padding  *
        *  because we have to generate model values in this region rather    *
        *  than just set them to 0.                                          *
        \********************************************************************/

        df->naxes[1] = (xmax - xmin + 1) * input.sampfac + 
		       2 * NINT(cpar->psfsz[1]/2.);
        df->naxes[2] = (ymax - ymin + 1) * input.sampfac + 
		       2 * NINT(cpar->psfsz[2]/2.);
    };
}
