#include <stdlib.h>
#include <string.h>
#include "structs.h"

int (*pfunc)(const char *, ...);

void errorcheck (struct inpars *input, struct image *data, 
		struct image *badpix, struct image *psf, struct image *sigma)
{
    int xmin, xmax, ymin, ymax, errn=0;
    float yhalfbox, xhalfbox;
    long *dnaxes, *pnaxes;

    /*******************************************************************\
    *  First check to see that fitting region is within image boundary  *
    \*******************************************************************/

    dnaxes = data->naxes;
    pnaxes = psf->naxes;

    sscanf (input->imgsect, "[%d:%d,%d:%d]", &xmin, &xmax, &ymin, &ymax);
    if (xmin < 1) {
        xmin = 1;
        errn = 1;
    };
    if (ymin < 1) { 
        ymin = 1;
        errn = 1;
    };
    if (xmax > dnaxes[1]) {
        xmax = dnaxes[1];
        errn = 1;
    };
    if (ymax > dnaxes[2]) {
        ymax = dnaxes[2];
        errn = 1;
    };

    sprintf (input->imgsect, "[%d:%d,%d:%d]", xmin, xmax, ymin, ymax);
    strcpy (data->imgsect, input->imgsect);
    strcpy (badpix->imgsect, input->imgsect);
    strcpy (sigma->imgsect, input->imgsect);

    if (xmax == xmin || ymax == ymin) {
        pfunc ("\n");
        pfunc ("The fitting region has either zero width or is completely outside of image.\n");  
        pfunc ("Quitting now.\n\n");
	exit (2);
    } else if (errn == 1) {
        pfunc ("\n");
	pfunc ("-- WARNING: Fitting box exceeds image boundary; vignetting to fit.\n");
    };

    /******************************************************************\
    *  Check to see that convolution box is larger than the PSF used   *
    \******************************************************************/

    if (input->convbox[1] < pnaxes[1]/input->sampfac || 
	input->convbox[2] < pnaxes[2]/input->sampfac) {

        pfunc ("\n");
        pfunc ("-- WARNING: Convolution PSF kernel exceeds the convolution box size.\n");
    };   

}
