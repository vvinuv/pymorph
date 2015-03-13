/******************************************************************\
*  This subroutine takes as input a PSF, the size of the image     *
*  being convolved, and a flag for whether or not to add 0         *
*  padding.  Then this subroutine outputs a Fourier transformed    *
*  version of the PSF which is used direcly by convolve.c, and	   *
*  the size of the FFT working image needed.                       *
* 								   *
*  If add0pad = 1, the model image is padded additionally by a     *
*  buffer zone of 0's the size of the PSF in both x and y before   *
*  the image is extended up to the next 2^n, where n is an         *
*  integer.  This is done so that convolution doesn't corrupt      *
*  parts of the model image.  Sometimes, however, the model        *
*  image has already been padded with *real* model values          *
*  instead of 0, then setting add0pad = 0, the model image will    *
*  not have a buffer zone of 0 padding, except for what's needed   *
*  to extend the image up to the next 2^n -- a requirement of      *
*  the FFT algorithm.                                              *
\******************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "nrutil.h"
#include "structs.h"
#include "mymath.h"

#define STRLEN 100

void cdft2d(int, int, int, double **, double *, int *, double *);
void writefits(char *phead, struct image *img, char *object, int append) ;

double **psfprep (struct image psf, long model_naxes[], int add0pad,
    long complex_naxes[])
{
    extern int (*pfunc)(const char *, ...);
    unsigned long x, y;
    float **temp, sum=0., psfnorm=0.;
    double *w;
    int *ip, n, nip;
    unsigned long osize[3], ohalf[3];
    double **fftpsf;

    osize[1] = psf.naxes[1]; 
    osize[2] = psf.naxes[2];

    /***************************************************************\
    *  Figure out how big of a complex image working space we need  *
    *  for convolution with 0 padding and all.                      *
    \***************************************************************/

    n = 0;
    if (add0pad)
        while (pow (2, ++n) < model_naxes[1] + 2*NINT(psf.naxes[1]/2.));
    else if (!add0pad)
        while (pow (2, ++n) < model_naxes[1]);
    complex_naxes[1] = pow(2, n);

    n = 0;
    if (add0pad)
        while (pow (2, ++n) < model_naxes[2] + 2*NINT(psf.naxes[2]/2.));
    else if (!add0pad)
        while (pow (2, ++n) < model_naxes[2]);
    complex_naxes[2] = pow(2, n);

    fftpsf = dmatrix (0, complex_naxes[2]-1, 0, complex_naxes[1]*2-1);

    for (y=0; y < complex_naxes[2]; ++y)
        for (x=0; x < complex_naxes[1] * 2; ++x)
            fftpsf[y][x] = 0.;

    /************************************************************\
    * Create space for a FFT PSF image with space multiple of 2. * 
    * "fftpsf.z" has 2*osize[1] storage spaces along the x and  *
    * osize[2] spaces along the y to hold both the real and      *
    * imaginary parts of the complex array.                      *
    \************************************************************/

    temp = matrix (1, complex_naxes[2], 1, complex_naxes[1]); 

    for (y=1; y <= psf.naxes[2]; y++)
        for (x=1; x <=psf.naxes[1]; x++)
	    sum += psf.z[y][x];

    for (y=1; y <= complex_naxes[2]; ++y)
	for (x=1; x <= complex_naxes[1]; ++x)
	    temp[y][x] = 0.;

    /*********************************************************\
    *                                                         *   
    *  Now distribute the PSF, centered on "half" to corners  *
    *  of the image.                                          *
    *                                                         *   
    \*********************************************************/

    ohalf[1] = NINT(osize[1]/2.);
    ohalf[2] = NINT(osize[2]/2.);

    /* corner 1 -- top right of orig to bottom left of new */
    for (y=ohalf[2]+1; y <= osize[2]; y++) {
	for (x=ohalf[1]+1; x <= osize[1]; x++) {
            temp[y-ohalf[2]][x-ohalf[1]] = psf.z[y][x]/sum;
	};
    };

    /* corner 2 -- bottom left of orig to top right of new */
    for (y=ohalf[2]; y >= 1; y-- ) {
	for (x=ohalf[1]; x >= 1; x--) {
            temp[complex_naxes[2]-ohalf[2]+y][complex_naxes[1]-ohalf[1]+x] = 
		psf.z[y][x]/sum;
	};
    };

    /* corner 3 -- bottom right of orig to top left of new */
    for (y=ohalf[2]; y >= 1; y--) {
        for (x=ohalf[1]+1; x <= osize[1]; ++x){
	    temp[complex_naxes[2]-ohalf[2]+y][x - ohalf[1]] = psf.z[y][x]/sum;
	};
    };

    /* corner 4 -- top left of orig to bottom right of new */
    for (y=ohalf[2]+1; y <= osize[2]; ++y){
	for (x=ohalf[1]; x>= 1; x--){
	    temp[y-ohalf[2]][complex_naxes[1]-ohalf[1]+x] = psf.z[y][x]/sum;
	};
    };

    for (y=1; y <= complex_naxes[2]; ++y) {
	for (x=1; x <= complex_naxes[1]; ++x)
	    fftpsf[y-1][(x-1)*2] = temp[y][x];
    };

    free_matrix(temp, 1, complex_naxes[2], 1, complex_naxes[1]);

#if DEBUG
  sprintf (fftpsf.name, "prefft-psf.fits");  /* Output image */
  writefits ("test.fits", fftpsf, "", 0);
#endif

    /* Do FFT now */

    nip = IMAX (complex_naxes[2], complex_naxes[1]);
    ip = ivector (0, 2 + (int)sqrt(nip + 0.5)-1);
    n = IMAX (complex_naxes[2], complex_naxes[1] * 2) * 3 / 2;
    w = dvector (0, n-1);
    ip[0] = 0;
    
    cdft2d (complex_naxes[2], complex_naxes[1]*2, 1, fftpsf, NULL, ip, w);

    /* Now normalize the PSF */
    for (y=0; y < complex_naxes[2]; ++y)
	for (x=0; x < complex_naxes[1]*2; ++x)
	    fftpsf[y][x] = fftpsf[y][x] / (complex_naxes[2] * 
			   complex_naxes[1]);

    free_ivector (ip, 0, 2 + (int)sqrt(nip + 0.5)-1);
    free_dvector (w, 0, n-1);

#if DEBUG 
  for (y=0; y < complex_naxes[2]; ++y) {
     for (x=0; x < complex_naxes[1]; ++x) {
	    psfnorm += fftpsf[y][x];
	};
  };

  printf ("PSF normalization: %f\n", psfnorm);

  #if DEBUG 
    sprintf (psf.name, "orig-psf.fits");    /* Output image */
    writefits ("test.fits", psf, "", 0);
    sprintf (fftpsf->name, "postfft-psf.fits");  /* Output image */
    writefits ("test.fits", fftpsf, "", 0);
  #endif
#endif

    return (fftpsf);
}

/*===========================================================================*/

/***************************************************************************\
*  Embed a PSF image that has odd number of pixels in length/width into an  *
*  image that is even in both directions..                                  *
\***************************************************************************/

void psfcheck (struct image *psf)
{
    int i, j, xoff=0, yoff=0;
    long new_naxes[3];
    float **newimg;

    if ((psf->naxes[1] % 2) != 0) {
	new_naxes[1] = psf->naxes[1] + 1;
	xoff = 1;
    } else
	new_naxes[1] = psf->naxes[1];

    if ((psf->naxes[2] % 2) != 0) {
	new_naxes[2] = psf->naxes[2] + 1;
	yoff = 1;
    } else
	new_naxes[2] = psf->naxes[2];

    if (xoff == 1 || yoff == 1) {
	newimg = matrix (1, new_naxes[2], 1, new_naxes[1]);

        for (j=1; j <= new_naxes[2]; j++) {
	    for (i=1; i <= new_naxes[1]; i++) {

		if ( (i==1 && xoff == 1) ||
		    (j==1 && yoff == 1 ) )

		    newimg[j][i] = 0.;
		else 
		    newimg[j][i] = psf->z[j-yoff][i-xoff];
	    };
        };

	free_matrix (psf->z, 1, psf->naxes[2], 1, psf->naxes[1]);
	psf->naxes[1] = new_naxes[1];
	psf->naxes[2] = new_naxes[2];
	psf->z = newimg;
    };
}


/*===========================================================================*/

/***************************************************************************\
*  Get the charge diffusion kernel of the PSF.  This kernel is convolved    *
*  with image model *after* GALFIT convolves the oversampled PSF with the   *
*  model and rebins the model down to science resolution.                   *
\***************************************************************************/

struct image *getkernel (char *input_kernel)
{
    FILE *kfile;
    char string[STRLEN], s2[] = " ", *p;
    struct image *kernel;
    int nx=0, ny=0, nc, i, j, ix, iy;

    kfile = fopen (input_kernel, "r");
    if (kfile != (FILE *) 0) {

	/****************************************************************\
	*  Figure out how many rows and columns there are in the kernel  *
	\****************************************************************/

        while ( ! feof (kfile) ) {
            if ( fscanf (kfile, " %[^\n]", string) > 0 && 
		strncmp (string, " ", 1) != 0 && 
		strncmp (string, "\n", 1) != 0) {

	        ++ny;
	        nx = 1;
		strtok (string, s2);
		while (strtok (NULL, s2) != NULL)
		    ++nx;

		/*******************************************\
		*  Error: number of columns isn't the same  *
		\*******************************************/

		if (ny > 1 && nx != nc) {
		    close (kfile);
		    return (NULL);
		} else if (ny==1)
		    nc = nx;
	    };
	};
	close (kfile);

	if (ny<=1 && nx<=1)   /*  No kernel */
	    return (NULL);

	/************************************************\
	*  Now reopen the file and read in the kernel.   *
	\************************************************/

        kfile = fopen (input_kernel, "r");
	kernel = (struct image *) malloc ((size_t)(sizeof (struct image)));
	kernel->naxes[1] = NINT(nx/2.)*2;
	kernel->naxes[2] = NINT(ny/2.)*2;
	kernel->z = matrix (1, kernel->naxes[2], 1, kernel->naxes[1]);

	for (j = 1; j <= kernel->naxes[2]; j++) {
	    for (i = 1; i <= kernel->naxes[1]; i++)
		kernel->z[j][i] = 0.;
	};

	if (ny%2 != 0)
	    iy = 1;
	else 
	    iy = 0;

        while ( ! feof (kfile) ) {

	    if (nx%2 != 0)
	        ix = 2;
	    else
	        ix = 1;
            if ( fscanf (kfile, " %[^\n]", string) > 0 ) {
	        ++iy;
		kernel->z[iy][ix] = atof (strtok(string, s2));
		while ((p=strtok (NULL, s2)) != NULL) {
		    ++ix;
		    kernel->z[iy][ix] = atof (p);
		};
	    };
	};

	fclose (kfile);
	return (kernel);

    } else
	return (NULL);

}
