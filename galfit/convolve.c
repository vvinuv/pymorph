/******************************************************************\
*  This subroutine convolves a model image with the PSF.  The      *
*  input PSF must already be Fourier transformed using the         *
*  subroutine psfprep.c.  Complex_naxes[] stores the size of the   *
*  FFT image, and the size is the number of real+complex pairs.    *
*  So the physical storage in the x-direction is twice as large    *
*  as the number in complex_naxes[2].                              *
\******************************************************************/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "structs.h"
#include "nrutil.h"
#include "mymath.h"

void cdft2d(int, int, int, double **, double *, int *, double *);
void copymat (float **, double **, long naxes[], int n);
double **psfprep (struct image psf, long model_naxes[], int add0pad,
    long complex_naxes[]);

void convolve (struct image *model, double **fftpsf, long complex_naxes[])
{

    void complmultipl(double *ar1, double *ar2, unsigned long ntot);
    void mkimg (float **array, long naxes[], char *outname, char *tname);

    long i, x, y, cx, ix, iy;
    struct fitpars *fptr;
    unsigned long ntot;
    int xmin, xmax, ymin, ymax, xupbound, xlobound, yupbound, ylobound;
    int *ip, n, nip;
    double *w, **fftwork;
    float **complex_img;

    /**************************************\
    *  Create FFT work space and clear it  *
    \**************************************/

    fftwork = dmatrix (0, complex_naxes[2]-1, 0, complex_naxes[1]*2-1);
    complex_img = matrix (1, complex_naxes[2], 1, complex_naxes[1]*2);

    for (y=0; y < complex_naxes[2]; ++y)
        for (x=0; x < complex_naxes[1] * 2; ++x)
	    fftwork[y][x] = 0.;

    /******************************************************\
    *  Now copy the region for convolution into a complex  *
    *  number array.                                       *
    \******************************************************/

    for (y=1; y <= model->naxes[2]; y++)
	for (x=1; x <= model->naxes[1]; x++) {
            cx = x * 2;
	    complex_img[y][cx-1] = model->z[y][x];
	    complex_img[y][cx] = 0.;
	};

    /**************************************\
    *  Next, blank out the padding region  *
    \**************************************/

    yupbound = complex_naxes[2];
    xlobound = IMIN(2*model->naxes[1] + 1, 2 * complex_naxes[1]);
    xupbound = 2*complex_naxes[1];

    for (y = 1; y <= yupbound; y++)
	for (x = xlobound; x <= xupbound; x++)
	    complex_img[y][x] = 0.;

    xupbound = IMIN ((model->naxes[1] + 1) * 2, 2 * complex_naxes[1]);
    ylobound = IMIN (model->naxes[2] + 1, complex_naxes[2]);
    yupbound = complex_naxes[2];

    for (y = ylobound; y <= yupbound; y++)
	for (x = 1; x <= xupbound; x++)
	    complex_img[y][x] = 0.;

    /************************************\
    *  Initialize convolution variables  *
    \************************************/

    ntot = complex_naxes[1] * complex_naxes[2];
    nip = IMAX (complex_naxes[2], complex_naxes[1]);
    ip = ivector (0, 2 + (int)sqrt(nip + 0.5)-1);
    n = IMAX (complex_naxes[2], complex_naxes[1] * 2) * 3 / 2;
    w = dvector (0, n-1);
    ip[0] = 0;

    /***************************************\
    *  Finally do the actual convolution    *
    \***************************************/

    copymat (complex_img, fftwork, complex_naxes, 1) ;
    cdft2d (complex_naxes[2], complex_naxes[1]*2, 1, fftwork, NULL, ip, w);
    complmultipl (*fftwork, *fftpsf, ntot); 
    cdft2d (complex_naxes[2], complex_naxes[1]*2, -1, fftwork, NULL, ip, w);
    copymat (complex_img, fftwork, complex_naxes, -1);

    /********************************************************************\
    *  Now copy the convolution results back into the input model array  *
    \********************************************************************/

    for (y = 1; y <= model->naxes[2]; y++){
        for (x = 1; x <= model->naxes[1]; x++){
            cx = x * 2;
            model->z[y][x] = complex_img[y][cx-1];
        };
    };

    /**************************************\
    *  Now free up convolution workspaces  *
    \**************************************/

    free_dmatrix (fftwork, 0, complex_naxes[2]-1, 0, complex_naxes[1]*2-1);
    free_matrix (complex_img, 1, complex_naxes[2], 1, complex_naxes[1]*2);
    free_ivector (ip, 0, 2 + (int)sqrt(nip + 0.5)-1);
    free_dvector (w, 0, n-1);

}


/****************************************************************************/

void complmultipl(double ar1[], double ar2[], unsigned long n) {

    double realpart, aimagpart;
    unsigned long i, k, j, ntot;

    ntot = 2 * n;
    for (i = 0; i< ntot; i+=2){
        k = i;
        j = i+1;
        realpart = ar1[k] * ar2[k] - ar1[j] * ar2[j];
        aimagpart = ar1[k] * ar2[j] + ar1[j] * ar2[k];
        ar1[k] = realpart;
        ar1[j] = aimagpart;
    };
}

/****************************************************************************/

void writefits(char *phead, struct image *img, char *object, int append);

void mkimg (float **array, long naxes[], char *outname, char *tname)
{

    struct image outimg;

    outimg.z = array;
    strcpy (outimg.name, outname);
    outimg.naxes[1] = naxes[1];
    outimg.naxes[2] = naxes[2];

    writefits (tname, &outimg, "", 0);
}


/****************************************************************************/

/******************************************************\
*  Now that we're done with the convolution, copy the  *
*  convolved region back to the mini-size image.       *
\******************************************************/

void copy_convregion (struct convpars *cpar, float **model, 
    struct image convolved_img, int os)
{
    extern struct image *kernel;

    int xlobound, xupbound, ylobound, yupbound, x, y, ix, iy, xmin, xmax, ymin,
	ymax, xbox, ybox, xx, yy;
    long kern_complex_naxes[3];
    struct image tmodel;
    double **fftkern;

    sscanf (convolved_img.imgsect,"[%d:%d,%d:%d]",&xmin,&xmax,&ymin,&ymax); 
    xbox = xmax - xmin + 1;
    ybox = ymax - ymin + 1;

    xlobound = (int) (cpar->psfsz[1]/2.);
    ylobound = (int) (cpar->psfsz[2]/2.);
    xupbound = xlobound + xbox * os;
    yupbound = ylobound + ybox * os;

    tmodel.z = matrix (1, ybox, 1, xbox);
    tmodel.naxes[1] = xbox;
    tmodel.naxes[2] = ybox;

    for (y = ylobound; y < yupbound; y+=os){
        for (x = xlobound; x < xupbound; x+=os){
	    ix = (x - xlobound)/os + 1;
            iy = (y - ylobound)/os + 1;

	    tmodel.z[iy][ix] = 0.;
	    for (yy = y; yy < y+os; yy++) {
		for (xx = x; xx < x+os; xx++) 
	            tmodel.z[iy][ix] += convolved_img.z[yy+1][xx+1];
	    };
	};
    };

    /*  Apply charge diffusion kernel  */
    if (kernel != NULL) {
        fftkern = psfprep (*kernel, tmodel.naxes, 1, kern_complex_naxes);
        convolve (&tmodel, fftkern, kern_complex_naxes);

        free_dmatrix (fftkern, 0, kern_complex_naxes[2]-1, 0,
            kern_complex_naxes[1]*2-1);
    };

    for (iy = 1; iy <= tmodel.naxes[2]; iy++) {
        for (ix = 1; ix <= tmodel.naxes[1]; ix++) {
	    x = xmin + ix - 1;
	    y = ymin + iy - 1;
	    model[y][x] = tmodel.z[iy][ix];
	};
    };

    free_matrix (tmodel.z, 1, ybox, 1, xbox);
}

