#include "fitsio.h"

#define OFFSET 1
#define NPARS 10
#define REAL 1
#define COMPLEX 2

/* Holds the primary input parameters specified by user */

struct inpars {
    char initparfile[FLEN_FILENAME];    /* name of the prepared input file */
    char data[FLEN_FILENAME];
    char output[FLEN_FILENAME];
    char sigma[FLEN_FILENAME];
    char psf[FLEN_FILENAME];
    char kernel[FLEN_FILENAME];		/* Diffusion kernel */
    char constraints[FLEN_FILENAME];    /* Name of the file with constraints */
    int  sampfac;
    char badpix[FLEN_FILENAME];
    char imgsect[70];			/* image section to fit */
    float cboxcent[3];			/* convolution box center */
    int convbox[3];			/* convolution box size */
    float magzpt;			/* magnitude zeropoint */
    float d[3];				/* plate scale */
    int coerr;				/* estimate error by parabol. method */
    char device[10];                    /* scrolling display type */
    int create;                         /* A boolean to see whether to
					   just create a lens model */
    int interact; 			/* A boolean to see whether to see
					   if the user wants command line 
					   interaction                      */ 

};


  /*******************************************************************\
  * NOTE: "**z" and "*flux" will point to the SAME memory address.    *
  *       Therefore if the flux changes in vector *z, so will **flux  *
  *       and vice versa.   ALL the arrays will use a one-offset      *
  *       convention                                                  *
  \*******************************************************************/

/* All the arrays will use a ONE-OFFSET Fortran convention */

struct image {
    char name[FLEN_FILENAME];
    float **z;          /* flux at pixel (x,y) */
    long naxes[3];      /* X and Y dimensions.  Do not use naxes[0] */
    char imgsect[70];   /* section cut-out from the original image */
    float exptime;      /* image exposure time of the original image */
    float magzpt;	/* magnitude zeropoint of the original image */
    float muzpt;	/* mag/arcsec^2 zeropoint */
    float dp[3];	/* arcsec per pixel.  Do not use dp[0] */
    float gain;         /* gain of the image in electrons/adu  */
    float rdnoise;      /* readnoise of the image in electrons */
    float ncombine;     /* number of images combined to make up the current */
    int err;            /* Error code, when the image is blank.           */
};


struct sampars {
    int nobjs;		/* Number of objects being fitted (including sky) */
    long nmask;         /* Total number of pixels masked out.  */
};

struct fitpars {
    char objtype[8];    /* Sersic, Nuker, Expdisk, deVauc,                 */
			/* Gaussian, Moffat, sky                           */
    float a[NPARS+OFFSET];   /* parameters                                 */
    float sig[NPARS+OFFSET]; /* uncertainties                              */
    int ia[NPARS+OFFSET];   /* 1 = fit parameter, 0 = hold parameter fixed */
			    /* -1 = not used                               */
    int outtype;        /* How to deal with the output of this object:     */
			/* 0 = create normal residual image                */
			/* 1 = Do not remove this component in the resid.  */
    struct fitpars *next;
};

/* Parameter constraints */

struct cons {
    int comp[3]; 	/* Component(s) involved in the constraints        */
    char op;		/* operation 					   */
    int par;		/* Parameter to constrain			   */
    float cval[3];	/* Constraint values (min, max)			   */
    struct cons *next;
};

/* Convolution parameters */

struct convpars {
    float cent[3];      /* x and y position of the convolution box    */
    int boxsz[3];       /* This is the valid length and width of the      */
			/* convolution region.   It's equal to            */
			/* inpars.convbox if it can fit within the image. */
    char imgsect[70];   /* Image section (relative to cut-out stamp      */
			/* coordinates) used in convolution.         */
    int sampfac;        /* PSF sampling factor */
    long psfsz[3];      /* x and y dimensions of the PSF image.      */
};

/* Arrays holding the derivative of the flux at each pixel  */
/* with respect to each of the parameters being fitted.     */

struct derivs {
    char imgsect[70];       /*  The section in the fitting region where  */
			    /*  one should copy the derivative images    */
			    /*  back into.                               */
    long naxes[3];          /*  The dimension of the derivative images.  */
    float **dpm[NPARS+1];   /*  2-D image                                */
    struct derivs *next;
};


/* Storage space for coefficients used to interpolate the PSF   */
/* using bi-cubic interpolation.                                */

struct bicubic_grid {
    float **c;
};
