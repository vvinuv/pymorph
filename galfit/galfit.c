#define VERSION "2.0.3c -- Oct. 22, 2005"

/***********************************************************************\
* Note that all functions whose prototypes are declared before the      *
* start of the first function code will be found in an external file    *
* (usually because more than one function calls them).  Functions whose *
* prototypes are declared inside another function will be found in the  *
* same file as the calling function (usually because they are used only *
* by one calling function)						*
\***********************************************************************/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <curses.h>
#include "fitsio.h"
#include "nrutil.h"
#include "structs.h"
#include "sersic.h"
#include "debug.h"

void read_input (struct inpars *input, struct fitpars *fpar);
void assign_pars (struct inpars input, struct image *, struct image *, 
		  struct image *, struct image *, struct fitpars *,
		  struct convpars *);
int readfits (struct image *);
float sigcheck (struct image sigma);
void psfcheck (struct image *psf);
struct image *getkernel (char *kernel_file);
float psf_fwhm (struct image psf);
int badpixlist (struct image *badpix);
int read_cons (struct inpars input, struct cons *constr);
void errorcheck (struct inpars *input, struct image *data, 
	struct image *badpix, struct image *psf, struct image *sigma);
void getmu (float *, float *, float *, int);
void search_min (struct image *, struct image *, struct image *, 
		struct image *, struct fitpars *, struct cons *constr,
							struct convpars *);
void initcurses (char device[]);
void writefits (char *, struct image *, char *object, int);

struct inpars input = {"", "", "", "none", "", "", "none", 1, "none", "", 
                       {0., 0., 0.}, {0,0,0}, 0., {0., 0.}, 0, "both", 0, 0};
struct sampars sample = {0, 0};
float fwhmpsf, sigval;

int (*pfunc)(const char *, ...);
struct image psf, *kernel;
char device[10];

main (int argc, char *argv[])
{
    void blank_img (struct image *data, struct inpars input);

    char phead[FLEN_FILENAME], gname[FLEN_FILENAME];
    struct image data, sigma, badpix;
    struct fitpars fpar;
    struct convpars cpar;
    struct cons constr;
    int returnval;

    printf ("\nGALFIT Version %s\n\n", VERSION);

    if (argc > 1)
        strcpy (input.initparfile, argv[1]);

    /***************************************\
    *  Read in data, psf, and sigma images  *
    *  As well as initial fit parameters    *
    \***************************************/

    read_input (&input, &fpar); 
    assign_pars (input, &data, &sigma, &psf, &badpix, &fpar, &cpar);

    strcpy (phead, data.name);     /* Output images will have the same  */
				   /* header as the primary input image */
  
    /******************************************************\
    *  Set up output to either be in regular window mode,  *
    *  curses mode, or combination of both.                *
    \******************************************************/

    strcpy (device, input.device);
    if (strncmp (device, "regular", 7)!= 0) {
        initcurses(device);  /* Initialize the curses window */
        if (strncmp (device, "curses", 6)== 0)
            pfunc = printw;
        else
            pfunc = printf;
    } else
        pfunc = printf;

    pfunc ("\n\n");

    /*********************************************************************\
    *  Check to see that the input images exist, else deal with failures  *
    \*********************************************************************/

    if ((returnval = readfits (&data)) == 1) {
        pfunc ("-- WARNING: No input data image found.  Generating a model based\n ");
        pfunc ("  on input parameters, with exposure time of 1 second.\n\n");
	input.create = 1;
        sprintf (input.data, "none");
        blank_img (&data, input);
    } else if (returnval == 2) {
	pfunc ("-- WARNING: The exposure time header is missing.  Default to 1 second.\n\n");
    };

    if (readfits (&psf) == 1) {
	pfunc ("-- PSF image not found.  No convolution done! \n"); 
        sprintf (psf.name, "none");
    } else {
	/*  Make sure the PSF naxes are even  */
	psfcheck (&psf);
	if (input.sampfac > 1)
	    kernel = getkernel (input.kernel);

	if (input.sampfac > 1 && kernel ==  NULL)
	    pfunc ("-- No CCD charge diffusion kernel applied.\n");
	
        /* Keep track of PSF size to figure out minimal padding later  */
        cpar.psfsz[1] = psf.naxes[1];
        cpar.psfsz[2] = psf.naxes[2];
	fwhmpsf = psf_fwhm (psf);
    };

    if (readfits (&sigma) == 1) {
	pfunc ("-- No sigma image.  Creating one using: GAIN=%.2f, RDNOISE=%.2f, NCOMBINE=%.1f.\n", data.gain, data.rdnoise, data.ncombine);
    } else {
	if (sigval = sigcheck(sigma)) 
	    pfunc ("-- WARNING: Sigma image has 0. values.  Set to a low value of %.2e. \n", sigval);
    };

    /*************************************************************\
    *  Look for a badpix mask image.  Even if none found, create  *
    *  an array of 0s (by badpixlist) anyway for later use as     *
    *  a lookup table to prevent the same pixel from being        *  
    *  used twice.                                                *
    \*************************************************************/
   
    if (readfits (&badpix) == 1 || badpix.naxes[1] != data.naxes[1] ||
				   badpix.naxes[2] != data.naxes[2]) {
	if (badpix.err == 0) {
	    pfunc ("-- WARNING:  The pixel mask is not the same size as the data.\n\n");
	} else {
            badpix.naxes[1] = data.naxes[1];
            badpix.naxes[2] = data.naxes[2];
	    if (badpixlist(&badpix) == 1) {
	        pfunc ("-- No masking image used.\n");
	        sprintf (badpix.name, "none");
	    };
	};
    };

    /*******************************************************\
    *  Look for a file that contains parameter constraints  *
    \*******************************************************/

    if (read_cons (input, &constr) == 1)
	pfunc ("-- No constraint file found, or error in constraint file.\n");

    /******************************************************\
    *  Catch problems with the user input, specifically    *
    *  with image region sizes, that might cause the       *
    *  program to crash or cause unanticipated problems.   *
    \******************************************************/

    errorcheck (&input, &data, &badpix, &psf, &sigma);

    /***************************************************************\
    *  For Sersic profile only -- interpolates to get the relation  *
    *                           between the powerlaw index and mu   *
    \***************************************************************/

    getmu (TWOn, MU, Y2, NINTERP);

    /*****************************\
    *  Now the real fun begins!!  *
    \*****************************/

    search_min (&data, &sigma, &psf, &badpix, &fpar, &constr, &cpar); 

    /************************************************\
    *  DONE!  Now free up all the input matrices.    *
    \************************************************/

    free_matrix (data.z, 1, data.naxes[2], 1, data.naxes[1]);
    if (sigma.err == 0)
        free_matrix (sigma.z, 1, sigma.naxes[2], 1, sigma.naxes[1]);
    if ( strncmp (psf.name, "none", 4) != 0)
        free_matrix (psf.z, 1, psf.naxes[2], 1, psf.naxes[1]);
    if (kernel != NULL) 
	free_matrix (kernel->z, 1, kernel->naxes[2], 1, kernel->naxes[1]);

    free_matrix (badpix.z, 1, badpix.naxes[2], 1, badpix.naxes[1]);

    /* Beep when done */
    printf ("%c", 7);

    if (strncmp (device, "regular", 7)!= 0)
        endwin ();

}


void blank_img (struct image *data, struct inpars input)
{
    int xmin, xmax, ymin, ymax;

    sscanf (input.imgsect, "[%d:%d,%d:%d]", &xmin, &xmax, &ymin, &ymax);
    sprintf (data->imgsect, "[%d:%d,%d:%d]", &xmin, &xmax, &ymin, &ymax);
    data->naxes[1] = xmax;
    data->naxes[2] = ymax;
    data->z = matrix (1, data->naxes[2], 1, data->naxes[1]);
}
