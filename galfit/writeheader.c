#include <stdlib.h>
#include <string.h>
#include <curses.h>
#include "fitsio.h"
#include "structs.h"

void writeheader (struct image *img, struct fitpars *fpar, double chisq,
					int ndof, int hdu, int offx, int offy)
{
    void objects_head (struct fitpars *fpar, fitsfile *fptr, 
							int offx, int offy);
    extern struct inpars input;
    char param[FLEN_CARD];
    fitsfile *fptr;
    int status = 0;
    float chi2nu;

    if (fits_open_file (&fptr, img->name, READWRITE, &status)) {
        pfunc ("\n Can't append to image header %s!\n", img->name);
        printerror( status );
    };

    if (fits_movabs_hdu (fptr, hdu, IMAGE_HDU, &status)) {
        pfunc ("\n Can't append to image header %s!\n", img->name);
        printerror( status );
    };

    fits_write_comment (fptr, " ", &status);
    fits_write_comment (fptr, "========== GALFIT Input Parameters ==========", &status);
    fits_write_comment (fptr, "", &status);

    fits_update_key(fptr, TSTRING, "INITFILE", input.initparfile,
                                          "GALFIT input file", &status);

    fits_update_key(fptr, TSTRING, "DATAIN", input.data,
                                          "Input data image", &status);

    fits_update_key(fptr, TSTRING, "SIGMA", input.sigma,
                                          "Input sigma image", &status);

    fits_update_key(fptr, TSTRING, "PSF", input.psf,
                                          "Convolution PSF", &status);

    fits_update_key(fptr, TSTRING, "CONSTRNT", input.constraints,
                                          "Parameter constraint file", &status);

    fits_update_key(fptr, TSTRING, "MASK", input.badpix,
                                          "Input mask image", &status);

    fits_update_key(fptr, TSTRING, "FITSECT", input.imgsect,
                                          "Image section fitted", &status);

    sprintf (param, "%-.2f, %-.2f", input.cboxcent[1], input.cboxcent[2]);

    fits_update_key(fptr, TSTRING, "CBOXCENT", param,
                                          "Convolution box center", &status);

    sprintf (param, "%i, %i", input.convbox[1], input.convbox[2]);

    fits_update_key(fptr, TSTRING, "CONVBOX", param,
                                          "Convolution box size", &status);

    fits_update_key(fptr, TFLOAT, "MAGZPT", &input.magzpt,
                                          "Magnitude zeropoint", &status);

    fits_write_comment (fptr, "========== GALFIT Final Parameters ==========", &status);

    objects_head (fpar, fptr, offx, offy);

    fits_update_key(fptr, TDOUBLE, "Chisq", &chisq,
                                          "Chi^2 of fit", &status);
    fits_update_key(fptr, TINT, "NDOF", &ndof,
                                          "Degrees of Freedom", &status);

    chi2nu = chisq/ndof;
    fits_update_key(fptr, TFLOAT, "Chi2nu", &chi2nu,
                                          "Reduced Chi^2", &status);

    fits_write_comment (fptr, "=============================================", &status);
    fits_write_comment (fptr, " ", &status);

    fits_close_file(fptr, &status);                   /* close the file */

    return;
}



void objects_head (struct fitpars *fpar, fitsfile *fptr, int offx, int offy)
{
    char par[FLEN_CARD], skyc[FLEN_CARD], comp[FLEN_KEYWORD], 
	 key[FLEN_KEYWORD];
    float parval;
    extern float xskycent, yskycent;
    int status = 0, i, j;

    j = 1;
    while (fpar != NULL) {
        sprintf (comp, "COMP_%i", j);
        fits_update_key(fptr, TSTRING, comp, fpar->objtype,
                                              "Object type", &status);
        for (i=1; i<= 10; i++) {
	    if (i==1 && strncmp(fpar->objtype, "sky", 3) != 0) 
		parval = fpar->a[i] + offx;
	    else if (i==2 && strncmp(fpar->objtype, "sky", 3) != 0)
		parval = fpar->a[2] + offy;
	    else
		parval = fpar->a[i];
	    
	    sprintf (par, "%-.4f +/- %-.4f", parval, fpar->sig[i]);

            if (strncmp(fpar->objtype, "sky", 3) != 0) {

	        if (i==1) {
	            sprintf (key, "%i_XCENT", j);
                    fits_update_key(fptr, TSTRING, key, par, "X center [pixel]", 
								&status);
	        };

	        if (i==2) {
	            sprintf (key, "%i_YCENT", j);
                    fits_update_key(fptr, TSTRING, key, par, "Y center [pixel]", 
								&status);
	        };

                /*******\
                * Nuker *
                \*******/

	        if (strncmp(fpar->objtype, "nuker", 5) == 0) {
		    if (i==3) {
	                sprintf (key, "%i_MU", j);
                        fits_update_key(fptr, TSTRING, key, par,
                                          "Surface brightness at Rb", &status);
		    };

		    if (i==4) {
	                sprintf (key, "%i_Rb", j);
                        fits_update_key(fptr, TSTRING, key, par,
                                          "Break radius Rb [pixels]", &status);
		    };

		    if (i==5) {
	                sprintf (key, "%i_alpha", j);
                        fits_update_key(fptr, TSTRING, key, par,
                                          "alpha", &status);
		    };

		    if (i==6) {
	                sprintf (key, "%i_beta", j);
                        fits_update_key(fptr, TSTRING, key, par,
                                          "beta", &status);
		    };

		    if (i==7) {
	                sprintf (key, "%i_gamma", j);
                        fits_update_key(fptr, TSTRING, key, par,
                                          "gamma", &status);
		    };
	        };

                /**************************\
                * Sersic and deVaucouleurs *
                \**************************/

	        if (strncmp(fpar->objtype, "sersic", 6) == 0 ||
				strncmp(fpar->objtype, "devauc", 6) == 0) {
		    if (i==3) {
	                sprintf (key, "%i_MAG", j);
                        fits_update_key(fptr, TSTRING, key, par,
                                          "Total magnitude", &status);
		    };

		    if (i==4) {
	                sprintf (key, "%i_Re", j);
                        fits_update_key(fptr, TSTRING, key, par,
                                          "Effective radius Re [pixels]", &status);
		    };

		    if (i==5 && strncmp(fpar->objtype, "sersic", 6) == 0) {
	                sprintf (key, "%i_n", j);
                        fits_update_key(fptr, TSTRING, key, par,
                                          "Sersic index", &status);
		    };
	        };

                /******************\
                * Exponential Disk *
                \******************/

	        if (strncmp(fpar->objtype, "expdisk", 7) == 0) {
		    if (i==3) {
	                sprintf (key, "%i_MAG", j);
                        fits_update_key(fptr, TSTRING, key, par,
                                          "Total magnitude", &status);
		    };

		    if (i==4) {
	                sprintf (key, "%i_Rs", j);
                        fits_update_key(fptr, TSTRING, key, par,
                                          "Scalelength [pixels]", &status);
		    };
		};

                /*********************\
                * Gaussian and Moffat *
                \*********************/

	        if (strncmp(fpar->objtype, "gauss", 5) == 0 || 
		    strncmp(fpar->objtype, "moffat", 6) == 0 ||
		    strncmp(fpar->objtype, "psf", 3) == 0 ) {

		    if (i==3) {
	                sprintf (key, "%i_MAG", j);
                        fits_update_key(fptr, TSTRING, key, par,
                                          "Total magnitude", &status);
		    };

		    if (i==4 && strncmp(fpar->objtype, "psf", 3) != 0) {
	                sprintf (key, "%i_FWHM", j);
                        fits_update_key(fptr, TSTRING, key, par,
                                          "FWHM [pixels]", &status);
		    };

		    if (i==5 && strncmp(fpar->objtype, "moffat", 6) == 0) {
	                sprintf (key, "%i_C", j);
                        fits_update_key(fptr, TSTRING, key, par,
                                          "Powerlaw", &status);
		    };
		};


                /*******************\
                * Common parameters *
                \*******************/

	        if (i==8 && strncmp(fpar->objtype, "psf", 3) != 0) {
	            sprintf (key, "%i_AR", j);
                    fits_update_key(fptr, TSTRING, key, par,
                                      "Axis ratio (b/a)", &status);
	        };

	        if (i==9 && strncmp(fpar->objtype, "psf", 3) != 0) {
	            sprintf (key, "%i_PA", j);
                    fits_update_key(fptr, TSTRING, key, par,
                      "Position Angle (PA) [Degrees: Up=0, Left=90]", &status);
	        };

	        if (i==10 && strncmp(fpar->objtype, "psf", 3) != 0) {
	            sprintf (key, "%i_C", j);
                    fits_update_key(fptr, TSTRING, key, par,
                                  "Diskiness (<0) or Boxiness (>0)", &status);
	        };
            };

            /*****\
            * Sky *
            \*****/

            if (strncmp(fpar->objtype, "sky", 3) == 0) {
		if (i==1) {

	            sprintf (skyc, "%-.4f", xskycent + offx);
	            sprintf (key, "%i_XCENT", j);
                    fits_update_key(fptr, TSTRING, key, skyc, "X center [pixel]", 
								&status);

	            sprintf (skyc, "%-.4f", yskycent + offy);
	            sprintf (key, "%i_YCENT", j);
                    fits_update_key(fptr, TSTRING, key, skyc, "Y center [pixel]", 
								&status);

	            sprintf (key, "%i_SKY", j);
                    fits_update_key(fptr, TSTRING, key, par, "Sky background [ADUs]", 
								&status);
		};

		if (i==2) {
	            sprintf (key, "%i_DSKYDX", j);
                    fits_update_key(fptr, TSTRING, key, par, "x sky gradient [ADUs]", 
								&status);
		};

		if (i==3) {
	            sprintf (key, "%i_DSKYDY", j);
                    fits_update_key(fptr, TSTRING, key, par, "y sky gradient [ADUs]", 
								&status);
		};
	    };
        };

	fits_write_comment (fptr, "------------------------", &status);
	fpar = fpar->next;
	j++;
    };
}
