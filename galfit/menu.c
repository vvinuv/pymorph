#include <stdio.h>
#include <string.h>
#include "structs.h"

int totnobj;

void menu (struct inpars *input, struct fitpars *fpar, FILE *ftype, 
			float xoffset, float yoffset, float chisq, int ndof)
{
    int box[5]={0, 0, 0, 0, 0}, i = 0;
    float x, y;
    char object[8];

    if (chisq > 0.)
        fprintf (ftype, "\n#  Chi^2/nu = %.3f,  Chi^2 = %.3f,  Ndof = %d",
						chisq/ndof, chisq, ndof);
    fprintf (ftype, "\n\n");
    fprintf (ftype, "================================================================================\n");
    fprintf (ftype, "# IMAGE and GALFIT CONTROL PARAMETERS\n");
    fprintf (ftype, "A) %-14s      # Input data image (FITS file)\n", input->data);
    fprintf (ftype, "B) %-14s      # Output data image block\n", input->output);
    fprintf (ftype, "C) %-14s      # Sigma image name ", input->sigma);
    fprintf (ftype, "(made from data if blank or \"none\") \n");
    fprintf (ftype, "D) %-10s %-8s # Input PSF image and ", input->psf, input->kernel);
    fprintf (ftype, "(optional) diffusion kernel\n");
    fprintf (ftype, "E) %-14d      # PSF oversampling factor relative to data \n", input->sampfac);
    fprintf (ftype, "F) %-14s      # Bad pixel mask (FITS image or ASCII coord list)\n", input->badpix);
    fprintf (ftype, "G) %-14s      # File with parameter constraints (ASCII file) \n", input->constraints);
    sscanf (input->imgsect, "[%d:%d,%d:%d]", &box[1], &box[2], &box[3], &box[4]);
    fprintf (ftype, "H) %-4d %-4d %-4d %-4d # Image region to fit (xmin xmax ymin ymax)\n", 
					box[1], box[2], box[3], box[4]);
    fprintf (ftype, "I) %-6d %-6d       # Size of the convolution box (x y)\n", 
					input->convbox[1], input->convbox[2]); 
    fprintf (ftype, "J) %-14.3f      # Magnitude photometric zeropoint \n", input->magzpt);
    fprintf (ftype, "K) %-6.3f %-6.3f       # Plate scale ", input->d[1], input->d[2]);
    fprintf (ftype, "(dx dy) \n");
    fprintf (ftype, "O) %-14s      # Display type (regular, curses, both)\n", input->device);
    fprintf (ftype, "P) %-14d      # Create ouput only? (1=yes; 0=optimize) \n", input->create);
    fprintf (ftype, "S) %-14d      # Modify/create objects interactively?\n", input->interact);
    fprintf (ftype, "\n");
    fprintf (ftype, "# INITIAL FITTING PARAMETERS\n");
    fprintf (ftype, "#\n");
    fprintf (ftype, "#   For object type, the allowed functions are: \n");
    fprintf (ftype, "#       nuker, sersic, expdisk, devauc, king, psf, gaussian, and moffat.\n");
    fprintf (ftype, "\n");
    fprintf (ftype, "# Objtype:      Fit?         Parameters \n");
    fprintf (ftype, "\n");
    while ( fpar != NULL && strncmp (fpar->objtype, "none", 4) != 0 ) {
        fprintf (ftype, "# Object number: %d\n", ++i);
	fprintf (ftype, " 0) %10s         #    Object type\n", fpar->objtype);
        if (strncmp (fpar->objtype, "sky", 3) != 0) {
            x = fpar->a[1] + xoffset;
            y = fpar->a[2] + yoffset;
	};

        if (strncmp (fpar->objtype, "nuker", 5) == 0) {

            fprintf (ftype, " 1) %-7.4f %-7.4f %d %d  #   position x, y\n", x, 
			                y, fpar->ia[1], fpar->ia[2]);
            fprintf (ftype, " 3) %-7.4f    %-5d   #     mu(Rb)\n", 
                        			    fpar->a[3], fpar->ia[3]);
            fprintf (ftype, " 4) %-7.4f    %-5d   #      Rb\n", 
						fpar->a[4], fpar->ia[4]);
            fprintf (ftype, " 5) %-7.4f    %-5d   #     alpha \n", 
		                            fpar->a[5], fpar->ia[5]);
            fprintf (ftype, " 6) %-7.4f    %-5d   #      beta \n", 
                      			    fpar->a[6], fpar->ia[6]);
            fprintf (ftype, " 7) %-7.4f    %-5d   #     gamma \n", 
			                            fpar->a[7], fpar->ia[7]);
	};

        if (strncmp (fpar->objtype, "king", 4) == 0) {

            fprintf (ftype, " 1) %-7.4f %-7.4f %d %d  #   position x, y\n", x, 
			                y, fpar->ia[1], fpar->ia[2]);
            fprintf (ftype, " 3) %-7.4f    %-5d   #     mu(0) \n", 
                        			    fpar->a[3], fpar->ia[3]);
            fprintf (ftype, " 4) %-7.4f    %-5d   #      Rc \n", 
						fpar->a[4], fpar->ia[4]);
            fprintf (ftype, " 5) %-7.4f    %-5d   #      Rt \n", 
		                            fpar->a[5], fpar->ia[5]);
            fprintf (ftype, " 6) %-7.4f    %-5d   #    alpha \n", 
                      			    fpar->a[6], fpar->ia[6]);
            fprintf (ftype, " 7) %-7.4f    %-5d   #     ----- \n", 
			                            fpar->a[7], fpar->ia[7]);
	};

        if (strncmp (fpar->objtype, "expdisk", 7) == 0) {

            fprintf (ftype, " 1) %-7.4f %-7.4f %d %d  #   position x, y\n", x, 
			                y, fpar->ia[1], fpar->ia[2]);
            fprintf (ftype, " 3) %-7.4f    %-5d   #  total magnitude \n", 
                        			    fpar->a[3], fpar->ia[3]);
            fprintf (ftype, " 4) %-7.4f    %-5d   #      Rs \n", 
						fpar->a[4], fpar->ia[4]);
            fprintf (ftype, " 5) %-7.4f    %-5d   #     ----- \n", 
		                            fpar->a[5], fpar->ia[5]);
            fprintf (ftype, " 6) %-7.4f    %-5d   #     ----- \n", 
                      			    fpar->a[6], fpar->ia[6]);
            fprintf (ftype, " 7) %-7.4f    %-5d   #     ----- \n", 
			                            fpar->a[7], fpar->ia[7]);
	};

        if (strncmp (fpar->objtype, "sersic", 6) == 0) {

            fprintf (ftype, " 1) %-7.4f %-7.4f %d %d  #   position x, y\n", x, 
			                y, fpar->ia[1], fpar->ia[2]);
            fprintf (ftype, " 3) %-7.4f    %-5d   #  total magnitude \n", 
                        			    fpar->a[3], fpar->ia[3]);
            fprintf (ftype, " 4) %-7.4f    %-5d   #      R_e\n", 
						fpar->a[4], fpar->ia[4]);
            fprintf (ftype, " 5) %-7.4f    %-5d   #  exponent (de Vaucouleurs = 4) \n", 
		                            fpar->a[5], fpar->ia[5]);
            fprintf (ftype, " 6) %-7.4f    %-5d   #     ----- \n", 
                      			    fpar->a[6], fpar->ia[6]);
            fprintf (ftype, " 7) %-7.4f    %-5d   #     ----- \n", 
			                            fpar->a[7], fpar->ia[7]);
	};

        if (strncmp (fpar->objtype, "devauc", 6) == 0 ) {

            fprintf (ftype, " 1) %-7.4f %-7.4f %d %d  #   position x, y\n", x, 
			                y, fpar->ia[1], fpar->ia[2]);
            fprintf (ftype, " 3) %-7.4f    %-5d   #  total magnitude \n", 
                        			    fpar->a[3], fpar->ia[3]);
            fprintf (ftype, " 4) %-7.4f    %-5d   #      R_e\n", 
						fpar->a[4], fpar->ia[4]);
            fprintf (ftype, " 5) %-7.4f    %-5d   #     ----- \n", 
		                            fpar->a[5], fpar->ia[5]);
            fprintf (ftype, " 6) %-7.4f    %-5d   #     ----- \n", 
                      			    fpar->a[6], fpar->ia[6]);
            fprintf (ftype, " 7) %-7.4f    %-5d   #     ----- \n", 
			                            fpar->a[7], fpar->ia[7]);
	};

	if (strncmp (fpar->objtype, "gaussian", 8) == 0) {

            fprintf (ftype, " 1) %-7.4f %-7.4f %d %d  #   position x, y \n", 
	                x, y, fpar->ia[1], fpar->ia[2]);
            fprintf (ftype, " 3) %-7.4f    %-5d   #   magnitude \n", 
        			                fpar->a[3], fpar->ia[3]);
            fprintf (ftype, " 4) %-7.4f    %-5d   #     FWHM           \n", 
                            			 fpar->a[4], fpar->ia[4]);
            fprintf (ftype, " 5) %-7.4f    %-5d   #     ----- \n", 
		                            fpar->a[5], fpar->ia[5]);
            fprintf (ftype, " 6) %-7.4f    %-5d   #     ----- \n", 
                      			    fpar->a[6], fpar->ia[6]);
            fprintf (ftype, " 7) %-7.4f    %-5d   #     ----- \n", 
			                            fpar->a[7], fpar->ia[7]);
	};

	if (strncmp (fpar->objtype, "psf", 3) == 0) {

            fprintf (ftype, " 1) %-7.4f %-7.4f %d %d  #   position x, y \n", 
	                x, y, fpar->ia[1], fpar->ia[2]);
            fprintf (ftype, " 3) %-7.4f    %-5d   #   magnitude \n", 
        			                fpar->a[3], fpar->ia[3]);
            fprintf (ftype, " 4) %-7.4f    %-5d   #     ----- \n", 
		                            fpar->a[4], fpar->ia[4]);
            fprintf (ftype, " 5) %-7.4f    %-5d   #     ----- \n", 
		                            fpar->a[5], fpar->ia[5]);
            fprintf (ftype, " 6) %-7.4f    %-5d   #     ----- \n", 
                      			    fpar->a[6], fpar->ia[6]);
            fprintf (ftype, " 7) %-7.4f    %-5d   #     ----- \n", 
			                            fpar->a[7], fpar->ia[7]);
	};

	if (strncmp (fpar->objtype, "moffat", 6) == 0 ) {

            fprintf (ftype, " 1) %-7.4f %-7.4f %d %d  #   position x, y \n", 
	                x, y, fpar->ia[1], fpar->ia[2]);
            fprintf (ftype, " 3) %-7.4f    %-5d   #   magnitude        \n", 
        			                fpar->a[3], fpar->ia[3]);
            fprintf (ftype, " 4) %-7.4f    %-5d   #     FWHM           \n", 
                            			 fpar->a[4], fpar->ia[4]);
            fprintf (ftype, " 5) %-7.4f    %-5d   #   powerlaw \n", 
						fpar->a[5], fpar->ia[5]);
            fprintf (ftype, " 6) %-7.4f    %-5d   #     ----- \n", 
                      			    fpar->a[6], fpar->ia[6]);
            fprintf (ftype, " 7) %-7.4f    %-5d   #     ----- \n", 
			                            fpar->a[7], fpar->ia[7]);
	};

	if (strncmp (fpar->objtype, "sky", 3) == 0)  {
            fprintf (ftype, " 1) %-7.4f    %-5d   #  sky background at center of fitting region  \n", fpar->a[1], fpar->ia[1]);
	    fprintf (ftype, " 2) %-7.4f    %-5d   #  dsky/dx (sky gradient in x) \n", fpar->a[2], fpar->ia[2]);
	    fprintf (ftype, " 3) %-7.4f    %-5d   #  dsky/dy (sky gradient in y) \n", fpar->a[3], fpar->ia[3]);
	} else {
            fprintf (ftype, " 8) %-7.4f    %-5d   #  axis ratio (b/a)  \n", 
                        			    fpar->a[8], fpar->ia[8]);
            fprintf (ftype, " 9) %-7.4f    %-5d   #  position angle (PA) \n", 
			                            fpar->a[9], fpar->ia[9]);
            fprintf (ftype, "10) %-7.4f    %-5d   #  diskiness(-)/boxiness(+)\n", 
			                            fpar->a[10], fpar->ia[10]);
	};

        fprintf (ftype, " Z) %-6d             #  Output option (0 = residual, 1 = Don't subtract) \n", 
			                            fpar->outtype);

        fpar = fpar->next;
        fprintf (ftype, "\n");
    };
    fprintf (ftype, "================================================================================\n");
    fprintf (ftype, "\n");

    return;
}


