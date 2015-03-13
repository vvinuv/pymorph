#include <stdio.h>
#include <math.h>
#include "structs.h"

/* #define CSQ(i) (sqrt(chi2nu * covar[(i)][(i)]))  */


/**********************************************************************\
*  The parameter uncertainty SIG(i) is scaled by chi2/nu.              *
\**********************************************************************/

#define SIG(i) (sqrt(chi2nu) * sig[i])

float pa (float instrpa);

void fprintpar (int lowx, int lowy, struct fitpars *fpar, FILE *ftype, 
	char newout[], double *sig, float *coerr, double chisq, int ndof)
{
    extern struct inpars input;
    extern float xskycent, yskycent;
    float x, y, chi2nu;
    int noff = 0;

    chi2nu = chisq / ndof;

    fprintf (ftype, "-----------------------------------------------------------------------------\n\n");
    fprintf (ftype, "Input image     : %s%s \n", input.data, input.imgsect);
    fprintf (ftype, "Init. par. file : %s \n", input.initparfile);
    fprintf (ftype, "Restart file    : %s \n", newout);
    fprintf (ftype, "Output image    : %s \n\n", input.output);

    while (fpar != NULL) {
        x = fpar->a[1] + lowx - 1.;
	y = fpar->a[2] + lowy - 1.;
        fpar->a[9] = pa(fpar->a[9]);     /* make sure  -90 < PA <= 90 */
        
        if (strncmp(fpar->objtype, "sersic", 6)==0 || 
                                    strncmp(fpar->objtype, "moffat", 6)==0) {
            fprintf (ftype, " %-9s: (%-.2f, %-.2f)   %-7.4f   %-7.4f   %-7.4f   %-7.4f  %-7.4f  %-7.4f\n",
		fpar->objtype, x, y, fpar->a[3], fpar->a[4], fpar->a[5],
				 fpar->a[8], fpar->a[9], fpar->a[10]);

            fprintf (ftype, " %-12s (%-.2f, %-.2f)   %-7.4f   ", "", 
				SIG(1+noff),SIG(2+noff), SIG(3+noff));
	    fprintf (ftype, "%-7.4f   %-7.4f   %-7.4f  %-7.4f  %-7.4f\n",
		    SIG(4+noff), SIG(5+noff), SIG(8+noff), SIG(9+noff), SIG(10+noff));
	};

        if (strncmp(fpar->objtype, "devauc", 6)==0 || 
				strncmp(fpar->objtype, "expdisk", 7)==0) {
            fprintf (ftype, " %-9s: (%-.2f, %-.2f)   %-8.4f  %-8.4f  %-8.4f  %-8.4f  %-8.4f\n",
		fpar->objtype, x, y, fpar->a[3], fpar->a[4], 
				fpar->a[8], fpar->a[9], fpar->a[10]);

            fprintf (ftype, " %-12s (%-.2f, %-.2f)   %-8.4f  %-8.4f  %-8.4f  %-8.4f  %-8.4f\n",
			 "", SIG(1+noff), SIG(2+noff), SIG(3+noff), 
			SIG(4+noff), SIG(8+noff), SIG(9+noff), SIG(10+noff));
	};

        if (strncmp(fpar->objtype, "nuker", 5)==0) {
            fprintf (ftype, " %-9s: (%-.2f, %-.2f)  %-6.4f  %-8.4f  %-6.4f %-6.4f %-6.4f  %-6.4f  %-7.4f %-7.4f\n",
	        fpar->objtype, x, y, fpar->a[3], fpar->a[4], fpar->a[5], 
                fpar->a[6], fpar->a[7], fpar->a[8], fpar->a[9], fpar->a[10]);

            fprintf (ftype, " %-12s (%-.2f, %-.2f)  %-6.4f  %-8.4f %-6.4f  %-6.4f  %-6.4f  %-6.4f  %-7.4f  %-7.4f\n",
	            "", SIG(1+noff), SIG(2+noff), SIG(3+noff), SIG(4+noff), 
			SIG(5+noff), SIG(6+noff), SIG(7+noff), SIG(8+noff), 
						  SIG(9+noff), SIG(10+noff));
        };

        if (strncmp(fpar->objtype, "king", 4)==0) {
            fprintf (ftype, " %-9s: (%8.4f, %8.4f)  %7.4f  %8.4f  %8.4f %7.4f %7.4f %8.4f %8.4f\n",
	        fpar->objtype, x, y, fpar->a[3], fpar->a[4], fpar->a[5], 
                fpar->a[6], fpar->a[8], fpar->a[9], fpar->a[10]);

            fprintf (ftype, " %-12s (%-.2f, %-.2f)     %-8.4f  %-8.4f  %-7.4f %-7.4f %-7.4f  %-7.4f  %-7.4f\n",
			 "", SIG(1+noff), SIG(2+noff), SIG(3+noff), SIG(5+noff), 
			SIG(6+noff), SIG(8+noff), SIG(9+noff), SIG(10+noff));
	};

        if (strncmp(fpar->objtype, "gaussian", 8)==0) {
            fprintf (ftype, " %-9s: (%-.2f, %-.2f)   %-8.4f   %-8.4f   %-8.4f   %-8.4f  %-8.4f\n",
	        fpar->objtype, x, y, fpar->a[3], fpar->a[4], fpar->a[8], 
					fpar->a[9], fpar->a[10]);

            fprintf (ftype, " %-12s (%-.2f, %-.2f)   %-8.4f   %-8.4f   %-8.4f   %-8.4f  %-8.4f\n",
	    		"", SIG(1+noff), SIG(2+noff), SIG(3+noff), SIG(4+noff),
				    SIG(8+noff), SIG(9+noff), SIG(10+noff));
        };


        if (strncmp(fpar->objtype, "psf", 3)==0) {
            fprintf (ftype, " %-9s: (%-.2f, %-.2f)   %-8.4f   \n",
	        fpar->objtype, x, y, fpar->a[3]);

            fprintf (ftype, " %-12s (%-.2f, %-.2f)   %-8.4f\n",
	    		"", SIG(1+noff), SIG(2+noff), SIG(3+noff));
        };

	if (strncmp(fpar->objtype, "penny", 5)==0) {
            fprintf (ftype, " %-9s: (%-.2f, %-.2f)   %-8.4f   %-8.4f  %-8.4f\n",
				fpar->objtype, x, y, 
					fpar->a[3], fpar->a[4], fpar->a[5]);

            fprintf (ftype, " %-12s (%-.2f, %-.2f)   %-8.4f   %-8.4f  %-8.4f\n",
			"", SIG(1+noff), SIG(2+noff), SIG(3+noff), 
						SIG(4+noff), SIG(5+noff));
	};

	if (strncmp(fpar->objtype, "sky", 3) ==0) {
	    fprintf (ftype, " sky      : [%-.2f, %-.2f]   %-7.4f   %-7.4e   %-7.4e\n", 
				xskycent + lowx - 1, yskycent + lowy - 1, 
					fpar->a[1], fpar->a[2], fpar->a[3]);
            fprintf (ftype, "                             %-7.4f   %-7.4e   %-7.4e\n",
					SIG(1+noff), SIG(2+noff), SIG(3+noff));
	};

	fpar = fpar->next;
	noff += NPARS;
    };

    fprintf (ftype, " Chi^2 = %.7f,  ndof = %d\n", chisq, ndof);
    fprintf (ftype, " Chi^2/nu = %.5f \n\n", chisq / ndof);
    fprintf (ftype, "-----------------------------------------------------------------------------\n");

    return;
}

