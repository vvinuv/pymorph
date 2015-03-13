/***************************************************************************\
*  This subroutine is controls the sanity of the parameter steps.           *
*  Sometimes parameter steps can be way off whack even though the           *
*  sentiment is in the right direction.  This subroutine enforces changes   *
*  that are not too big.                                                    *
\***************************************************************************/

#include <math.h>
#include "nrutil.h"
#include "structs.h"

#define DELMAG 10.                /* max change in mag                 */
#define DELAX 100.                /* max change in axis ratio          */
#define MAXLENGRATIO 10.          /* max change in scale length ratio  */
#define MAXNRATIO 10.             /* maximum change in n (ratio)       */
#define MAXDPA 500.               /* maximum change in PA (degrees)    */
#define MAXDELTA 100000.          /* maximum change in position (pix)  */

#define MINAXRATIO 0.01		  /* minimum allwed axis ratio         */
#define MINLENG 0.01              /* minimum allowed scale length      */


void par_monitor (int ma, float a[], int ia[], float atry[], double da[],
    struct fitpars *fpar)
{

    extern int (*pfunc)(const char *, ...);
    extern char device[];
    int j, l;

    for (j=0, l=1; l<=ma; l++) {
	if (ia[l] == 1) {
	    ++j;
            atry[l] = a[l] + da[j];

            /***********************************************\
            *  AXIS RATIO: Constrain to be between 0 and 1  *
            \***********************************************/

            if ( (l+2) % 10 == 0) {
                if (atry[l] > 1.) {
                    /* Change the PA, but no need to change dA for PA */
                    if (ia[l+1] == 1)
                        atry[l+1] = atry[l+1] + 90.;   

                    /* Invert ax-ratio */
                    if (ia[l] == 1)
                        atry[l] = 1./atry[l];

                }

                if (fabs(a[l] / atry[l]) >= DELAX)
	            atry[l] = a[l] / DELAX;

                if (atry[l] < 0.)               /* Make sure the axis ratio */
                    atry[l] = fabs (atry[l]);   /* is not less than 0.      */
            };

            /*****************************************************\
            *  SCALELENGTHS: Constrain to be greater than 1e-2    *
            \*****************************************************/

            if ( (l+6) % 10 == 0) {

                /* Make sure length doesn't change so much at once */

                if (atry[l] / a[l] >= MAXLENGRATIO)
                    atry[l] = a[l] * MAXLENGRATIO; 
                else if (a[l] / atry[l] >= MAXLENGRATIO)
	            atry[l] = a[l] / MAXLENGRATIO;

	        if ( a[l] <= 0.00 )
	            a[l] = MINLENG;

	        if ( atry[l] <= 0.00 )
	            atry[l] = MINLENG;
            }; 

            /******************************************************\
            *  POWERLAW: Constrain alpha and n to be less than 20  *
            *  and greater than 0.01.                              *
            \******************************************************/

            if ( (l+5) % 10 == 0 ) {
		if (strncmp (fpar->objtype, "king", 4) != 0) {
                    if (fabs(atry[l] / a[l]) >= MAXNRATIO)
                        atry[l] = a[l] * MAXNRATIO;
                    else if (a[l] / atry[l] >= MAXNRATIO)
                        atry[l] = a[l] / MAXNRATIO;

                    if (atry[l] > 20.)
 	                atry[l] = 1e30;
	            else if (atry[l] < 0.0)
	                atry[l] = 1e30;
	        };
            };

            /************************************************************\
            *  PA, and (X, Y) POSITIONS: can be allowed to be negative,  *
            *  but if the change is too big, the initial guess is        *
            *  probably way too bad.  Hold it fixed for now....          *
            \************************************************************/

            if ( (l+9) % 10 == 0 || (l+8) % 10 == 0 )  {
                if (atry[l] - a[l] >= MAXDELTA)
                    atry[l] = a[l] + MAXDELTA;
                else if (a[l] - atry[l] >= MAXDELTA)
                    atry[l] = a[l] - MAXDELTA; 
            };

            /**************************************************************\
            *  DISKYNESS/BOXYNESS: constrain to between -2 and positive 2  *
            \**************************************************************/

            if (l % 10 == 0) {
	         if (atry[l] <= -2.0)
                    atry[l] = -2.0;
            };

            /********************************************\
            *  Constrain everything else to be positive  *
            \********************************************/

            if ((l+3)%10 == 0 || (l+4)%10 == 0)
                if (atry[l] < 0.) atry[l] = 0.0;
        };
	if (l%10 == 0) fpar = fpar->next;
    };
}
