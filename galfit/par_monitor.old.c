/***************************************************************************\
*  This subroutine is controls the sanity of the parameter steps.           *
*  Sometimes parameter steps can be way off whack even though the           *
*  sentiment is in the right direction.  This subroutine enforces changes   *
*  that are not too big.                                                    *
\***************************************************************************/

#include <math.h>
#include "nrutil.h"
#if DEBUG
    #include <curses.h>
#endif

#define DELMAG 10.                 /* max change in mag                 */
#define DELAX 10.                 /* max change in axis ratio          */
#define MAXLENGRATIO 10.          /* max change in scale length ratio  */
#define MAXNRATIO 10.		  /* maximum change in n (ratio)       */
#define MAXDPA 500.		  /* maximum change in PA (degrees)    */
#define MAXDELTA 10.		  /* maximum change in position (pix)  */

#define MINAXRATIO 0.1		  /* minimum allwed axis ratio         */
#define MINLENG 0.1               /* minimum allowed scale length      */


void par_monitor (int l, int j, float a[], int ia[], float atry[], double da[])
{

    extern int (*pfunc)(const char *, ...);
    extern char device[];

    atry[l] = a[l] + da[j];

#if DEBUG 
    pfunc ("%d atry[%d] = %f,  da = %f\n", (int)(l/10)+1, l%10, a[l], da[j]);
    if (strncmp (device, "regular", 7) != 0)
        refresh();
#endif 

    /***********************************************\
    *  AXIS RATIO: Constrain to be between 0 and 1  *
    \***********************************************/

    if ( (l+2) % 10 == 0) {
        if (atry[l] > 1. && atry[l] < 10.) {
            /* Change the PA, but no need to change dA for PA */
            if (ia[l+1] == 1)
                atry[l+1] = atry[l+1] + 90.;   

            /* Invert ax-ratio */
            if (ia[l] == 1)
                atry[l] = 1./atry[l];

            /* Change Rs, but no need to change dA for Rs */
            if (ia[l-4] == 1)
                atry[l-4] = atry[l-4]/atry[l]; 

        } else if (atry[l] > 10.)
            atry[l] = a[l];         /* Change is way to big... hold */
                                    /* it fixed for now....         */

        if (atry[l] <= 0. || a[l] / atry[l] >= DELAX)
	    atry[l] = a[l] / DELAX;

        if (atry[l] <= MINAXRATIO)        /* Make sure the axis ratio */
            atry[l] = MINAXRATIO;         /* is not less than 0.01    */

        return;
    };

    /****************************************************\
    *  MAGNITUDE:  Constrain to be between -50 and +50   *
    \****************************************************/

    if ( (l+7) % 10 == 0) {

         /* If the next step is way too big big, damp out change */

         if (atry[l] - a[l] >= DELMAG)
	     atry[l] = a[l] + DELMAG;
         else if (a[l] - atry[l] >= DELMAG)
             atry[l] = a[l] - DELMAG;

	 if (atry[l] >= 50.)
	     atry[l] = 50.;
	 else if (atry[l] <= -50.)
	     atry[l] = -50.;

         return;
    }; 

    /*****************************************************\
    *  SCALELENGTHS: Constrain to be greater than 1e-2    *
    \*****************************************************/

    if ( (l+6) % 10 == 0) {

       /* Make sure length doesn't change so much at once */

/*      if (atry[l] / a[l] >= MAXLENGRATIO)
            atry[l] = a[l] * MAXLENGRATIO; 
        else if (a[l] / atry[l] >= MAXLENGRATIO)
	    atry[l] = a[l] / MAXLENGRATIO;
*/
	if ( a[l] <= 0.1 ) a[l] = 0.1;

        return;
    }; 

    /******************************************************\
    *  POWERLAW: Constrain alpha and n to be less than 20  *
    *  and greater than 0.01.                              *
    \******************************************************/

    if ( (l+5) % 10 == 0 ) {
        if (fabs(atry[l] / a[l]) >= MAXNRATIO)
            atry[l] = a[l] * MAXNRATIO;
        else if (atry[l] <= 0. || a[l] / atry[l] >= MAXNRATIO)
            atry[l] = a[l] / MAXNRATIO;

        if (atry[l] > 20.)
	    atry[l] = 20.;
        else if (atry[l] < 0.01)
	    atry[l] = 0.01;

        return;
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

        return;
    };

    if ( (l+1) % 10 ==0 ) {
        if (fabs(da[j]) >= MAXDPA)
            atry[l] = a[l];

        return;
    };

    /**************************************************************\
    *  DISKYNESS/BOXYNESS: constrain to between -2 and positive 2  *
    \**************************************************************/

    if (l % 10 == 0) {
        if (atry[l] <= -2.0)
            atry[l] = -1.99;

        return;
    };

    /********************************************\
    *  Constrain everything else to be positive  *
    \********************************************/

    if ( ((l+3)%10 == 0 || (l+4)%10 == 0) && a[l] < 0.) a[l] = 0.01;

    if ( ((l+3)%10 == 0 || (l+4)%10 == 0) && atry[l] <= 0.) {
        atry[l] = 0.01;
    }

    return;

}
