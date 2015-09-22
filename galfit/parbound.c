/*************************************************************************\
*  This subroutine checks to see if parameters in atry[] have exceeded    *
*  the constraint boundaries.  If so, we have to hold those parameters    * 
*  fixed and change the tia[] flag.  The tia[] flag is analogous to the   *
*  ia[] flags, and tia[] will be used to help construct a smaller         *
*  covariance matrix, from which new da[] will be calculated.             *
\*************************************************************************/

#include "structs.h"

int parbound (struct cons *constr, float atry[], int tia[], float *aorig, 
							int ma, int mfit)
{
    int tmfit, ncomp, i, j;
    float bpos, bneg;

    ncomp = ma / 10;
    tmfit = mfit;
    while (constr != NULL && constr->comp[1] != 0) {

	if (constr->comp[1] <= ncomp && constr->comp[2] <= ncomp) {

	/****************************************\
	*  Constrain deviation from input value  *
	\****************************************/

	    if (constr->op == 0) {
	        bneg = constr->cval[1];
	        bpos = constr->cval[2];
	        i = (constr->comp[1]-1) * 10 + constr->par;
/*	        j = track[constr->par][constr->comp[1]]; */

		/******************************************************\
		*  If the parameter is pegged, recalculate gradient    *
                *  with that parameter held fixed to the pegged value  *
		\******************************************************/

		if (atry[i] == aorig[i] + bpos || atry[i] == aorig[j] + bneg) {
		    tia[i] = 0;
		    --tmfit;
		};

	        if (atry[i] > aorig[i] + bpos)
	            atry[i] = aorig[i] + bpos;

	        if (atry[i] < aorig[i] + bneg)
		    atry[i] = aorig[i] + bneg;
	    } else if (constr->op == 't') {

	/*****************************************\
	*  Constrain parameter by absolute range  *
	\*****************************************/

	        bneg = constr->cval[1];
	        bpos = constr->cval[2];
	        i = (constr->comp[1]-1) * 10 + constr->par;
/*	        j = track[constr->par][constr->comp[1]]; */

		/******************************************************\
		*  If the parameter is pegged, recalculate gradient    *
                *  with that parameter held fixed to the pegged value  *
		\******************************************************/

		if (atry[i] == bpos || atry[i] == bneg) {
		    tia[i] = 0;
		    --tmfit;
		};

	        if (atry[i] > bpos)
	            atry[i] = bpos;

	        if (atry[i] < bneg)
		    atry[i] = bneg;
	    } else if (constr->op == '-') {

	/************************************************************\
	*  Constrain difference between parameters of two to         *
	*  different components to be within a range bpos and bneg.  *
	\************************************************************/

/*	        comp1 = constr->comp[1];
	        comp2 = constr->comp[2]; */
	        bpos = constr->cval[2];
	        bneg = constr->cval[1];
/*	        x = track[constr->par][comp1];
	        y = track[constr->par][comp2]; */

	        i = (constr->comp[1]-1) * 10 + constr->par;
	        j = (constr->comp[2]-1) * 10 + constr->par;

		/********************************************************\
		*  If the parameter difference is pegged, recalculate    *
                *  the *joint* gradient by holding the second parameter  *
                *  fixed.                                                *
		\********************************************************/

		if (atry[j] == aorig[i] + bpos || atry[j] == aorig[i] + bneg) {
		    tia[j] = 0;
		    --tmfit;
		}

	        if (atry[j] > atry[i] + bpos)
		    atry[j] = atry[i] + bpos; 
	        else if (atry[j] < atry[i] + bneg)
		    atry[j] = atry[i] + bneg; 
	    } else if (constr->op == '/') {

	/********************************************\
	*  Constrain ratio between 2 components to   *
	*  be within a range bpos and bneg           *
	\********************************************/

/*	        comp1 = constr->comp[1];
	        comp2 = constr->comp[2]; */
	        bpos = constr->cval[2];
	        bneg = constr->cval[1];
/*	        x = track[constr->par][comp1];
	        y = track[constr->par][comp2];  */

	        i = (constr->comp[1]-1) * 10 + constr->par;
	        j = (constr->comp[2]-1) * 10 + constr->par;

		/******************************************************\
		*  If the parameter ratio is pegged, recalculate the   *
                *  *joint* gradient by holding the second parameter    *
                *  fixed.                                              *
		\******************************************************/

		if (atry[j] == aorig[i] / bpos || atry[j] == aorig[i] / bneg) {
		    tia[j] = 0;
		    --tmfit;
		};

	        if (atry[j] < atry[i] / bpos)
		    atry[j] = atry[i] / bpos;
	        else  if (atry[j] > atry[i] / bneg)
		    atry[j] = atry[i] / bneg;
	    }
	};

	constr = constr->next;

    };

    return (tmfit);

}
