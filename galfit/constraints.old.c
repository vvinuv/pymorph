#include <math.h>
#include "structs.h"

void constraints (int ma, int mfit, float atry[], float a[], int ia[], 
	struct cons *constr, double *da, float *aorig)

{
    float bpos, bneg, R;
    int i, j, k, l, x, y, comp1, comp2, ncomp;

    ncomp = ma / 10;
    while (constr != NULL && constr->comp[1] != 0) {

	/****************************************\
	*  Constrain deviation from input value  *
	\****************************************/

	if (constr->comp[1] <= ncomp && constr->comp[2] <= ncomp) {
	    if (constr->op == 0) {
	        bneg = constr->cval[1];
	        bpos = constr->cval[2];
	        i = (constr->comp[1]-1) * 10 + constr->par;

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

	        if (atry[i] > bpos)
	            atry[i] = bpos;

	        if (atry[i] < bneg)
		    atry[i] = bneg;

	    } else if (constr->op == '-') {

	/************************************************************\
	*  Constrain difference between parameters of two to         *
	*  different components to be within a range bpos and bneg.  *
	\************************************************************/

	        comp1 = constr->comp[1];
	        comp2 = constr->comp[2];
	        bpos = constr->cval[2];
	        bneg = constr->cval[1];

	        i = (constr->comp[1]-1) * 10 + constr->par;
	        j = (constr->comp[2]-1) * 10 + constr->par;

	        if (atry[j] > atry[i] + bpos)
		    atry[j] = atry[i] + bpos; 
	        else if (atry[j] < atry[i] + bneg)
		    atry[j] = atry[i] + bneg; 

	    } else if (constr->op == '/') {

	/************************************************\
	*  Constrain ratio R between 2 components to     *
	*  be within a range bpos and bneg               *
	\************************************************/

	        comp1 = constr->comp[1];
	        comp2 = constr->comp[2];
	        bpos = constr->cval[2];
	        bneg = constr->cval[1];

	        i = (constr->comp[1]-1) * 10 + constr->par;
	        j = (constr->comp[2]-1) * 10 + constr->par;

	        if (atry[j] < atry[i] / bpos)
		    atry[j] = atry[i] / bpos;
	        else  if (atry[j] > atry[i] / bneg)
		    atry[j] = atry[i] / bneg;
	    };
	};

	constr = constr->next;

    };
}
