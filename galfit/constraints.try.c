#include <math.h>
#include "structs.h"

int constraints (int ma, int mfit, int *tmfit, float atry[], float a[],
	int ia[], int tia[], struct cons *constr, double beta[], float *aorig)

{
    float bpos, bneg;
    int i, j, k, l, comp1, comp2, ncomp, hitwall=0;

    ncomp = ma / 10;
    while (constr != NULL && constr->comp[1] != 0) {
	if (constr->comp[1] <= ncomp && constr->comp[2] <= ncomp) {

	    if (constr->op == '-') {

	/************************************************************\
	*  Constrain difference between parameters of two to         *
	*  different components to be within a range bpos and bneg.  *
	\************************************************************/

	        comp1 = constr->comp[1];
	        comp2 = constr->comp[2];
	        bpos = constr->crange[2];
	        bneg = constr->crange[1];

	        i = (constr->comp[1]-1) * 10 + constr->par;
	        j = (constr->comp[2]-1) * 10 + constr->par;

	        if (atry[j] > atry[i] + bpos) {

		    /*  Deal with discontinuity in the transition   */
		    /*  between	no constraints and constraints      */

		    if (fabs(atry[i] - a[i]) < fabs(atry[j] - a[j])) {
			atry[i] = a[i];
		        atry[j] = atry[i] + bpos; 
		    } else {
			atry[j] = a[j];
			atry[i] = atry[j] - bpos;
		    };

	        } else if (atry[j] < atry[i] + bneg) {
		    if (fabs(atry[i] - a[i]) < fabs(atry[j] - a[j])) {
			atry[i] = a[i];
		        atry[j] = atry[i] + bneg; 
		    } else {
			atry[j] = a[j];
			atry[i] = atry[j] - bneg; 
		    };
		};

	        if (atry[j] == atry[i] + bpos || atry[j] == atry[i] + bneg) {
		    constr->hitwall = 1;		
		    tia[j] = 0;
		    hitwall = 1;
		} else
		    constr->hitwall = 0;

	    } else if (constr->op == '/') {

	/************************************************\
	*  Constrain ratio R between 2 components to     *
	*  be within a range bpos and bneg               *
	\************************************************/

	        comp1 = constr->comp[1];
	        comp2 = constr->comp[2];
	        bpos = constr->crange[2];
	        bneg = constr->crange[1];

	        i = (constr->comp[1]-1) * 10 + constr->par;
	        j = (constr->comp[2]-1) * 10 + constr->par;

	        if (atry[j] < atry[i] / bpos) {

		    /*  Deal with discontinuity in the transition   */
		    /*  between	no constraints and constraints      */

		    if (fabs(atry[i] - a[i]) < fabs(atry[j] - a[j])) {
			atry[i] = a[i];
		        atry[j] = atry[i] / bpos;
		    } else {
			atry[j] = a[j];
		        atry[i] = atry[j] * bpos;
		    };
		    constr->val = bpos;

	        } else  if (atry[j] > atry[i] / bneg) {
		    if (fabs(atry[i] - a[i]) < fabs(atry[j] - a[j])) {
			atry[i] = a[i];
		        atry[j] = atry[i] / bneg;
		    } else {
			atry[j] = a[j];
		        atry[i] = atry[j] * bneg;
		    };
		    constr->val = bneg;
		};

		if (atry[j] == atry[i] / bpos || atry[j] == atry[i] / bneg) {
		    constr->hitwall = 1;		
		    constr->hitwall = 1;
		    tia[j] = 0;
		    hitwall = 1;
		} else
		    constr->hitwall = 0;

	    } else if (constr->op == 0) {

	/****************************************\
	*  Constrain deviation from input value  *
	\****************************************/
	    
	        bneg = constr->crange[1];
	        bpos = constr->crange[2];
	        i = (constr->comp[1]-1) * 10 + constr->par;

	        if (atry[i] > aorig[i] + bpos)
	            atry[i] = aorig[i] + bpos;

	        if (atry[i] < aorig[i] + bneg)
		    atry[i] = aorig[i] + bneg;

		if (atry[i] == aorig[i] + bpos || atry[i] == aorig[i] + bneg) {
		    constr->hitwall = 1;
		    tia[i] = 0;
		    hitwall = 1;
		} else
		    constr->hitwall = 0;

	    } else if (constr->op == 't') {

	/*****************************************\
	*  Constrain parameter by absolute range  *
	\*****************************************/

	        bneg = constr->crange[1];
	        bpos = constr->crange[2];
	        i = (constr->comp[1]-1) * 10 + constr->par;

	        if (atry[i] > bpos)
	            atry[i] = bpos;

	        if (atry[i] < bneg)
		    atry[i] = bneg;

		if (atry[i] == bpos || atry[i] == bneg) {
		    constr->hitwall = 1;
		    tia[i] = 0;
		    hitwall = 1;
		} else
		    constr->hitwall = 0;
	    }
	};
	constr = constr->next;
    };

    return (hitwall);
}


/*****************************************************************************/

void couple_params (float *a, struct cons *constr)
{
    int comp, par, i;

    while (constr != NULL) {
	i = ((constr->comp[2])-1)*10 + constr->par;

	if (constr->hitwall) {
	    if (constr->op == '/')
	        a[i] =  constr->comp[1] / constr->val;
	    else if (constr->op == '-')
	        a[i] =  constr->comp[1] - constr->val;
	};

	constr = constr->next;
    };
}
