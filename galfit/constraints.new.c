#include <math.h>
#include "structs.h"

void constraints (int ma, int mfit, float atry[], float a[], int ia[], 
	struct cons *constr, double *da, float *aorig)

{
    float hibound, lobound, R;
    int i, j, k, l, x, y, comp1, comp2, ncomp;

    ncomp = ma / 10;
    while (constr != NULL && constr->comp[1] != 0) {

	/****************************************\
	*  Constrain deviation from input value  *
	\****************************************/

	if (constr->comp[1] <= ncomp && constr->comp[2] <= ncomp) {
	    if (constr->op == 0) {
	        lobound = constr->cval[1];
	        hibound = constr->cval[2];
	        i = (constr->comp[1]-1) * 10 + constr->par;

	        if (atry[i] > aorig[i] + hibound)
	            atry[i] = aorig[i] + hibound;

	        if (atry[i] < aorig[i] + lobound)
		    atry[i] = aorig[i] + lobound;

	    } else if (constr->op == 't') {

	/*****************************************\
	*  Constrain parameter by absolute range  *
	\*****************************************/

	        lobound = constr->cval[1];
	        hibound = constr->cval[2];
	        i = (constr->comp[1]-1) * 10 + constr->par;

	        if (atry[i] > hibound)
	            atry[i] = hibound;

	        if (atry[i] < lobound)
		    atry[i] = lobound;

	    } else if (constr->op == '-') {

	/************************************************************\
	*  Constrain difference between parameters of two different  *
	*  components to be within a range hibound and lobound.      *
	\************************************************************/

	        comp1 = constr->comp[1];
	        comp2 = constr->comp[2];
	        hibound = constr->cval[2];
	        lobound = constr->cval[1];

	        i = (constr->comp[1]-1) * 10 + constr->par;
	        j = (constr->comp[2]-1) * 10 + constr->par;

	        if (atry[i] - atry[j] >= hibound) {
		    if (fabs(atry[i] - a[i]) <= fabs(atry[j] - a[j]))
			atry[i] = atry[j] + hibound;
		    else
		        atry[j] = atry[i] - hibound; 
	        } else if (atry[i] - atry[j] <= lobound) {
		    if (fabs(atry[i] - a[i]) <= fabs(atry[j] - a[j]))
			atry[i] = atry[j] + lobound;
		    else
		        atry[j] = atry[i] - lobound; 
		};

	    } else if (constr->op == '/') {

	/************************************************\
	*  Constrain ratio R between 2 components to     *
	*  be within a range hibound and lobound         *
	\************************************************/

	        comp1 = constr->comp[1];
	        comp2 = constr->comp[2];
	        hibound = constr->cval[2];
	        lobound = constr->cval[1];

	        i = (constr->comp[1]-1) * 10 + constr->par;
	        j = (constr->comp[2]-1) * 10 + constr->par;

	        if (atry[i] / atry[j] >= hibound) {
		    if (fabs(atry[i] - a[i]) <= fabs(atry[j] - a[j]))
			atry[i] = fabs(atry[j]) * hibound;
		    else
		        atry[j] = fabs(atry[i]) / hibound;
	        } else  if (atry[i] / atry[j] <= lobound) {
		    if (fabs(atry[i] - a[i]) <= fabs(atry[j] - a[j]))
		        atry[i] = fabs(atry[j]) * lobound;
		    else
		        atry[j] = fabs(atry[i]) / lobound;
		};
	    };
	};

	constr = constr->next;

    };
}
