
/*===============================================*/
/*                                               */
/*  Compute the 2nd derivative of 2n-vso-mu and  */
/*  store the results into Y2                    */
/*                                               */
/*===============================================*/

void spline (float *, float *, int, float, float, float *);

void getmu (float *TWOn, float *MU, float *Y2, int ninterp)
{

/* Hardwire the derivative of the 1st point (=0.05) and the last
   point (=1)                                                     */

    spline (TWOn, MU, ninterp, 0.05, 1., Y2);
    return;
}


