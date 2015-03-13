/****************************************************************************\
*  From Numerical Recipes.  This subroutine calculates the mathematical Beta *
*  function, which is ultimately used to compute the integrated magnitude of *
*  the generalized elliptical Sersic, Gaussian, Expdisk, etc. functions.     *
*									     *
*  The generalized ellipse is defined as:                                    *
*       (|x|/a)^c + (|y|/b)^c = 1					     *
*  where a true ellipse has power of c^2.  			             *
*									     *
*  The function "ratio" returns the ratio of the generalized ellipse area    *
*  over the true ellipse with the same axis ratio.  We need to do this       *
*  because while we can determine the total flux of elliptical models easily *
*  a generalized ellipse is more complicated.                                *
\****************************************************************************/

#include <math.h>
#include "const.h"

float gammln(float xx);

float ratio (float c) 
{
    float f;
    float beta (float, float);

    f = PI * c / (4. * beta(1./c, 1+1./c));
    return (f);
}

float beta(float z, float w)
{
	return exp(gammln(z)+gammln(w)-gammln(z+w));
}

