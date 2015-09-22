#include <math.h>
#include "mymath.h"
#include "nrutil.h"
#include "structs.h"
#include "const.h"

#define KERNEL 30

struct bicubic_grid **cg;

double **psfprep (struct image psf, long model_naxes[], int add0pad,
    long complex_naxes[]);
void bcuint(float **c, float x1l, float x1u, float x2l, float x2u, 
    float x1, float x2, float *ansy, float *ansy1, float *ansy2);
void convolve (struct image *model, double **fftpsf, long complex_naxes[]);

float alpha_win, tau_win;

void shift_psf (struct image psf, struct image *tpsf, float dx, float dy)
{

    float sinc (float x);
    extern float fwhmpsf;
    float ansy1, ansy2;
    int i, j, xn, yn;
    float **sincfunc, xc, yc, x, y, radius;
    float yi[5], y1i[5], y2i[5], y12i[5];
    double **fftsinc;
    long complex_naxes[3];
    struct image sinc_img;

    xc = yc = NINT(KERNEL / 2.+0.5);
    radius = fwhmpsf * 2.5;

    if (fwhmpsf >= 1.2 && fwhmpsf <= 2.2) {
        alpha_win = 1;
        tau_win = 10.;
    } else if (fwhmpsf > 2.2 && fwhmpsf <= 4) {
	alpha_win = 2.;
	tau_win = 20.;
    };

    if (fwhmpsf >= 1.2 && fwhmpsf <= 4) {

        /***************************\
        *  Duplicate the input PSF  *
        \***************************/

        for (j=1; j <= psf.naxes[2]; j++) 
            for (i=1; i <= psf.naxes[1]; i++)
	        tpsf->z[j][i] = psf.z[j][i];

        /***************************************************************\
        *  Create a 2-D Sinc function, tapered by the Kaiser function,  *
        *  and shifted by dx and dy.                                    *
        \***************************************************************/

        sincfunc = matrix (1, KERNEL, 1, KERNEL);

        for (j=1 ; j <= KERNEL; j++)
            for (i=1; i <= KERNEL; i++) {
	        x = i - xc - dx;
	        y = j - yc - dy;
	        sincfunc[j][i] = sinc (x) * sinc (y);
	    };

        /****************************************************************\
        *  Now convolve the shifted, tapered Sinc function with the psf  *
        \****************************************************************/

        sinc_img.naxes[1] = KERNEL;
        sinc_img.naxes[2] = KERNEL;
        sinc_img.z = sincfunc;
        fftsinc = psfprep (sinc_img, psf.naxes, 1, complex_naxes);
        convolve (tpsf, fftsinc, complex_naxes);

        free_dmatrix (fftsinc, 0, complex_naxes[2]-1, 0, complex_naxes[1]*2-1);
        free_matrix (sincfunc, 1, KERNEL, 1, KERNEL);

        /*************************************************\
        *  Now bi-cubic interpolate to get outer regions  *
        \*************************************************/

        xc = (int)(psf.naxes[1]/2.) + 1;
        yc = (int)(psf.naxes[2]/2.) + 1;
        for (j=1; j <= psf.naxes[2]; j++)
            for (i=1; i <= psf.naxes[1]; i++) {
	        x = i - xc + dx + 1;
	        y = j - yc + dy + 1;
	        if (sqrt(x*x+y*y) >= radius) {
		    xn = (int)(i-dx);
		    yn = (int)(j-dy);

		    if (yn >= 0 && xn >= 0 && xn <= psf.naxes[1] &&
			yn <= psf.naxes[2])

		        bcuint (cg[yn][xn].c, xn, xn+1, yn, yn+1, 
			    i-dx, j-dy, &tpsf->z[j][i], &ansy1, &ansy2);
		    else
			tpsf->z[j][i] = 0.;
		
		};
	};
    } else {
        for (j=1 ; j <= psf.naxes[2]; j++)
            for (i=1; i <= psf.naxes[1]; i++) {
	        xn = (int)(i-dx);
	        yn = (int)(j-dy);


		if (yn >= 0 && xn >= 0 && xn <= psf.naxes[1] && 
		    yn <= psf.naxes[2])

	            bcuint (cg[yn][xn].c, xn, xn+1, yn, yn+1, 
		        i-dx, j-dy, &tpsf->z[j][i], &ansy1, &ansy2);
		else
		    tpsf->z[j][i] = 0.;
	    };
    };
}

/**************************************************************************\
*  This is a Sinc function tapered by the Kaiser function.                 *
\**************************************************************************/

float sinc (float x)
{
    float func, win;
    float bessi0(float x);

    if (x==0.) return (1.);
    if (fabs(x) > tau_win) return (0.);
    win = alpha_win * pow (1-x/tau_win * x/tau_win ,0.5);
    func = sin (x * PI)/(x * PI) * bessi0 (win) / bessi0 (alpha_win); 

    return (func);
}

/**************************************************************************/

float bessi0(float x)
{
	float ax,ans;
	double y;

	if ((ax=fabs(x)) < 3.75) {
		y=x/3.75;
		y*=y;
		ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
			+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
	} else {
		y=3.75/ax;
		ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
			+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
			+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
			+y*0.392377e-2))))))));
	}
	return ans;
}

