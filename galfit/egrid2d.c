#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define FLOAT_MIN 1.e-7
#define PI 3.14159265352982
#define ASTEP 0.01
#define DELPHI 1.0

float rdist (float x, float y, float a[]);
float egrid (float (*nuker)(float a[], float x, float y), float xmin, 
			float xmax, float ymin, float ymax, float a[]);
float xy2phi (float x, float y, float a[], int *quad);

/***************************************************************************/

float egrid (float (*func)(float a[], float x, float y), float xmin, 
			float xmax, float ymin, float ymax, float a[]) 
							
{
    void minmaxphi (float xmin, float xmax, float ymin, float ymax, float a[],
				float *minphi, float *maxphi, int quad[]);
    void minmaxrad (float xmin, float xmax, float ymin, float ymax, float a[],
						float *minr, float *maxr);
    float phi2xy (float phi, float r, float a[], float *x, float *y);
    void swap (float *min, float *max);
    float traparea (float a[], float cphi, float r);

    float xf[33], yf[3], minphi, maxphi, minr, maxr, x, y,
	  theta, r, cphi, flux, dphi, dflux, temp;
    int i, done = 0, centpix = 0, quad[5]={0, 0, 0, 0, 0};

    a[9] = -(a[9] - PI/2.);   /* A kludge, due to the definition of PA in
                                 nuker.c                                  */

    if (xmin > xmax) 
        swap (&xmin, &xmax);

    if (ymin > ymax)
        swap (&ymin, &ymax);

    minmaxphi (xmin, xmax, ymin, ymax, a, &minphi, &maxphi, quad);
    minmaxrad (xmin, xmax, ymin, ymax, a, &minr, &maxr);

    if ((xmin <= a[1] && xmax >= a[1]) && (ymin <= a[2] && ymax >= a[2]) || 
								quad[0]==1) {
        minr = 0.;
        centpix = 1;
        if ((xmin < a[1] && xmax > a[1]) && (ymin < a[2] && ymax > a[2])) {
            minphi = 0.;
            maxphi = 2. * PI;
	};
    }

    maxr = maxr / a[8];       /* Make sure we integrate out       */
    minr = minr * a[8];       /* and in far enough.               */

    if (minr == 0.) centpix = 1;  /* If the center is on the boundary */

    flux = 0.;
    dphi = DELPHI /180. * PI;
    for (cphi = minphi; cphi <= maxphi; cphi += dphi) {
        r = maxr;
        done = 0;
        while (!done &&  r >= minr) {
            r = r / (1+ASTEP);
            phi2xy (cphi, r, a, &x, &y);
            if (x >= xmin && x <= xmax && y >= ymin && y <= ymax) {
                dflux = func (a, x, y) * traparea (a, cphi, r);
                flux += dflux;
                if (centpix && (dflux/flux) < 1.e-5)  done = 1;
	    } else if (r < minr)
                done = 1;
        };
    };

    a[9] = (-a[9] + PI/2.);   /* Convert back to the regular def. of PA  */

    return (flux);

}

/***************************************************************************/

float traparea (float a[], float cphi, float r)
{
    float phi2xy (float phi, float r, float a[], float *x, float *y);

    float x[5], y[5], area, delphi;

    delphi = DELPHI / 180. * PI;

    phi2xy (cphi - delphi/2., (r*(1+ASTEP) + r)/2., a, &x[1], &y[1]);
    phi2xy (cphi + delphi/2., (r*(1+ASTEP) + r)/2., a, &x[2], &y[2]);
    phi2xy (cphi + delphi/2., (r/(1+ASTEP) + r)/2., a, &x[3], &y[3]);
    phi2xy (cphi - delphi/2., (r/(1+ASTEP) + r)/2., a, &x[4], &y[4]);

    area = 0.5 * (fabs((x[1]-x[2]) * (y[3]-y[2]) - (y[1]-y[2]) * (x[3]-x[2])) +
           fabs((x[1]-x[4]) * (y[3]-y[4]) - (y[1]-y[4]) * (x[3]-x[4])));

    return (area);

}

/***************************************************************************/

float phi2xy (float phi, float r, float a[], float *x, float *y)
{
    float cosphi, sinphi, sinPA, cosPA, theta, dx, dy, pa;

    cosphi = cos(phi);
    sinphi = sin(phi);
    sinPA = sin(a[9]);
    cosPA = cos(a[9]);

    dx = r * (-cosphi * sinPA - a[8] * sinphi * cosPA);
    dy = r * (cosphi * cosPA - a[8] * sinphi * sinPA);
    *x = dx + a[1];
    *y = dy + a[2];

    theta = atan (dy/dx);
    return (theta);
}

/***************************************************************************/

void minmaxrad (float xmin, float xmax, float ymin, float ymax, float a[],
						float *minr, float *maxr)
{
    int i;
    float rd[5];
    float cdist (float x, float y, float a[]), near;

    rd[1] = rdist (xmin, ymax, a);
    rd[2] = rdist (xmin, ymin, a);
    rd[3] = rdist (xmax, ymin, a);
    rd[4] = rdist (xmax, ymax, a);

    *minr = 1.e6;
    *maxr = 0;

    /* If the integration box doesn't cross the x or y-axis centered on
       the object... */

    for (i=1; i <= 4; i++) {
        if (rd[i] < *minr) *minr = rd[i];
        if (rd[i] > *maxr) *maxr = rd[i];
    };

    /* But if the integration box does cross the x or y-axis centered on
       the object... */

    if ((a[1] < xmax && a[1] > xmin)) {
        if (ymax - a[2] < 0.) 
	    near = ymax;
        else
            near = ymin;
        *minr = rdist (a[1], near, a);
    };

    if ((a[2] < ymax && a[2] > ymin)) {
        if (xmax - a[1] < 0.) 
	    near = xmax;
        else
            near = xmin;
        *minr = rdist (near, a[2], a);
    }

    /* If the center is on the border of the integration box */

    if (a[1] < xmax && a[1] > xmin && (a[2] == ymin || a[2] == ymax)) *minr=0;
    if (a[2] < ymax && a[2] > ymin && (a[1] == xmin || a[1] == xmax)) *minr=0;

}

/***************************************************************************/

void minmaxphi (float xmin, float xmax, float ymin, float ymax, float a[],
				float *minphi, float *maxphi, int quad[])
{
    int quadrant[5], i;
    float phi[5];

    phi[1] = xy2phi (xmin, ymax, a, &quadrant[1]);
    phi[2] = xy2phi (xmin, ymin, a, &quadrant[2]);
    phi[3] = xy2phi (xmax, ymin, a, &quadrant[3]);
    phi[4] = xy2phi (xmax, ymax, a, &quadrant[4]);

    for (i=1; i<= 4; i++)
        quad[quadrant[i]] = 1;

    *minphi = 2 * PI;
    *maxphi = 0.;

    /* gotta watch out if the integration region crosses phi = 0 */
    
    if (quad[1] && quad[4]) {
        for (i=1; i<= 4; i++) {
            if ((phi[i] - PI) >= -1.e-5 && phi[i] < *minphi && quadrant[i])
		*minphi = phi[i];
            if (phi[i] < PI && phi[i] > *maxphi && quadrant[i])
                *maxphi = phi[i];
        };
        *maxphi += (2*PI);
    } else {
        for (i=1; i<= 4; i++) {
            if (phi[i] > *maxphi && quadrant[i])
                *maxphi = phi[i];
            if (phi[i] < *minphi && quadrant[i])
                *minphi = phi[i];
        };
    };

/*    printf ("\n minphi = %f maxphi = %f\n", 
			*minphi/PI * 180, *maxphi/PI * 180);  */

}

/***************************************************************************/

float rdist (float x, float y, float a[]) 
{
    float r, c, xd, yd, Ib, pa, cosPA, sinPA, xrot, absxrot, yrot, yrot8,
	  absyrot, fabsr;

    c = a[10] + 2.;
    pa = a[9];

    xd = x - a[1];
    yd = y - a[2];

    sinPA = sin(pa);
    cosPA = cos(pa);

    xrot = -xd * sinPA + yd * cosPA;
    absxrot = pow (fabs(xrot), c);
    yrot = xd * cosPA + yd * sinPA;
    yrot8 = yrot / a[8];
    absyrot = pow (fabs(yrot8), c);
    fabsr = absxrot + absyrot;

    r = pow (fabsr, 1./c);
    return (r);
}

/***************************************************************************/

float xy2phi (float x, float y, float a[], int *quadrant)
{
    float phi, thetap, theta, xd, yd, galPA, sintheta, costheta, cosqr, l, 
	  pixPA; 

    xd = x - a[1];
    yd = y - a[2];

    galPA = a[9];		      /* PA is the PA of the galaxy with 0 deg
				  	 up, and -90 degrees to the right. */

    /* theta is the angle spanned by pos (x,y), the center of the galaxy,
       and North of the galaxy */

    if (xd == 0.) {
        if (yd < 0.)
            pixPA = PI;
        else if (yd >= 0.)
            pixPA = 0.;
    } else if (yd == 0.) {
        if (xd > 0.) 
            pixPA = -PI/2.;
        else if (xd < 0.) 
            pixPA = PI/2.;
    } else {
        l = pow ( xd * xd + yd * yd, 0.5);
        xd = xd / l;
        yd = yd / l;
        pixPA = atan(yd / xd) + PI/2.;
    };

    /* Identify which quadrant the pixel is in relative to the center
       of the object.  Adjust so that theta is between 0 and 360.      */

    if ( (xd > 0. && yd < 0.) || (xd > 0. && yd > 0.))
        pixPA = pixPA + PI;

    theta = pixPA - galPA;
    if (theta < 0.) theta = theta + 2*PI;

    if (xd == 0. && yd == 0.)
        *quadrant = 0.;
    else
        *quadrant = (int) (theta / (PI/2.) + 1);
    if (*quadrant == 5) *quadrant = 1;

    if (fabs(theta - PI/2.) < FLOAT_MIN) 
	phi = PI/2.;
    else if (fabs (theta + PI/2.) < FLOAT_MIN)
        phi = -PI/2.;
    else
        phi = atan (1/a[8] * tan(theta));

    if ( (*quadrant == 2 || *quadrant == 3) && fabs(theta - PI/2.) > FLOAT_MIN)
        phi += PI;
    if ( *quadrant == 4 ) phi += (2*PI);

/*  printf ("\n x = %2.1f  y = %2.1f  pixPA = %5.2f  theta = %5.2f  phi = %5.2f  quad = %i\n", x, y, pixPA/PI * 180.,	theta/PI * 180., phi / PI * 180., *quadrant); */
    return (phi);
}

/***************************************************************************/

void swap (float *min, float *max) 
{
    float temp;

    if (*min > *max) {
        temp = *max;
        *max = *min;
        *min = temp;
    };
}

