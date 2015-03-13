/* note #undef's at end of file */
#define NRANSI
#include "nrutil.h"

void spline(float x[], float y[], int n, float yp1, float ypn, float y2[])
{
	int i,k;
	float p,qn,sig,un,*u;

	u=vector(1,n-1);
	if (yp1 > 0.99e30)
		y2[1]=u[1]=0.0;
	else {
		y2[1] = -0.5;
		u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
	}
	for (i=2;i<=n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
	}
	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	for (k=n-1;k>=1;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
	free_vector(u,1,n-1);
}


void splint(float xa[], float ya[], float y2a[], int n, float x, float *y)
{
	void nrerror(char error_text[]);
	int klo,khi,k;
	float h,b,a;

	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) nrerror("Bad xa input to routine splint");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}


void splie2(float x1a[], float x2a[], float **ya, int m, int n, float **y2a)
{
	void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
	int j;

	for (j=1;j<=m;j++)
		spline(x2a,ya[j],n,1.0e30,1.0e30,y2a[j]);
}


void splin2(float x1a[], float x2a[], float **ya, float **y2a, int m, int n,
	float x1, float x2, float *y)
{
	void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
	void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
	int j;
	float *ytmp,*yytmp;

	ytmp=vector(1,m);
	yytmp=vector(1,m);
	for (j=1;j<=m;j++)
		splint(x2a,ya[j],y2a[j],n,x2,&yytmp[j]);
	spline(x1a,yytmp,m,1.0e30,1.0e30,ytmp);
	splint(x1a,yytmp,ytmp,m,x1,y);
	free_vector(yytmp,1,m);
	free_vector(ytmp,1,m);
}
#undef NRANSI
