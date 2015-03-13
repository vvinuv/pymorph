/* note #undef's at end of file */
#define NRANSI
#include "nrutil.h"

void bcuint(float **c, float x1l, float x1u, float x2l, float x2u, 
    float x1, float x2, float *ansy, float *ansy1, float *ansy2)
{
	int i;
	float t,u,d1,d2;

	d1=x1u-x1l;
	d2=x2u-x2l;
	if (x1u == x1l || x2u == x2l) nrerror("Bad input in routine bcuint");
	t=(x1-x1l)/d1;
	u=(x2-x2l)/d2;
	*ansy=(*ansy2)=(*ansy1)=0.0;
	for (i=4;i>=1;i--) {
		*ansy=t*(*ansy)+((c[i][4]*u+c[i][3])*u+c[i][2])*u+c[i][1];
		*ansy2=t*(*ansy2)+(3.0*c[i][4]*u+2.0*c[i][3])*u+c[i][2];
		*ansy1=u*(*ansy1)+(3.0*c[4][i]*t+2.0*c[3][i])*t+c[2][i];
	}
	*ansy1 /= d1;
	*ansy2 /= d2;
}
#undef NRANSI
