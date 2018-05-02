#include <stdio.h>				/* printf */
#include <time.h>
#include <iostream>
#include <math.h>
//using namespace std;

void zrssg0(double const, double *, int const, double *, int const,
	            double *, double *);

void main()
{
	double *x, *x0, *zr1, *zr2;
	int ix, ix0, nx, nx0;
	double dx, dx0, df;
//
	df = 0.03;
//
//  vza
	nx = 11;
	dx = 2.0/(nx-1);
	x  = new double [nx];
	x[0] = -1.0;
	for (ix = 1; ix < nx; ix++)
		x[ix] = x[ix-1] + dx;
//
//	sza
	nx0 = 21;
	dx0 = 1.0/(nx0-1);
	x0 = new double [nx0];
	x0[0] = 0.0;
	for (ix0 = 1; ix0 < nx0; ix0++)
		x0[ix0] = x0[ix0-1] + dx0;
//
	zr1 = new double [nx*nx0];
	zr2 = new double [nx*nx0];
//
	for (ix0 = 0; ix0 < nx0; ix0++) {
		zrssg0(df, x0, nx0, x, nx, zr1, zr2);
		for (ix = 0; ix < nx; ix++)
			printf("%12.6f", zr1[ix0*nx+ix]);
		printf("\n");
		for (ix = 0; ix < nx; ix++)
			printf("%12.6f", zr2[ix0*nx+ix]);
		printf("\n");
	}
//
	printf("Done! \n");
//	std::cin >> cvar;
//
}