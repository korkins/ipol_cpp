#include <stdio.h>				/* printf */
#include <time.h>
#include <iostream>
#include <math.h>
#define PI 3.141592653589793
//using namespace std;

void gauszw(double const, double const, int const, double *, double *);
void print2d(double **, int const, int const, char *);

void main()
{
	double *z, *w, **aout;
	int i, n, nlin, ncol, ilin;
	double x1, x2;
//
	n = 13;
	z  = new double [n];
	w = new double [n];
//
	nlin = n;
	ncol = 2+3;
//
	aout = new double *[nlin];
	for (ilin = 0; ilin < nlin; ilin++)
		aout[ilin] = new double[ncol];
//
	x1 = 0.0;
	x2 = PI;
	gauszw(x1, x2, n, z, w);
//
	for (i = 0; i < n; i++) {
		aout[i][0] = i;
		aout[i][1] = x1;
		aout[i][2] = x2;
		aout[i][3] = z[i];
		aout[i][4] = w[i];

	}
//
	print2d(aout, nlin, ncol, "gauszw0pi.txt");
	printf("Done! \n");
//
}