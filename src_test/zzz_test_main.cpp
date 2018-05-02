#include <stdio.h>				/* printf */
#include <time.h>
#include <iostream>
#include <math.h>
//using namespace std;

double sum1d(const double *, int const);

void main()
{
	int ilin, icol;
	double x, **y2d, *y1d;
//
	y1d = new double[3];
//	3x4 array
	y2d = new double *[3]; // ilin
	for (ilin = 0; ilin < 3; ilin++) {
		y1d[ilin] = 1.0*ilin;
		y2d[ilin] = new double[4];
		for (icol = 0; icol < 4; icol++)
			y2d[ilin][icol] = 1.0*(ilin + icol);
		x = sum1d(y2d[ilin], 4);
		printf("i, x = %i, %f \n", ilin, x);
	}
	x = sum1d(y1d, 3);
	printf("sum1(1d) %f \n", x);
	printf("Done! \n");
//
}

double sum1d(const double *x1d, int const nx) {
	double s = 0.0;
	for (int i = 0; i < nx; i++) s += x1d[i];
	return s;
}