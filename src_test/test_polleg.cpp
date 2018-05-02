#include <stdio.h>				/* printf */
#include <time.h>
#include <iostream>
#include <math.h>
//using namespace std;

void polleg(int const, double const *, int const, double *);

void main()
{
	double *x, *pk, *pkx;
//	char cvar;
	int ix, ik, k, k1, nx;
	double dx, p0, p1, p2, p3, p6, p7, p10,
		    dp0, dp1, dp2, dp3, dp6, dp7, dp10,
		    x2, x3, x4, x5, x6, x7,x8,x9, x10; 
//
	nx = 2001;
	k = 1000;
	k1 = k+1;
	dx = 2.0/(nx-1);
	x  = new double [nx];
	pk = new double [(k+1)*nx];
	pkx = new double[k1];
//
	x[0] = -1.0;
	for (ix = 1; ix < nx; ix++)
		x[ix] = x[ix-1] + dx;
//
	double start = (double)clock() /(double) CLOCKS_PER_SEC;
	polleg(k, x, nx, pk);
/*	
	for (ix = 0; ix < nx; ix++) {
		polleg(k1, x[ix], pkx);
		for (ik = 0; ik < k1; ik++) {
			//printf("%i %i %12.2e \n", ix, ik, pkx[ik]);
			pk[ix*k1 + ik] = pkx[ik];
		}
	}
*/
	
	double end = (double)clock() / (double) CLOCKS_PER_SEC;

	dp0 = 0.0;
	dp1 = 0.0;
	dp2 = 0.0;
	dp3 = 0.0;
	dp6 = 0.0;
	dp7 = 0.0;
	dp10 = 0.0;

	for (ix = 0; ix < nx; ix++){
		p0 = 1.0;
		p1 = x[ix];
		x2 = p1*p1;
		x3 = x2*p1;
		x4 = x3*p1;
		x5 = x4*p1;
		x6 = x5*p1;
		x7 = x6*p1;
        x8 = x7*p1;
		x9 = x8*p1;
		x10 = x9*p1;
		p2 = 0.5*(3.0*x2 - 1.0);
		p3 = 0.5*(5.0*x3 - 3.0*x[ix]);
		p6 = (1.0/16.0)*(231.0*x6 - 315.0*x4 + 105.0*x2 - 5.0);
        p7 = (1.0/16.0)*(429.0*x7 - 693.0*x5 + 315.0*x3 - 35.0*x[ix]);
       p10 = (1.0/256.0)*(46189.0*x10 - 109395.0*x8 + 90090.0*x6 - 30030.0*x4 + 3465.0*x2 - 63.0);
//
		dp0 += fabs(p0 - pk[ix*k1+0]);
		dp1 += fabs(p1 - pk[ix*k1+1]);
		dp2 += fabs(p2 - pk[ix*k1+2]);
		dp3 += fabs(p3 - pk[ix*k1+3]);
		dp6 += fabs(p6 - pk[ix*k1+6]);
		dp7 += fabs(p7 - pk[ix*k1+7]);
		dp10 += fabs(p10 - pk[ix*k1+10]);
	}
	printf("%12.2e %12.2e %12.2e %12.2e %12.2e %12.2e %12.2e \n", dp0/nx, dp1/nx, dp2/nx, dp3/nx, dp6/nx, dp7/nx, dp10/nx);
//
	printf("\n");
	printf("cpu time %fs\n", end - start);
	printf("Done! \n");
//	std::cin >> cvar;
//
}