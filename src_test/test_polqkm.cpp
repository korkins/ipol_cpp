#include <stdio.h>				/* printf */
#include <time.h>
#include <iostream>
#include <math.h>
//using namespace std;

void polqkm(int const, int const, double *, int const, double *);

void main()
{
	double *x, *qk1, *qk2, *qk3, *z, *q256;
//	char cvar;
	int ix, ik, k, k1, nx;
	double x2, dx, q11, q21, q31, q22, q32, q42, q33, q43,
		   dp1, dp2, dp3, dp4, dp5, dp6, dp7, dp8, z1, z2;
//
	nx = 2001;
	k = 1000;
	k1 = k+1;
	dx = 2.0/(nx-1);
	x  = new double [nx];
	qk1 = new double [(k+1)*nx];
	qk2 = new double [(k+1)*nx];
	qk3 = new double [(k+1)*nx];
//
	x[0] = -1.0;
	for (ix = 1; ix < nx; ix++)
		x[ix] = x[ix-1] + dx;
	x[nx-1] = 1.0;
//
	double start = (double)clock() /(double) CLOCKS_PER_SEC;
/*
	polqkm(1, k, x, nx, qk1); // m = 1
	polqkm(2, k, x, nx, qk2); // m = 2
	polqkm(3, k, x, nx, qk3); // m = 3
	*/
	double end = (double)clock() / (double) CLOCKS_PER_SEC;

	dp1 = 0.0;
	dp2 = 0.0;
	dp3 = 0.0;
	dp4 = 0.0;
	dp5 = 0.0;
	dp6 = 0.0;
	dp7 = 0.0;
	dp8 = 0.0;

	for (ix = 0; ix < nx; ix++){
		x2 = x[ix]*x[ix];
//
//		m = 1
		q11 = sqrt( 0.5*(1.0 - x2) );                          // k = 1
		q21 = 3.0*x[ix]*sqrt( (1.0 - x2)/6.0 );                // k = 2
		q31 = (3.0/4.0)*(5.0*x2 - 1.0)*sqrt( (1.0 - x2)/3.0 ); // k = 3
//
//		m = 2
		q22 = 3.0/(2.0*sqrt(6.0))*(1.0 - x2);		            // k = 2
		q32 = 15.0/sqrt(120.0)*x[ix]*(1.0 - x2);	            // k = 3
		q42 = 15.0/(2.0*sqrt(360.0))*(7.0*x2 - 1.0)*(1.0 - x2); // k = 4
//
//		m = 3
		q33 = 15.0/sqrt(720.0)*(1.0 - x2)*sqrt(1.0 - x2);       // k = 3
		q43 = 105.0/sqrt(5040.0)*(1.0 - x2)*x[ix]*sqrt(1.0 - x2); // k = 4
//
//		printf("%12.3f %12.2e %12.2e \n", qk3[ix*k1+0], qk3[ix*k1+1], qk3[ix*k1+2]);

		dp1 += fabs(q11 - qk1[ix*k1+1]); // m = 1, k = 1
		dp2 += fabs(q21 - qk1[ix*k1+2]); // m = 1, k = 2
		dp3 += fabs(q31 - qk1[ix*k1+3]); // m = 1, k = 3

		dp4 += fabs(q22 - qk2[ix*k1+2]); // m = 2, k = 2
		dp5 += fabs(q32 - qk2[ix*k1+3]); // m = 2, k = 3
		dp6 += fabs(q42 - qk2[ix*k1+4]); // m = 2, k = 4

		dp7 += fabs(q33 - qk3[ix*k1+3]); // m = 2, k = 3
		dp8 += fabs(q43 - qk3[ix*k1+4]); // m = 2, k = 4

/*
		dp2 += fabs(p2 - pk[ix*k1+2]);
		dp3 += fabs(p3 - pk[ix*k1+3]);
		dp6 += fabs(p6 - pk[ix*k1+6]);
		dp7 += fabs(p7 - pk[ix*k1+7]);
		dp10 += fabs(p10 - pk[ix*k1+10]);
		*/
	}
	//printf("%12.2e %12.2e %12.2e %12.2e %12.2e %12.2e %12.2e %12.2e \n", dp1/nx, dp2/nx, dp3/nx, dp4/nx, dp5/nx, dp6/nx, dp7/nx, dp8/nx);
//	printf("%12.2e %12.2e %12.2e %12.2e %12.2e %12.2e %12.2e %12.2e \n", dp1/nx, dp2/nx, dp3/nx, dp4/nx, dp5/nx, dp6/nx, dp5/nx, dp6/nx);
//
	z  = new double [7];
	q256 = new double [(512+1)*7];
	z[0] = -1.00;
	z[1] = -0.50;
	z[2] =  0.00;
	z[3] =  0.25;
	z[4] =  0.50;
	z[5] =  0.75;
	z[6] =  1.00;
	polqkm(256, 512, z, 7, q256); // m = 1
	for (ix = 0; ix < 7; ix++)
		z[ix] = q256[ix*(512+1)+512]; // k = 512 (513 total)

//
	printf("\n");
	printf("cpu time %fs\n", (end - start)/3.0);
	printf("Done! \n");
//	std::cin >> cvar;
//
}