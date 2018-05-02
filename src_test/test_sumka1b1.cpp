#include <stdio.h>				/* printf */
#include <time.h>
#include <iostream>
#include <math.h>
#define PI 3.141592653589793
//using namespace std;

void sumka1b1(double *, int const, double *, double *, int const, double *, double *);
void read2d(char *, int const, int const, double **);
void print2d(double **, int const, int const, char *);
void print1d(double *, int const, char *);

void main()
{
	int const nx = 10;
	double a1[nx], b1[nx];
	double x[] = { -1.0, -0.75, -0.5, -0.25, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };
	int ix, nk, nlin, ncol, ilin, icol, nlin_pm, ncol_pm, nlin_xk, ncol_xk, nmu;
	double r1, r2, **aout, *a1k, *b1k, **ain, **pmat, **xk, *mu, *x1, *x2, *f1, *f2, *x1k, *x2k;
//
	nk = 3;
	a1k  = new double [nk];
	b1k = new double [nk];
	a1k[0] = 1.0;
	a1k[1] = 0.0;
	a1k[2] = 0.5;
	b1k[0] = 0.0;
	b1k[1] = 0.0;
	b1k[2] = 0.5*sqrt(6.0);
//
	nlin = nx;
	ncol = 5;
	aout = new double *[nlin];
	for (ilin = 0; ilin < nlin; ilin++)
		aout[ilin] = new double[ncol];
//
	sumka1b1(x, nx, a1k, b1k, nk, a1, b1);
//
	for (ix = 0; ix < nx; ix++) {
		aout[ix][0] = x[ix];
		aout[ix][1] = a1[ix];
		aout[ix][2] = b1[ix];
		r1 = 0.75*(1.0 + x[ix] * x[ix]);
		r2 = 0.75*(x[ix] * x[ix] - 1.0);
		aout[ix][3] = 100.0*(1.0 - a1[ix] / r1);
		aout[ix][4] = r2 - b1[ix];
	}
	delete a1k;
	delete b1k;
	print2d(aout, nlin, ncol, "a1b1_r.txt");
	delete aout;
//-------------------------------------------------------------------
/*
	nk = 2000;
	a1k = new double[nk];
	b1k = new double[nk];
//
	for (ix = 0; ix < nx; ix++) {
		sumka1b1(x[ix], a1k, b1k, nk, a1b1);
		aout[ix][0] = x[ix];
		aout[ix][1] = a1b1[0];
		aout[ix][2] = a1b1[1];
		aout[ix][3] = -999.0;
		aout[ix][4] = -999.0;
*/

/*
	nlin = 4;
	ncol = 3;
	ain = new double *[nlin];
	aout = new double *[nlin];
	for (ilin = 0; ilin < nlin; ilin++) {
		aout[ilin] = new double[ncol];
		ain[ilin] = new double[ncol];
	}
	read2d("ain.txt", nlin, ncol, ain);
	print2d(ain, nlin, ncol, "aout.txt");
*/
	printf("2nd test >>> \n");
	nmu = 1900;
	nlin_pm = nmu;
	ncol_pm = 5;
	
	pmat = new double *[nlin_pm];
	for (ilin = 0; ilin < nlin_pm; ilin++)
		pmat[ilin] = new double[ncol_pm];
	printf("read pmat >>> \n");
	read2d("pm.txt", nlin_pm, ncol_pm, pmat);
	printf("done read pmat \n");
//
	for (icol = 0; icol < ncol_pm; icol++)
		printf("%12.8f", pmat[nmu-1][icol]);
	printf("\n last line of pmat printed ok \n");
//
	mu = new double[nmu];
	x1 = new double[nmu];
	x2 = new double[nmu];
	f1 = new double[nmu];
	f2 = new double[nmu];
	for (ilin = 0; ilin < nlin_pm; ilin++) {
		mu[ilin] = cos(pmat[ilin][0] * 0.01745329251994329576923690768489);
		x1[ilin] = pmat[ilin][1];
		x2[ilin] = pmat[ilin][2];
	}
//
	printf("read xk \n");
	nk = 932;
	nlin_xk = nk;
	ncol_xk = 7;

	xk = new double *[nlin_xk];
	for (ilin = 0; ilin < nlin_xk; ilin++)
		xk[ilin] = new double[ncol_xk];
	printf("read xk >>> \n");
	read2d("xk.txt", nlin_xk, ncol_xk, xk);
	printf("done read xk \n");
	//
	for (icol = 0; icol < ncol_xk; icol++)
		printf("%12.8f", xk[10][icol]);
	printf("\n 10th line of xk printed ok \n");
//
	x1k = new double[nk];
	x2k = new double[nk];
	for (ilin = 0; ilin < nlin_xk; ilin++) {
		x1k[ilin] = xk[ilin][1];
		x2k[ilin] = xk[ilin][5];
	}
//
	sumka1b1(mu, nmu, x1k, x2k, nk, f1, f2);
//
	print1d(f1, nmu, "f1.txt");
	print1d(f2, nmu, "f2.txt");
//
	printf("Done! \n");
//
}