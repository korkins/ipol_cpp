#include <stdio.h>				/* printf */
#include <time.h>
#include <iostream>
#include <math.h>
//using namespace std;

void print1d(double *, int const, char *);
void print2d(double **, int const, int const, char *);

void main()
{
	int nm, nn, nk, ix, iy;
	double **a, **b, *v;
//
	double start = (double)clock() / (double)CLOCKS_PER_SEC;
//
	nm = 2;
	nn = 3;
	nk = 4;
//
//  1d in
	v = new double [nk];
	for (ix = 0; ix < nk; ix++)
		v[ix] = 1.0*ix;
	print1d(v, nk, "v.txt");
//
//  2d a[nm][nn]
	a = new double *[nm];        // rows first
	for (iy = 0; iy < nm; iy++)
		a[iy] = new double [nn]; // columns next
	for (iy = 0; iy < nm; iy++)
		for (ix = 0; ix < nn; ix++)
			a[iy][ix] = 1.0*ix + iy;
	print2d(a, nm, nn, "a.txt");
//
//  2d b[nn][nk]
	b = new double *[nn];
	for (iy = 0; iy < nn; iy++)
		b[iy] = new double[nk];
	for (iy = 0; iy < nn; iy++)
		for (ix = 0; ix < nk; ix++)
			b[iy][ix] = 1.0*ix*iy + ix + iy;
	print2d(b, nn, nk, "b.txt");
//
	double end = (double)clock() / (double) CLOCKS_PER_SEC;
//
	delete(a);
	delete(b);
	delete(v);
//
	printf("\n");
	printf("cpu time %fs\n", end - start);
	printf("Done! \n");
//	std::cin >> cvar;
//
}