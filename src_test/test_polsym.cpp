#include <stdio.h>				/* printf */
#include <time.h>
#include <iostream>
#include <math.h>

void polleg(int const, double const *, int const, double *);
void polqkm(int const, int const, double const *, int const, double *);
void polrtm(int const, int const, double const *, int const, double *, double *);

void main()
{
	double *px, *mx, *pk_plus, *pk_minus, *r_plus, *r_minus, *t_plus, *t_minus;
	int ix, ik, nx, nk, im, m;
	double dx, p_plus, p_minus, dp, csym, csymR, csymT,
		   sum_dPk, sum_dQkm, sum_dRkm, sum_dTkm;
//
	double start = (double)clock() / (double)CLOCKS_PER_SEC;
//
	nk = 21;
	nx = 5;
	dx = 0.2;
	pk_plus = new double[nk*nx];
	pk_minus = new double[nk*nx];
	r_plus = new double[nk*nx];
	r_minus = new double[nk*nx];
	t_plus = new double[nk*nx];
	t_minus = new double[nk*nx];
//
	px = new double[nx];
	mx = new double[nx];

//
	px[0] =  1.0;
	mx[0] = -1.0;
	printf("ix = %i   +x =%5.2f   -x = %5.2f \n", 1, px[0], mx[0]);
	for (ix = 1; ix < nx; ix++) {
		px[ix] = px[ix - 1] - dx;
		mx[ix] = mx[ix - 1] + dx;
		printf("ix = %i   +x =%5.2f   -x = %5.2f \n", ix+1, px[ix], mx[ix]);
	}
	printf("\n \n");
//------------------------------------------------------------------------------
//
	polleg(nk-1, px, nx, pk_plus);
	polleg(nk-1, mx, nx, pk_minus);
	sum_dPk = 0.0;
	for (ix = 0; ix < nx; ix++) {
		for (ik = 0; ik < nk; ik++) {
			p_plus = pk_plus[ix*nk+ik];
			p_minus = pk_minus[ix*nk+ik];
			csym = pow(-1.0, ik);
			sum_dPk += fabs(p_plus - csym * p_minus);
		}
	}
	printf("Pk(x),  k = 0:%i,  sum_dp = %10.1e \n", nk - 1, sum_dPk);
//------------------------------------------------------------------------------
//
	sum_dQkm = 0.0;
	m = 19;
	for (im = 1; im < m; im++) {
		polqkm(im, nk - 1, px, nx, pk_plus);
		polqkm(im, nk - 1, mx, nx, pk_minus);
		for (ix = 0; ix < nx; ix++) {
			for (ik = 0; ik < nk; ik++) {
				p_plus = pk_plus[ix*nk + ik];
				p_minus = pk_minus[ix*nk + ik];
				csym = pow(-1.0, ik + im);
				sum_dQkm += fabs(p_plus - csym * p_minus);
			}
		}
	}
	printf("Qkm(x),  k = 0:%i,  m = 1:%i,  sum_dq = %10.1e \n", nk - 1, m, sum_dQkm);
//------------------------------------------------------------------------------
//
	sum_dRkm = 0.0;
	sum_dTkm = 0.0;
	m = 19;
	for (im = 1; im < m; im++) {
		polrtm(im, nk - 1, px, nx, r_plus, t_plus);
		polrtm(im, nk - 1, mx, nx, r_minus, t_minus);
		for (ix = 0; ix < nx; ix++) {
			for (ik = 0; ik < nk; ik++) {
//				Rkm(x)
				p_plus = r_plus[ix*nk + ik];
				p_minus = r_minus[ix*nk + ik];
				csymR = pow(-1.0, ik + im);
				sum_dRkm += fabs(p_plus - csymR * p_minus);
//				Tkm(x)
				p_plus = t_plus[ix*nk + ik];
				p_minus = t_minus[ix*nk + ik];
				csymT = -pow(-1.0, ik + im);
				sum_dTkm += fabs(p_plus - csymT * p_minus);
			}
		}
	}
	printf("Rkm(x),  k = 0:%i,  m = 1:%i,  sum_dr = %10.1e \n", nk - 1, m, sum_dRkm);
	printf("Tkm(x),  k = 0:%i,  m = 1:%i,  sum_dt = %10.1e \n", nk - 1, m, sum_dTkm);
//------------------------------------------------------------------------------
//
	double end = (double)clock() / (double) CLOCKS_PER_SEC;
//
	printf("\n");
	printf("cpu time %fs \n", end - start);
	printf("Done! \n");
	system("pause");
}