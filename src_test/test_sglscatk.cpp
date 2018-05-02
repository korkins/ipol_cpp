#include <stdio.h> /* printf */
#include <math.h>  /* sin, cos */
#include <time.h>
#include <iostream>
#define PI 3.141592653589793
#define RAD PI/180.0
//
void print2d(double **, int const, int const, char *);
void read2d(char *, int const, int const, double **);
void sglscatk(double const *, double const *, int const, int const, double const *, double const *,
	int const, double const *, double const *, int const, double const, double const *,
	int const, double const *, double const *, int const, double *, double *, double *);
//
void main()
{
	int ilin, iaz, imu0, imu, ilr, ik, isca, naz, nmu0, nmu, nup, nk, nlr, nlin, ncol, nsca;
	double tau0, di1, dq1, du1, sum_i1s, sum_q1s, sum_u1s, df, d, ssa05;

	double r1k[3], r2k[3];

	double *mu, *smu, *mu0, *smu0, *aza, *caz, *saz, *tau, *ssa, *a1k, *b1k, **xk_txt,
		*i1, *q1, *u1, *i14, *q14, *u14, **tau_ssa, **sord_test48, *theta;

//  Case 1: Hovenier et al, 1 layer ( +extra solar-view geometries )
//
//  Size of arrays
	naz = 3;
	nmu0 = 5;
	nmu = 7;
	nup = 3;
	nk = 117;
	nlr = 1;
	nsca = naz * nmu0 * nmu;
//
//  Allocate arrays
	aza = new double[naz]; caz = new double[naz]; saz = new double[naz];
	mu0 = new double[nmu0]; smu0 = new double[nmu0];
	mu = new double[nmu]; smu = new double[nmu];
	tau = new double[nlr]; ssa = new double[nlr];
	a1k = new double[nlr * nk]; b1k = new double[nlr * nk];
	i1 = new double[nsca]; q1 = new double[nsca]; u1 = new double[nsca];
//  moments from txt file
	nlin = nk;
	ncol = 7; // k, a1k, a2k, a3k, a4k, b1k, b2k
	xk_txt = new double *[nlin];
	for (ilin = 0; ilin < nlin; ilin++)
		xk_txt[ilin] = new double[ncol];
//
//  Input
	aza[0] = 0.0; aza[1] = 90.0; aza[2] = 180.0; // all from benchmark, bm
	mu0[0] = 1.0; mu0[1] = 0.75; mu0[2] = 0.6; mu0[3] = 0.4; mu0[4] = 0.2; // imu0_bm = 2
	mu[0] = -1.0; mu[1] = -0.6; mu[2] = -0.2; mu[3] = 0.1; mu[4] = 0.6; mu[5] = 0.8; mu[6] = 1.0; // imu_b = 1, 4
	tau[0] = 0.2; ssa[0] = 1.0;
//
	for (iaz = 0; iaz < naz; iaz++) {
		caz[iaz] = cos(aza[iaz] * RAD);
		saz[iaz] = sin(aza[iaz] * RAD);
	}
	for (imu0 = 0; imu0 < nmu0; imu0++) smu0[imu0] = sqrt(1.0 - mu0[imu0] * mu0[imu0]);
	for (imu = 0; imu < nmu; imu++) smu[imu] = sqrt(1.0 - mu[imu] * mu[imu]);
	tau0 = 0.0;
	for (ilr = 0; ilr < nlr; ilr++) tau0 += tau[ilr];
//
	read2d("Xk0117_0804.txt", nlin, ncol, xk_txt);
//	print2d(xk_txt, nlin, ncol, "xkout.txt");
	for (ilr = 0; ilr < nlr; ilr++)
		for (ik = 0; ik < nk; ik++) {
			a1k[ilr*nk + ik] = 0.5*ssa[ilr] * xk_txt[ik][1]; // a1k
			b1k[ilr*nk + ik] = 0.5*ssa[ilr] * xk_txt[ik][5]; // b1k
		}
//
	sglscatk(mu, smu, nmu, nup, mu0, smu0, nmu0, caz, saz, naz, tau0, tau, nlr, a1k, b1k, nk,
		i1, q1, u1); // scale_factor = 2pi;
//
	printf("TEST 1: 1 layer \n");
	for (iaz = 0; iaz < naz; iaz++)
		for (imu0 = 0; imu0 < nmu0; imu0++)
			for (imu = 0; imu < nmu; imu++) {
				isca = iaz * nmu0 * nmu + imu0 * nmu + imu;
				printf("%8.2f %6.2f %6.2f %12.6e % 12.6e %12.6e \n", aza[iaz], mu0[imu0], mu[imu],
					0.5*i1[isca], 0.5*q1[isca], 0.5*u1[isca]); // 0.5 to compare with Hovenier
			}
	printf("\n"); printf("\n"); printf("\n");
	//-----------------------------------------------------------------------------------------------------
	//
	// TEST 2:
	delete[] tau; delete[] ssa; delete[] a1k; delete[] b1k;
	nlr = 4;
	tau = new double[nlr]; ssa = new double[nlr];
	a1k = new double[nlr * nk]; b1k = new double[nlr * nk];

	tau[0] = 0.05; tau[1] = 0.03; tau[2] = 0.02; tau[3] = 0.1;
	for (ilr = 0; ilr < nlr; ilr++) {
		ssa[ilr] = 1.0;
		for (ik = 0; ik < nk; ik++) {
			a1k[ilr*nk + ik] = 0.5*ssa[ilr] * xk_txt[ik][1]; // a1k
			b1k[ilr*nk + ik] = 0.5*ssa[ilr] * xk_txt[ik][5]; // b1k
		}
	}
//
//  Reset output
	i14 = new double[nsca]; q14 = new double[nsca]; u14 = new double[nsca];
//
	sglscatk(mu, smu, nmu, nup, mu0, smu0, nmu0, caz, saz, naz, tau0, tau, nlr, a1k, b1k, nk,
		i14, q14, u14); // scale_factor = 2pi;
//
	printf("TEST 2: 4 layers MINUS 1 layer -> xi-deviation \n");
	sum_i1s = 0.0; sum_q1s = 0.0; sum_u1s = 0.0;
	for (iaz = 0; iaz < naz; iaz++)
		for (imu0 = 0; imu0 < nmu0; imu0++)
			for (imu = 0; imu < nmu; imu++) {
				isca = iaz * nmu0 * nmu + imu0 * nmu + imu;
				di1 = i1[isca] - i14[isca];
				dq1 = q1[isca] - q14[isca];
				du1 = u1[isca] - u14[isca];
				sum_i1s += di1 * di1;
				sum_q1s += dq1 * dq1;
				sum_u1s += du1 * du1;
				// printf("%8.2f %6.2f %6.2f %12.6e % 12.6e %12.6e \n", aza[iaz], mu0[imu0], mu[imu],
				// i1[isca] - i14[isca], q1[isca] - q14[isca], u1[isca] - u14[isca]);
				// 0.5*i14[isca], 0.5*q14[isca], 0.5*u14[isca]);
			}
	printf("%12.6e % 12.6e %12.6e \n", sqrt(sum_i1s)/nsca, sqrt(sum_q1s) / nsca, sqrt(sum_u1s) / nsca);
	printf("\n"); printf("\n"); printf("\n");

	delete[] tau; delete[] ssa; delete[] a1k; delete[] b1k;
	delete[] aza; delete[] caz; delete[] saz;
	delete[] mu0; delete[] smu0;
	delete[] mu; delete[] smu;
	delete[] i1; delete[] q1; delete[] u1;
	i1 = new double[nsca]; q1 = new double[nsca]; u1 = new double[nsca];
//-----------------------------------------------------------------------------------------------------
//
//  TEST 3: SORD_T048 - 30 optical layers, Rayleigh
//
//  Size of arrays
	naz = 37;
	nmu0 = 1;
	nmu = 36;
	nup = 18;
	nk = 3;
	nlr = 30;
	nsca = naz * nmu0 * nmu;
//
//  Allocate arrays
	aza = new double[naz]; caz = new double[naz]; saz = new double[naz];
	mu0 = new double[nmu0]; smu0 = new double[nmu0];
	mu = new double[nmu]; smu = new double[nmu]; theta = new double[nmu];
	tau = new double[nlr]; ssa = new double[nlr];
	a1k = new double[nlr * nk]; b1k = new double[nlr * nk];
	i1 = new double[nsca]; q1 = new double[nsca]; u1 = new double[nsca];
//
//	Rayleigh moments
	d = 0.03;
	df = (1.0 - d) / (1.0 + 0.5 * d);
	r1k[0] = 1.0; r1k[1] = 0.0; r1k[2] = 0.5*df;
	r2k[0] = 0.0; r2k[1] = 0.0; r2k[2] = sqrt(1.5)*df;
//
//	Read tau & ssa
	nlin = nlr;
	ncol = 2;
	tau_ssa = new double *[nlin];
	for (ilin = 0; ilin < nlin; ilin++)
		tau_ssa[ilin] = new double[ncol];
	read2d("tau_ssa_test48.txt", nlin, ncol, tau_ssa);
	for (ilr = 0; ilr < nlr; ilr++) {
		tau[ilr] = tau_ssa[ilr][0];
		ssa05 = 0.5*tau_ssa[ilr][1];
		for (ik = 0; ik < nk; ik++) {
			a1k[ilr*nk + ik] = ssa05 * r1k[ik]; // a1k
			b1k[ilr*nk + ik] = ssa05 * r2k[ik]; // b1k
		}
	}
//
//  Benchmark
	nlin = nsca;
	ncol = 8;
	sord_test48 = new double *[nlin];
	for (ilin = 0; ilin < nlin; ilin++)
		sord_test48[ilin] = new double[ncol];
	read2d("sord_test48.txt", nlin, ncol, sord_test48);
//
//  Solar-view geometry
	mu0[0] = 0.5;
	smu0[0] = sqrt(1.0 - mu0[0] * mu0[0]);
	for (imu = 0; imu < nup; imu++) {
		theta[imu] = imu * 5.0;
		theta[imu + nup] = (imu + 1 + nup) * 5.0;
		mu[imu] = -cos(RAD * imu * 5.0);
		mu[imu + nup] = -cos(RAD * (imu + 1 + nup) * 5.0);
	}
	for (imu = 0; imu < nmu; imu++) smu[imu] = sqrt(1.0 - mu[imu] * mu[imu]);
	for (iaz = 0; iaz < naz; iaz++) {
		aza[iaz] = iaz * 5.0;
		caz[iaz] = cos(RAD * aza[iaz]);
		saz[iaz] = sin(RAD * aza[iaz]);
	}
//
	tau0 = 0.0;
	for (ilr = 0; ilr < nlr; ilr++) tau0 += tau[ilr];
//  
	double start = (double)clock() /(double) CLOCKS_PER_SEC;
	for (int irun = 0; irun < 1000; irun++)
	sglscatk(mu, smu, nmu, nup, mu0, smu0, nmu0, caz, saz, naz, tau0, tau, nlr, a1k, b1k, nk,
		i1, q1, u1); // scale_factor = 2pi;
	double end = (double)clock() / (double)CLOCKS_PER_SEC;
//
	printf("TEST 3: IPRT Rayleigh 30 layers -> xi-deviation \n");
	sum_i1s = 0.0; sum_q1s = 0.0; sum_u1s = 0.0;
	for (iaz = 0; iaz < naz; iaz++)
		for (imu0 = 0; imu0 < nmu0; imu0++)
			for (imu = 0; imu < nmu; imu++) {
				isca = iaz * nmu0 * nmu + imu0 * nmu + imu;
				di1 = i1[isca] / (2.0 * PI) - sord_test48[isca][5];
				dq1 = q1[isca] / (2.0 * PI) - sord_test48[isca][6];
				du1 = u1[isca] / (2.0 * PI) - sord_test48[isca][7];
				sum_i1s += di1 * di1;
				sum_q1s += dq1 * dq1;
				sum_u1s += du1 * du1;
			}
	printf("%12.6e % 12.6e %12.6e \n", sqrt(sum_i1s) / nsca, sqrt(sum_q1s) / nsca, sqrt(sum_u1s) / nsca);
	printf("\n"); printf("\n"); printf("\n");
	printf("cpu time %fs\n", end - start);
	printf("Done!\n");
	system("pause");
//
}