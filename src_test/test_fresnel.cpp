#include <stdio.h> /* printf */
#include <math.h>  /* sin, cos */
#define PI 3.141592653589793
#define RAD PI/180.0

void fresnel_ip(double const, double const,
	            double const, double const,
			    double const, double const,
	            double *, double *, int const,
	            double *, double *, double *,
				double *, double *, double *,
				double *, double *, double *);
void print2d(double **, int const, int const, char *);

void main()
{
	int const nsz = 4, nvz = 5, naz = 11, nlin = nsz * nvz * naz, ncol = 12; // sza, vza, aza, fr11 ... fr33
	int ilin, icol, ivz, isz, iaz;
	double refre = 1.33, refim = 0.0;
	double cmui, cmur, smui, smur;
	double sza[] = { 0.0, 45.0,  60.0, 75.0 },
		   vza[] = { 0.0, 30.0,  45.0, 60.0,  80.0 },
		   aza[] = { 0.0, 45.0, -45.0, 90.0, 135.0, 180.0, 225.0, -225.0, 270.0, 315.0, 360.0 };
	double **aout;
	double *fr11, *fr12, *fr13, *fr21, *fr22, *fr23, *fr31, *fr32, *fr33, *caz, *saz;
//
	fr11 = new double [naz];
	fr12 = new double [naz];
	fr13 = new double [naz];
	fr21 = new double [naz];
	fr22 = new double [naz];
	fr23 = new double [naz];
	fr31 = new double [naz];
	fr32 = new double [naz];
	fr33 = new double [naz];
	caz = new double [naz];
	saz = new double [naz];
//
	aout = new double *[nlin];
	for (ilin = 0; ilin < nlin; ilin++)
		aout[ilin] = new double [ncol];
//
	for (ilin = 0; ilin < nlin; ilin++)
		for (icol = 0; icol < ncol; icol++)
			aout[ilin][icol] = 0.0;
//
	for (iaz = 0; iaz < naz; iaz++) {
		caz[iaz] = cos(aza[iaz] * RAD);
		saz[iaz] = sin(aza[iaz] * RAD);
	}
//
	ilin = -1;
	for (isz = 0; isz < nsz; isz++) {
		smui = sin(sza[isz] * RAD);
		cmui = cos(sza[isz] * RAD);
		for (ivz = 0; ivz < nvz; ivz++) {
			smur = sin(vza[ivz] * RAD);
			cmur = cos(vza[ivz] * RAD);
//
			fresnel_ip(refre, refim, smui, cmui, smur, -cmur, saz, caz, naz,
				        fr11, fr12, fr13,
				        fr21, fr22, fr23,
				        fr31, fr32, fr33);
//
			for (iaz = 0; iaz < naz; iaz++) {
				ilin += 1;
				aout[ilin][0] = sza[isz];
				aout[ilin][1] = vza[ivz];
				aout[ilin][2] = aza[iaz];
				aout[ilin][3] = fr11[iaz];
				aout[ilin][4] = fr21[iaz];
				aout[ilin][5] = fr31[iaz];
				aout[ilin][6] = fr12[iaz];
				aout[ilin][7] = fr22[iaz];
				aout[ilin][8] = fr32[iaz];
				aout[ilin][9] = fr13[iaz];
				aout[ilin][10] = fr23[iaz];
				aout[ilin][11] = fr33[iaz];
			}
		}
	}
	print2d(aout, nlin, ncol, "fresnel.txt");
//
	printf("Done! \n");
//	std::cin >> cvar;
//
}