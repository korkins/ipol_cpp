#include <stdio.h> /* printf */
#include <math.h>  /* sin, cos */
#define PI 3.141592653589793
#define RAD PI/180.0

void fresnel0_ip(double const, double const,
	          double const, double const,
			  double const, double const,
	          double *, double *, int const,
	          double *, double *, double *);
void print2d(double **, int const, int const, char *);
void rotator2(double const, double const,
	          double const, double const,
	          double const, double const,
	          double *);

void main()
{
	int const nsz = 4, nvz = 5, naz = 11, nlin = nsz * nvz * naz, ncol = 6; // sza, vza, aza, aout123
	int ilin, icol, ivz, isz, iaz;
	double refre = 1.33, refim = 0.0;
	double cmui, cmur, smui, smur;
	double sza[] = { 0.0, 45.0,  60.0, 75.0 },
		   vza[] = { 0.0, 30.0,  45.0, 60.0,  80.0 },
		   aza[] = { 0.0, 45.0, -45.0, 90.0, 135.0, 180.0, 225.0, -225.0, 270.0, 315.0, 360.0 };
	double **aout;
	double *fr1, *fr2, *fr3, *caz, *saz, *s2c2;
//
	fr1 = new double [naz];
	fr2 = new double [naz];
	fr3 = new double [naz];
	caz = new double [naz];
	saz = new double [naz];
	s2c2 = new double[2];
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
			fresnel0_ip(refre, refim, smui, cmui, smur, -cmur, saz, caz, naz, fr1, fr2, fr3);
//
			for (iaz = 0; iaz < naz; iaz++) {
				rotator2(smui, cmui, smur, -cmur, saz[iaz], caz[iaz], s2c2);
				ilin += 1;
				aout[ilin][0] = sza[isz];
				aout[ilin][1] = vza[ivz];
				aout[ilin][2] = aza[iaz];
				aout[ilin][3] = fr1[iaz];
				aout[ilin][4] = fr2[iaz];
				aout[ilin][5] = fr3[iaz];
			}
		}
	}
	print2d(aout, nlin, ncol, "fresnel0.txt");
//
	printf("Done! \n");
//	std::cin >> cvar;
//
}