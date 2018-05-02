#include <stdio.h> /* printf */
#include <math.h>  /* sin, cos */
#define PI 3.141592653589793
#define RAD PI/180.0

void roughsrf1(double const, double const, double const, double const, double const,
	           double *, int const, double *);
void print2d(double **, int const, int const, char *);

void main()
{
	int const nsz = 4, nvz = 5, naz = 11, nlin = nsz * nvz*naz, ncol = 4; // sza, vza, aza, rtls
	int ilin, icol, ivz, isz, iaz;
	double wsp = 10.0;
	double cmui, cmur, smui, smur;
	double sza[] = { 0.0, 45.0,  60.0, 75.0 },
		   vza[] = { 0.0, 30.0,  45.0, 60.0,  80.0 },
		   aza[] = { 0.0, 45.0, -45.0, 90.0, 135.0, 180.0, 225.0, -225.0, 270.0, 315.0, 360.0 };
	double **rtls;
	double *r, *caz;
//
	r = new double [naz];
	caz = new double [naz];
//
	rtls = new double *[nlin];
	for (ilin = 0; ilin < nlin; ilin++)
		rtls[ilin] = new double [ncol];
//
	for (ilin = 0; ilin < nlin; ilin++)
		for (icol = 0; icol < ncol; icol++)
			rtls[ilin][icol] = 0.0;
//
	for (iaz = 0; iaz < naz; iaz++) caz[iaz] = cos(aza[iaz] * RAD);
//
	ilin = -1;
	for (isz = 0; isz < nsz; isz++) {
		smui = sin(sza[isz] * RAD);
		cmui = cos(sza[isz] * RAD);
		for (ivz = 0; ivz < nvz; ivz++) {
			smur = sin(vza[ivz] * RAD);
			cmur = cos(vza[ivz] * RAD);
//
			roughsrf1(wsp, smui, cmui, smur, -cmur, caz, naz, r);
//
			for (iaz = 0; iaz < naz; iaz++) {
				ilin += 1;
				rtls[ilin][0] = sza[isz];
				rtls[ilin][1] = vza[ivz];
				rtls[ilin][2] = aza[iaz];
				rtls[ilin][3] = r[iaz];
			}
		}
	}
	print2d(rtls, nlin, ncol, "waves.txt");
//
	printf("Done! \n");
//	std::cin >> cvar;
//
}