#include <stdio.h> /* printf */
#include <math.h>  /* sin, cos */
#define PI 3.141592653589793
#define RAD PI/180.0

void rotator(double const, double const, double const, double const,
	         double const, double const, double *);
void rotator2(double const, double const, double const, double const,
			  double const, double const, double *);
void print2d(double **, int const, int const, char *);

void main()
{
	int const nsz = 7, nvz = 8, naz = 9, nlin = nsz * nvz*naz, ncol = 7; // sza, vza, aza, s1, c1, s2, c2
	int ilin, icol, ivz, isz, iaz;
	double cmui, cmus, smui, smus, saz, caz, s1, c1, s2, c2;
	double sza[] = { 0.0, 45.0, 60.0, 100.0, 150.0, 170.0, 180.0 },
		vza[] = { 0.0, 45.0, 60.0, 80.0, 110.0, 150.0, 170.0, 180.0 },
		aza[] = { 0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0, 360.0 };
	double **rot;
	double *s1c1s2c2, *s2c2;
//
	s1c1s2c2 = new double[4];
	s2c2 = new double[2];
//
	rot = new double *[nlin];
	for (ilin = 0; ilin < nlin; ilin++)
		rot[ilin] = new double [ncol];
//
	for (ilin = 0; ilin < nlin; ilin++)
		for (icol = 0; icol < ncol; icol++)
			rot[ilin][icol] = 0.0;
//
	s1 = 10.0, s2 = 10.0, c1 = 10.0, c2 = 10.0;
	ilin = -1;
	for (isz = 0; isz < nsz; isz++) {
		smui = sin(sza[isz] * RAD);
		cmui = cos(sza[isz] * RAD);
		for (ivz = 0; ivz < nvz; ivz++) {
			smus = sin(vza[ivz] * RAD);
			cmus = cos(vza[ivz] * RAD);
			for (iaz = 0; iaz < naz; iaz++) {
				saz = sin(aza[iaz] * RAD);
				caz = cos(aza[iaz] * RAD);
				rotator(smui, cmui, smus, cmus, saz, caz, s1c1s2c2);
				rotator2(smui, cmui, smus, cmus, saz, caz, s2c2);
				ilin += 1;
				rot[ilin][0] = sza[isz];
				rot[ilin][1] = vza[ivz];
				rot[ilin][2] = aza[iaz];
				rot[ilin][3] = s1c1s2c2[0];
				rot[ilin][4] = s1c1s2c2[1];
				rot[ilin][5] = s1c1s2c2[2]; // -s2c2[0];
				rot[ilin][6] = s1c1s2c2[3]; // -s2c2[1];
			}
		}
	}
	print2d(rot, nlin, ncol, "rotator.txt");
//
	printf("Done! \n");
//	std::cin >> cvar;
//
}