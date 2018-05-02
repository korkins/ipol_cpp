#include <cmath> /* sqrt */
//
void rotator2(double const, double const,
			  double const, double const,
	          double const, double const,
	          double *);
//
void fresnel0_ip(double const refre, double const refim,
	             double const smui, double const cmui,
	             double const smur, double const cmur,
			     double *saz, double *caz, int const naz,
	             double *fr1, double *fr2, double *fr3) {
/*------------------------------------------------------------------------------
PURPOSE:
	To compute the Fresnel reflection, with rotation, of the unpolarized light.

INPUT:
	refre   d[1]     refractive index: real part
	refim   d[1]     refractive index: imaginary part
	smui    d[1]     sin of incidence zenith, not equal 1.0
	cmui    d[1]     cos of incidence zenith, cmui > 0.0 (down) - as in RTE
	smur    d[1]     sin of reflection zenith, not equal 1.0
	cmur    d[1]     cos of reflection zenith, cmur < 0.0 (up) - as in RTE
	saz     d[naz]   sin(raz_rt)
	caz     d[naz]   cos(raz_rt); caz_rt = -caz_satellites (MODIS, POLDER)
	naz     i[1]     number of azimuths

OUTPUT:
	fr1, fr2, fr3   d[naz]   rotated 1st column of the Fresnel matrix

COMMENTS:
	This subroutine is a simplified version of [1], built on the basis of [2]
	but with: sin of sza, vza, and az on input. Refer to [3] for definition of
    the coordinate system.

REFERENCESS:
	1. fresnel_ip.cpp
	2. SORD\SRC\FRESNELR0_IP.f90
	3. Hovenier et al: green book

PROTOTYPE:
	void fresnel0_ip(double const, double const,
	                 double const, double const,
	                 double const,  double const,
				     double *, double *, int const,
				     double *,  double *,  double *);
------------------------------------------------------------------------------*/
//
//	LOCAL VARIABLES
	int iaz;
	double ca, cb, cinc, cinc2, cir, c2inc, refr2, rl2, rr2, sir, fr;
//
//	LOCAL ARRAYS
	double *rot;
//-----------------------------------------------------------------------------
	rot = new double[2];
//
	refr2 = refre * refre + refim * refim;
	cir = -cmui * cmur; // -cmur > 0
	sir = smui * smur;
	for (iaz = 0; iaz < naz; iaz++) {
//      cos(Scattering angle {between v and -v'}) = cos(2*Angle_of_incidence)
		c2inc = cir - sir * caz[iaz];
		cinc2 = 0.5 + 0.5 * c2inc;
		cinc = sqrt(cinc2);
		ca = sqrt(refr2 - 1.0 + cinc2);
		cb = refr2 * cinc;
		rl2 = (ca - cb) / (ca + cb);
		rr2 = (cinc - ca) / (cinc + ca);
		rl2 = rl2 * rl2;
		rr2 = rr2 * rr2;
		rotator2(smui, cmui, smur, cmur, saz[iaz], caz[iaz], rot);
		fr1[iaz] = 0.5 * (rl2 + rr2);
		fr = 0.5 * (rl2 - rr2);
		fr2[iaz] = fr * rot[1]; // cos(2*x)
		fr3[iaz] = fr * rot[0]; // sin(2*x)
	} // for iaz
//
	delete[] rot;
} // double fresnel0_ip
/*------------------------------------------------------------------------------
 18/03/17 - delete -> delete[]: http://www.cplusplus.com/doc/tutorial/dynamic/
 18/01/28 - First created, tested vs FRESNEL0.f90 for 1.33+i0.0
			sza[] = { 0.0, 45.0, 60.0, 75.0},
			vza[] = { 0.0, 30.0, 45.0, 60.0, 80.0},
			aza[] = { 0.0,   45.0, -45.0,  90.0, 135.0, 180.0,
			        225.0, -225.0, 270.0, 315.0, 360.0 };
			8 digits agreement.
			*** NOTE: for az < 0, 360 + az was used in the F90 subroutine ***
*/