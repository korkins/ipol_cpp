#include <cmath> /* sqrt */
//
void rotator(double const, double const,
			  double const, double const,
	          double const, double const,
	          double *);
//
void fresnel_ip(double const refre, double const refim,
	            double const smui, double const cmui,
	            double const smur, double const cmur,
			    double *saz, double *caz, int const naz,
	            double *fr11, double *fr12, double *fr13,
	            double *fr21, double *fr22, double *fr23,
				double *fr31, double *fr32, double *fr33) {
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
	fr11...fr33   d[naz]   rotated Fresnel matrix

COMMENTS:
	This subroutine is built on the basis of [1] but with: sin of sza, vza, and
	az on input. Refer to [3] for definition of the coordinate system.

REFERENCESS:
	1. SORD\SRC\FRESNELR0_IP.f90
	2. Hovenier et al: green book

PROTOTYPE:
	void fresnel_ip(double const, double const,
	                double const, double const,
	                double const,  double const,
				    double *, double *, int const,
				    double *,  double *,  double *,
					double *,  double *,  double *,
					double *,  double *,  double *);
------------------------------------------------------------------------------*/
//
//	LOCAL VARIABLES
	int iaz;
	double ca, cb, cinc, cinc2, cir, c2inc, f1, f2, f3, refr2, rl, rl2, rr, rr2,
		   sir;
//
//	LOCAL ARRAYS
	double *rot;
//-----------------------------------------------------------------------------
	rot = new double[4];
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
		rl = (cb - ca) / (cb + ca);
		rr = (cinc - ca) / (cinc + ca);
		rl2 = rl * rl;
		rr2 = rr * rr;
		rotator(smui, cmui, smur, cmur, saz[iaz], caz[iaz], rot); // s1 c1 s2 c2
		f1 = 0.5 * (rl2 + rr2);
		f2 = 0.5 * (rl2 - rr2);
		f3 = rl * rr;
//
		fr11[iaz] = f1;
		fr21[iaz] = f2 * rot[3]; // c2
		fr31[iaz] = f2 * rot[2]; // s2
//
		fr12[iaz] = f2 * rot[1]; // c1
		fr22[iaz] = f1 * rot[1] * rot[3] - f3 * rot[0] * rot[2]; // c1*c2 - s1*s2
		fr32[iaz] = f1 * rot[1] * rot[2] + f3 * rot[0] * rot[3]; // c1*s2 + s1*c2
//
		fr13[iaz] = -f2 * rot[0]; // -s1
		fr23[iaz] = -f1 * rot[0] * rot[3] - f3 * rot[1] * rot[2]; // -s1*c2 - c1*s2
		fr33[iaz] = -f1 * rot[0] * rot[2] + f3 * rot[1] * rot[3]; // -s1*s2 + c1*c2
	} // for iaz
//
	delete[] rot;
} // double fresnel_ip
/*------------------------------------------------------------------------------
 18/03/17 - delete -> delete[]: http://www.cplusplus.com/doc/tutorial/dynamic/
 18/02/17 - First created, tested vs FRESNEL0.f90 for 1.33+i0.0
			sza[] = { 0.0, 45.0, 60.0, 75.0},
			vza[] = { 0.0, 30.0, 45.0, 60.0, 80.0},
			aza[] = { 0.0,   45.0, -45.0,  90.0, 135.0, 180.0,
			        225.0, -225.0, 270.0, 315.0, 360.0 };
			8 digits agreement.
			*** NOTE: for az < 0, 360 + az was used in the F90 subroutine ***
*/