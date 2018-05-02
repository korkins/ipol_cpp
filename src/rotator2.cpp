#include <cmath>    /* abs, sqrt */
#define TINY 1.0e-8 // to compare doubles
//
void rotator2(double const smui, double const cmui,
	          double const smus, double const cmus,
	          double const saz, double const caz,
	          double &s2x, double &c2x) {
/*------------------------------------------------------------------------------
PURPOSE:
	To compute SIN and COS of the doubled rotation angles for the inverse (2nd)
    rotation in VRTE.

INPUT:
	smui   d[1]   sin(pi-incident_zenith); smui not equal 1.0
	cmui   d[1]   cos(pi-incident_zenith); cmui not equal 0.0; down: cmui > 0
	smus   d[1]   sin(pi-scattered_zenith); smus not equal 1.0
	cmus   d[1]   cos(pi-scattered_zenith); cmus not equal 0.0; down: cmus > 0
	saz    d[1]   sin(azimuth); azimuth = 0 - glint, forward scattering 
	caz    d[1]   cos(azimuth)

OUTPUT:
	s2x    d[1]   sin(2*x2)
	c2x    d[1]   cos(2*x2)

COMMENTS:
	THINKME: how to return 2 variabels, s2x and c2x, not via array?

	This is a reduced version of rotator.cpp created by removing computations
	related to the 1st rotation, s1 & c1

REFERENCESS:
	1. rotator.cpp

PROTOTYPE:
	void rotator(double const, double const, double const, double const,
	             double const, double const, double, double);
------------------------------------------------------------------------------*/
//
//	LOCAL VARIABLES
	double
		amui, amus, // absolute values of cosines
		ssca, csca, // sin and cos of the scattering angle
		sx2, cx2;   // same for the 2nd angle
//-----------------------------------------------------------------------------
//
//  By default: no rotation
	s2x = 0.0;
	c2x = 1.0;
//
//  Check for no rotation: BOTH beams travel close to Z or in the principal plane
	amui = abs(cmui);
	amus = abs(cmus);
	if (1.0 - amui < TINY && 1.0 - amus < TINY) return;
	if (abs(saz) < TINY) return;
//
	if (1.0 - amui < TINY) {
		cx2 = 1.0;
		sx2 = 0.0;
	} // if (1.0 - amui < TINY)
	else if (1.0 - amus < TINY) {
		if (1.0 + cmus < TINY) cx2 = -caz;
		if (1.0 - cmus < TINY) cx2 = caz;
		sx2 = -saz;
	} // else if (1.0 - amus < TINY)
	else {
//		General case: off the prinicpal plane, not along the Z-axis
		csca = cmui*cmus + smui*smus*caz;
		ssca = sqrt(1.0 - csca*csca);       // *** SLOW & UNSAFE OPERATION ***
		cx2 = (cmus*csca - cmui)/smus/ssca;
		sx2 = -smui*saz/ssca;
	} // if (1.0 - amui < TINY)
	s2x = 2.0*sx2*cx2;
	c2x = 2.0*cx2*cx2 - 1.0;
//
} // void rotator
/*------------------------------------------------------------------------------
 18/04/28 - rot[2] -> &s2x, &c2x. Tested in sglscatk.
 18/01/20 - First created and tested vs rotator.cpp for
			sza[] = { 0.0, 45.0, 60.0, 100.0, 150.0, 170.0, 180.0 },
			vza[] = { 0.0, 45.0, 60.0, 80.0, 110.0, 150.0, 170.0, 180.0 },
			aza[] = { 0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0, 360.0 };
			8 digits printed out - all agree: rotator - rotator2 = 0.00000000.
*/