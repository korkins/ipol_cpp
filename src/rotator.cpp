#include <cmath>    /* abs, sqrt */
#define TINY 1.0e-8 // to compare doubles
//
void rotator(double const smui, double const cmui,
	         double const smus, double const cmus,
	         double const saz, double const caz,
	         double *rot) {
/*------------------------------------------------------------------------------
PURPOSE:
	To compute SIN and COS of the doubled rotation angles for direct and inverse
    rotations in VRTE. Attn: mu = -cos(zenith) > 0 for the downward directions

INPUT:
	smui   d[1]   sin(pi-incident_zenith); smui not equal 1.0
	cmui   d[1]   cos(pi-incident_zenith); cmui not equal 0.0; down: cmui > 0
	smus   d[1]   sin(pi-scattered_zenith); smus not equal 1.0
	cmus   d[1]   cos(pi-scattered_zenith); cmus not equal 0.0; down: cmus > 0
	saz    d[1]   sin(azimuth); azimuth = 0 - glint, forward scattering 
	caz    d[1]   cos(azimuth)

OUTPUT:
	rot   d[4]   [s1 c1 s2 c2] where s & c means sin and cos,
	                           1 & 2-1st and 2nd rotation, e.g. c1 = cos(2*x1)

COMMENTS:
	THINKME: sqrt - to be removed to a calling subroutine, e.g. Fresnel ??

	The function is used for reflection of diffuse light from surface. In this
	case, cmui > 0 & cmus < 0.

    Positive mu are measured opposed to positive Z-direction, i.e. mu > 0(down)
    corresponds to zenith angle Theta = (pi/2:pi].

    The coordinate system is defined in [1, p.66]. 'Theta' means zenith angle.

    As defined in [1, p.11, Eq.(1.51)], the rotation matrix through an angle
    X > 0

        |1   0         0         0|   |1  0  0  0|
        |0   cos(2X)   sin(2X)   0|   |0  C  S  0|
    R = |0  -sin(2X)   cos(2X)   0| = |0 -S  C  0|                          (1)
        |0   0         0         1|   |0  0  0  1|

    In this subroutine, not the matrix R but the elements S, C are computed.
    The '-' in (1) is must be considered in the calling program. S1, C1 and S2,
    C2 are defined for X1 and X2, respectively. Both rotations are considered

    B = R(X2) A R(X1)                                                       (2)

    yet the first one, R(X1), is not necessary, for instance, in case of single
    scattering of natural light (use rotator2.cpp in the case)

    Equations are given in [1, pp.69-70, Eqs.(3.10)-(3.23)]. For the cases when
    aza = 0, pi, 2pi no rotation is performed. Two special cases must be
    considered using L'Hopital rule when [1, Eqs.(3.15)-(3.16)] can not be used
    [2].

    Scattering in the horizontal directions is not considered. Refer to [1, 2]
    for details.

    a) For the case az = 0, pi, 2pi no rotation needed.
    b) For the case when incident AND(!) scattered rays travel along the
    zenith/nadir direction (0, pi) no rotation is needed
    c) For the case when only incident ray travels along the zenith/nadir, the
    SECOND rotator, R2, is not needed (scattering plane contains the nadir/
    zenith direction and the direction of scattering).
    d) For the case when only scattered ray travels along the zenith/nadir
    direction, the FIRST rotator, R1, is not needed (scattering plane contains
    the nadir/zenith and the direction of incidence)

    All the cases (a)-(d) are included in this subroutine.

REFERENCESS:
	1. Hovenier JW, van der Mee C, Domke H. Transfer of polarized light in
       planetary atmospheres. Basic concepts and practical methods. Dordrecht:
       Kluwer Academic Publishers, 2004.
    2. Rozanov VV, Bremen University. Private communication.

PROTOTYPE:
	void rotator(double const, double const, double const, double const,
	             double const, double const, double *);
------------------------------------------------------------------------------*/
//
//	LOCAL VARIABLES
	double
		amui, amus, // absolute values of cosines
		ssca, csca, // sin and cos of the scattering angle
		sx1, cx1,   // sin & cos of the 1st angle of rotation
		sx2, cx2;   // same for the 2nd angle
//-----------------------------------------------------------------------------
//
//  By default: no rotation (beams travel along Z-axis or in the principal plane)
	rot[0] = 0.0, rot[1] = 1.0, rot[2] = 0.0, rot[3] = 1.0; // s1 c1 s2 c2
//
//  Check for no rotation: BOTH beams travel close to Z or in the principal plane
	amui = abs(cmui);
	amus = abs(cmus);
	if (1.0 - amui < TINY && 1.0 - amus < TINY) return;
	if (abs(saz) < TINY) return;
//
	if (1.0 - amui < TINY) {
		if (1.0 + cmui < TINY) cx1 = -caz;  // '-' [2], cos(Theta') = +1
		if (1.0 - cmui < TINY) cx1 = caz;   // '+' [2], cos(Theta') = -1
		sx1 = -saz;                         // [1, p.69, Eq.(3.13)]: -sin(fi-fi')
		cx2 = 1.0;
		sx2 = 0.0;
	} // if (1.0 - amui < TINY)
	else if (1.0 - amus < TINY) {
		cx1 = 1.0;
		sx1 = 0.0;
		if (1.0 + cmus < TINY) cx2 = -caz;  // '-' [2], cos(Theta) = +1
		if (1.0 - cmus < TINY) cx2 = caz;   // '-' [2], cos(Theta) = -1
		sx2 = -saz;                         // [1, p.69, Eq.(3.13)]: -sin(fi-fi')
	} // else if (1.0 - amus < TINY)
	else {
//		General case: off the prinicpal plane, not along the Z-axis
		csca = cmui*cmus + smui*smus*caz;
		ssca = sqrt(1.0 - csca*csca);       // *** SLOW & UNSAFE OPERATION ***
		cx1 = (cmui*csca - cmus)/smui/ssca; // [1, p.70, Eq.(3.20)]
		cx2 = (cmus*csca - cmui)/smus/ssca; // [1, p.70, Eq.(3.21)]
		sx1 = -smus*saz/ssca;               // [1, p.70, Eq.(3.23)]: -sin(fi-fi')
		sx2 = -smui*saz/ssca;               // [1, p.70, Eq.(3.23)]: -sin(fi-fi')
	} // if (1.0 - amui < TINY)
	rot[0] = 2.0*sx1*cx1;       // s1
	rot[1] = 2.0*cx1*cx1 - 1.0; // c1
	rot[2] = 2.0*sx2*cx2;       // s2
	rot[3] = 2.0*cx2*cx2 - 1.0; // c2
//
} // void rotator
/*------------------------------------------------------------------------------
 18/01/20 - Tested vs ROTATOR.F90 for
			sza[] = { 0.0, 45.0, 60.0, 100.0, 150.0, 170.0, 180.0 },
			vza[] = { 0.0, 45.0, 60.0, 80.0, 110.0, 150.0, 170.0, 180.0 },
			aza[] = { 0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0, 360.0 };
			8 digits printed out - all agree: F90 - cpp = 0.00000000.
 18/01/17 - First created from the SORD's ROTATOR.f90; completely rewritten
*/