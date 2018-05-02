#include <cmath> /* exp, sqrt */
#define NADIR 1.0e-8             // near-nadir directions
#define PI    3.1415926535897932 //
#define SQRPI 1.7724538509055160 // sqrt(pi)
#define SQR2  1.4142135623730950 // sqrt(2.0)
#define S1    5.12e-3            // sgm2 = s1*wsp + s2 
#define S2    3.0e-3             // sgm2 = s1*wsp + s2
//
void roughsrf1(double const wsp, double const smui, double const cmui,
		       double const smur, double const cmur, double *caz, int const naz,
	           double *wsrf) {
/*------------------------------------------------------------------------------
PURPOSE:
   To compute the RTLS BRDF for an array of azimuths

INPUT:
	wsp    d[1]     wind speed, m/s
	smui   d[1]     sin of incidence zenith, not equal 1.0
	cmui   d[1]     cos of incidence zenith, cmui > 0.0 (down) - as in RTE
	smur   d[1]     sin of reflection zenith, not equal 1.0
	cmur   d[1]     cos of reflection zenith, cmur < 0.0 (up) - as in RTE
	caz    d[naz]   cos(raz_rt); caz_rt = -caz_satellites (MODIS, POLDER)
	naz    i[1]     number of azimuths

OUTPUT:
	wsrf   d[naz]   waves

COMMENTS:
    This code was build using codes SHARM (Dr.Alexei Lyapustin, NASA GSFC) and
    ocean.phase.f (Dr.Michael Mishchenko, NASA GISS, available online at
    http://www.giss.nasa.gov/staff/mmishchenko/brf/).
 
    The complementary error function, REAL*8 DERFC, is computed using intrinsic
    function. If not available, derfc.f function should be downloaded from
    www.netlib.org with dependencies). ERFC is defined as
 
    erfc(x) = 1 - erf(x) = 2/sqrt(pi)*integral(exp(-t*t)dt, x, +inf)         (1)
 
    The bidirectional shadowing function is
 
    S = 1/(1 + F(mu) + F(mu')),                                              (2)
 
    where
 
    F = 1/2*{ exp(-v*v)/(sqrt(pi)*v) - erfc(v) },                            (3)
 
    erfc(x) is defined in Eq.(1) and v = mu/(S*sqrt(1 - mu2)), S is the slope
    distribution parameter usually defined as 'sigma'. Eq.(3) is taken from
    [1, Eq.(18)] with a corrected misprinting (compare with [2, Eq.(6)] or with
    [3, Eq.(88)]).
 
    NOTE: F(mu) = 0 if no shadowing is considered.
 
    NOTE: In this subroutine, S2 = sigma^2 is defined following code SHARM and
    differs from that defined in ocean.phase.f by the factor of 2!
 
    NOTE: Eq.(2) defined in [1, Eq.(17)] and in [3, Eq.(87)] differs from that
    defined in  [2, Eq.(7)].
 
    The probability distribution function is given in [2, Eq.(4)]
 
    P = a^2/(pi*mu*S^2)*exp((1 - 2a)/S^2),                                   (4)
    a = (1 + cos(2h))/(mu + mu')^2,                                          (5)
    cos(2h) = mu*mu' - smu*smu'*cos(fi-fi'),                                 (6)
 
    where 2h is the angle complementary to the scattering angle (h is the angle
    of reflection = angle of incidence).
 
    NOTE: the signs in Eqs.(4)-(6) are given for MUI > 0, and MUR > 0 [2].
 
    The BRDF for the rough sea surface is defined as [2, Eq.(3)]
 
    WSRF = P * S.                                                            (7)
 
    In order to bring the boundary condition to the same form as over land, the
    BRDF is defined following SHARM [4, Appendix, the last equation of section
    'C.Ocean', unnumbered]
 
    WSRF = P * S * pi/mu.                                                    (8)
 
    NOTE: 1/pi in Eq.(4) and *pi in Eq.(8) are cancelled.
 
    NOTE: the typical expression for 'sigma' is [1, Eq.(15); 2, p.4250]
 
    S^2 = 0.00534*U,                                                         (9)
 
    another definition was used in [3, Eq.18]
 
    2*S^2 = 0.00512*U + 0.003,                                              (10)
 
    and in SCIATRAN code. Definitions Eq.(9) and (10) are equivalent, i.e. give
    the same result on output if S^2 = 2*S^2.

REFERENCESS:
    1. Nakajima T, Tanaka M, JQSRT(1983), V.29, N.6, pp.521-537.
    2. Gordon HR, Wang M, Appl.Opt.(1992), V.31, N.21, pp.4247-4260.
    3. Mishchenko MI, Travis LD, JGR(1997), V.102, D14, pp.16989-17013.
    4. Lyapustin AI, Appl.Opt.(2005), V.44, N.36, pp.7764-7772.
    5. Cox C, Munk W, JOSA (1954), V.44, N.11, pp.838-850.
	6. SORD: ROUGHSRF.f90

PROTOTYPE:
	void roughsrf1(double const, double const, double const, double const,
	               double const, double *, int const, double *);
------------------------------------------------------------------------------*/
//
//	LOCAL VARIABLES
	int iaz;
	double a, cir, c2inc, p, sgm, sgm2, sfi, sfr, sir, x, x2;
//-----------------------------------------------------------------------------
//
	sgm2 = S1 * wsp + S2;
	sgm = sqrt(sgm2);
//
	if (smui < NADIR) // incident ray shadow factor
		sfi = 0.0;
	else {
		x = cmui / smui / sgm;
		sfi = 0.5 * (exp(-x * x) / (SQRPI * x) - erfc(x));
	}
//
	if (smur < NADIR) // reflected ray shadow factor
		sfr = 0.0;
	else {
		x = -cmur / smur / sgm; // -cmur > 0
		sfr = 0.5 * (exp(-x * x) / (SQRPI * x) - erfc(x));
	}
//
	x = cmui - cmur; // -cmur > 0
	x2 = x * x;
	cir = -cmui * cmur; // -cmur > 0
	sir = smui * smur;
	for (iaz = 0; iaz < naz; iaz++) {
//      cos(Scattering angle {between v and -v'}) = cos(2*Angle_of_incidence)
		c2inc = cir - sir * caz[iaz];
		a = (1.0 + c2inc) / x2;
		p = a * a * exp((1.0 - 2.0 * a) / sgm2) / (cir * sgm2);
		wsrf[iaz] = p / (1.0 + sfi + sfr);
	} // for iaz
//
} // double roughsrf1
/*------------------------------------------------------------------------------
 18/01/27 - First created, tested vs ROUGHSRF.f90 for wsp = 0.1, 1.0, 3.0, 10.0
			sza[] = { 0.0, 45.0, 60.0, 75.0},
			vza[] = { 0.0, 30.0, 45.0, 60.0, 80.0},
			aza[] = { 0.0,   45.0, -45.0,  90.0, 135.0, 180.0,
			        225.0, -225.0, 270.0, 315.0, 360.0 };
			8 digits agreement
*/