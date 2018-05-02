#include <cmath> /* abs, sqrt, acos */
#define TINY 1.0e-8             // to compare doubles
#define BR 1.0                  // Eq.(12) b/r
#define COS_MIN 0.03            // min value of zentih cos; depends on kv&kg
#define HB 2.0                  // Eq.(12) h/b
#define IP 0.31830988618379067  // 1/pi
#define PI2 1.5707963267948966  // pi/2
#define PI4 0.78539816339744831 // pi/4
//
void rtls1(double const kl, double const kv, double const kg,
		   double const smui, double const cmui,
		   double const smur, double const cmur,
	       double *caz, int const naz,
	       double *rtls) {
/*------------------------------------------------------------------------------
PURPOSE:
   To compute the RTLS BRDF for an array of azimuths

INPUT:
	kl     d[1]     lambertian kernel weight, kl = [0, 1]
	kv     d[1]     volumetric kernel weight, kv = [-1, 1] (?)
	kg     d[1]     geometric-optics kernel weight, kg = [-1, 1] (?)
	smui   d[1]     sin of incidence zenith, not equal 1.0
	cmui   d[1]     cos of incidence zenith, cmui > 0.0 (down) - as in RTE
	smur   d[1]     sin of reflection zenith, not equal 1.0
	cmur   d[1]     cos of reflection zenith, cmur < 0.0 (up) - as in RTE
	caz    d[naz]   cos(raz_rt); caz_rt = -caz_satellites (MODIS, POLDER)

OUTPUT:
	rtls   d[naz]   BRDF

COMMENTS:
    Neither this nor F90 subrotine was tested for BR & HB other than above.

	Refer to SHARM [3], BCondition.cpp -> LSRT for the original function.
	Unlike in SHARM, the following checks are NOT performed in rtls1:

		 cx for [-1...1] to compute x = acos(cx) and sx = sqrt(1.0 - cx*cx).

	THINKME1: for sqrt, one may use sx = sqrt(ABS(1.0 - cx*cx)).
	THINKME2: use if kv > 0 and if kg > 0 to save operations

	An example of possible LSRT parameters (SHARM):
		KL      KV      KG
		0.086   0.014   0.017
		0.175   0.108   0.041
		0.227   0.093   0.056
		0.330   0.053   0.066

	The RTLS code [1, 2] was adapted from the code SHARM [3]: BCondition.cpp,
	the LSRT subroutine. The following definition is used

	mu = cos(VZA) > 0, mu0 = cos(SZA) > 0, fi - azimuth                      (1)

	Note, min(mu) = min(mu0) = minimum_value [3]. BRDF is

	BRDF(SZA, VZA, AZA) = KL + KG*FG(mu0, mu, fi) + KV*FV(mu, mu, fi)        (2)

	FV = [(pi/2 - x)*cos(x) + sin(x)]/[mu0 + mu] - pi/4,                     (3)

	'Scattering' angle x is

	cos(x) = mu*mu0 + smu*smu0*cos(180-AZA)                                  (4)

	Note that sign at mu*mu0 here corresponds to the one used in SHARM's LSRT
	but differs from [3, Eq.(L-3)].
	???????????????????????????????????????????????????????????????????????????
	where AZA in RT methods differes from AZA for definition of surface
	reflection. Note, in Eq.(3) not-primed x, Eq.(4), is used. Standart RTLS

	FG = O(mu0, mu, fi) - 1/mu' - 1/mu0' + 0.5(1+cos(x'))/mu'/mu0',          (5)
	
	The O-function is

	O = 1/pi*(t - sin(t)cos(t))(1/mu' + 1/mu0')                              (7)

	and

	cos(t) = h/b sqrt( G'^2 + [tan(SZA')tan(VZA')sin(AZA')]^2 )/...
             ...(1/mu' + 1/mu0')                                             (8)

	with constraint that |cos(t)| <= 1. Finally,

	cos(x') = -cos(VZA')*cos(SZA') + sin(VZA')*sin(SZA')*cos(AZA')           (9)

	G = [tan2(SZA') + tan2(VZA') - 2*tan(SZA')*tan(VZA')*cos(180-AZA')]     (10)

	Note that '-' in Eq.(10) corresponds to the one used in SHARM, but differs
	from [3, Eq.(L-4)].

	The primed values are obtained using the following rule

	tan(ANG') = (b/r)*tan(ANG)                                              (11)

	The ratio of structural parameters is fixed following [2, P.981; 3, P.7770]

	h/b = 2, b/r = 1                                                        (12)

	Tables 1 and 2 from [1] give examples of input parameters. Note that f2
	[1, Eq.(8)] corresponds to Eq.(3) [3, Eq.(A16)] exactly except for a factor
	4/3pi. Thus, the values of k2 = [0.001...1.33] give a good idea of typical
	magnitude of the KV parameter. NOTE: k2 may exceed 1 in some cases!

	But f1 [1, Eq.(2)] differs from Eq.(5) [3, Eq.(A17)]. So, the values of
	k1 = [0...0.0073] gives only the idea for ratio between KG/KV as compared
	to k1/k2. From Tables 1 and 2 [1], k1/k2 << 1.

	It is noted in [3] after Eq.(A19) that FV and FG take both positive and
	negative values.

	CAUTION must be exercised when this model is used at zenith angles larger
	then 80 deg., when the BRDF may become negative or large [Ref.31 in 3].
	Note that MINC = 0.03 (see parameter below) corresponds to 88 deg.

	*** IMPORTANT: MINC depends on KV and KG and does not guarantee from
	*** negative values in all possible cases. Probably it would be better to
	*** replace negative values with the last positive one.

REFERENCESS:
    1. Roujean J-L, et al, JGR, 1992, V97, D18, P20455
    2. Lucht W, Schaaf, IEEE Trans Geosci Rem Sens, 2000, V38, N2, P977
    3. Lyapustin A, Appl Opt, 2005, V44, N36, P7764

PROTOTYPE:
	void rtls1(double const, double const, double const,
			   double const, double const,
			   double const, double const,
			   double *, int const,
			   double *);
------------------------------------------------------------------------------*/
//
//	LOCAL VARIABLES
	int iaz;
	double cc, ci, cr, ct, cx, g2, ic, o, si, sr, ss, sx, t, ti, ti2,
		   tr, tr2, tt, tt2, t2, x;
//-----------------------------------------------------------------------------
//
//  THINKME: using namespace std; maxnum = max(number1, number2) or C++ macros
	if (cmui < COS_MIN) ci = COS_MIN;
	else ci = cmui;
	if (-cmur < COS_MIN) cr = COS_MIN; // cmur < 0, cr > 0
	else cr = -cmur;
//
	si = smui;
	sr = smur;
	ti = si / ci;
	tr = sr / cr;
	cc = ci * cr;
	ss = si * sr;
	ic = 1.0 / (ci + cr);
//
//  THINKME: merge loops for fv & fg (note: same x in fv & fg; use x_v & x_g)
	for (iaz = 0; iaz < naz; iaz++) {
		cx = cc - ss * caz[iaz]; // "-" due to difference in definition of az
		sx = sqrt(1.0 - cx * cx);
		x = acos(cx);
		rtls[iaz] = kl + kv*(((PI2 - x) * cx + sx)*ic - PI4); // kl + kv*fv
	} // for iaz
//
//  Primed values: Eq.(10)
	ti *= BR;
	ti2 = ti * ti;
	ci = 1.0 / sqrt(ti2 + 1.0);
	si = sqrt(1.0 - ci * ci);
	tr *= BR;
	tr2 = tr * tr;
	cr = 1.0 / sqrt(tr2 + 1.0);
	sr = sqrt(1.0 - cr * cr);
	cc = ci * cr;
	ss = si * sr;
	ic = 1.0 / ci + 1.0 / cr; // ic redefined
	tt = ti2 + tr2;
	t2 = 2.0 * ti * tr;
	tt2 = ti2 * tr2;
//
	for (iaz = 0; iaz < naz; iaz++) {
		cx = cc - ss * caz[iaz];                    // - due to difference in az
		g2 = abs(tt + t2 * caz[iaz]);      // g*g, Eq.(8); note difference in az
		ct = HB * sqrt(g2 + tt2 * (1.0 - caz[iaz] * caz[iaz])) / ic;  // Eq.(7)
		if (ct > 1.0) ct = 1.0;           // ! SHARM [3], BCondition.cpp -> LSRT
		t = acos(ct);                // THINKME: check ct for [-1...1] as in [3]
		o = IP * (t - sqrt(1.0 - ct * ct) * ct) * ic;                  // Eq.(6)
		rtls[iaz] += kg * (o - ic + 0.5 * (1.0 + cx) / cc); // [kl+kv*fg]+kg*fg
	} // for iaz
//
} // double rtls1
/*------------------------------------------------------------------------------
 18/01/22 - First created, tested vs RTLS.f90 for [kl, kv, kg] = [0.3, 0.2, 0.1]
			sza[] = { 0.0, 45.0, 60.0, 75.0},
			vza[] = { 0.0, 30.0, 45.0, 60.0, 80.0},
			aza[] = { 0.0,   45.0, -45.0,  90.0, 135.0, 180.0,
			        225.0, -225.0, 270.0, 315.0, 360.0 };
			8 digits: 100*(1.0 - rtls/RTLS)% = 0.00000000.
			At some points, rtls.cpp = RTLS.f90 < 0.0

			Same for [kl, kv, kg] = [0.175, 0.108, 0.041]
*/