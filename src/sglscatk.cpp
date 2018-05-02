# include <math.h> /* fabs, exp, sqrt */
# define TINY 1.0e-8 // to compare doubles
# define C22 0.61237243569579452454932101867647 // 3/(2*sqrt(6)) for Q22(x)
//
void rotator2(double const, double const, double const, double const,
	          double const, double const, double &, double &);
//
void sglscatk(double const *mu, double const *smu, int const nmu, int const nup,
			  double const *mu0, double const *smu0, int const nmu0,
	          double const *caz, double const *saz, int const naz,
			  double const tau0, double const *tau, int const nlr,
			  double const *a1k, double const *b1k, int const nk,
	          double *i1, double *q1, double *u1) {                 // <- output
/*------------------------------------------------------------------------------
PURPOSE:
	To compute the single scattering path radiance on TOA and BOA using moments
	of the phase matrix (hence 'k')

INPUT:
	mu    d[nmu]      cos(vza), mu > 0 for downward
	smu   d[nmu]      sin(vza)
	nmu   i[1]        Number of view azimuths
	nup   i[1]        Number of upward directions
	mu0   d[nmu0]     cos(sza) > 0
	smu0  d[nmu0]     sin(sza)
	nmu0  i[1]        Number of solar zeniths
	caz   d[naz]      cos(aza), forward scattering: aza = 0
	saz   d[naz]      sin(aza)
	naz   i[1]        Number of azimuths
	tau0  d[1]        Total OT
	tau   d[nlr]      Optical thickness of layers
	nlr   i[1]        Number of optical layers
	a1k   d[nk*nlr]   ssa/2 & 2k+1 included: a1k[0] = ssa/2
	b1k   d[nk*nlr]   b1k[0:1] = 0.0, ssa/2 & 2k+1 included
	nk    i[1]        Number of moments, nk > 0

OUTPUT:
	i1, q1, u1   d[nmu*nmu0*naz]   the Stokes vector

COMMENTS:
	On input, nk is the fast dimension of a1k, b1k. On output, i1[mu->mu0->aza],
	where mu changes first (fastest). Upward directions, mu < 0 (if any), must
	come first.

REFERENCESS:
	1. Hovenier JW et al., 2004: Transfer of Polarized Light in Planetary
	   Atmosphere. Basic Concepts and Practical Methods, Dordrecht: Kluwer
	   Academic Publishers.

PROTOTYPE:
	void sglscatk(double *, double *, int const, int const, double *, double *,
		          int const, double *, double *, int const, double const,
				  double *, int const, double *, double *, int const,
				  double *, double *, double *)
------------------------------------------------------------------------------*/
//
//	LOCAL VARIABLES
	int iaz, ik, ilk, ilr, imu, imu0, isca, ixa, ixl, ix0;
	double cs, sn, cs0, sn0, x, x2, c0, ti, u, tai, tbi, ca, cb, sa, i1s, q1s,
		   pa1, pa2, pa3, pb1, pb2, pb3, a1, b1, et, s2x, c2x;
//
//	LOCAL ARRAYS
	double *ca1, *ca2, *cb1, *cb2, *e, *ta, *tb;
//-----------------------------------------------------------------------------
//
	e = new double[nlr * nmu * nmu0]; // e[nlr -> nmu -> nmu0]
	ta = new double[nlr];
	tb = new double[nlr];
//
//  OT above, ta, and below, tb, layer ilr
	ta[0] = 0.0;
	tb[0] = tau0 - tau[0];
	for (ilr = 1; ilr < nlr; ilr++) {
		ta[ilr] = ta[ilr - 1] + tau[ilr - 1];
		tb[ilr] = tb[ilr - 1] - tau[ilr];
	}
	tb[nlr - 1] = 0.0; // to avoid tb = -tiny
//
	for (imu0 = 0; imu0 < nmu0; imu0++) {
		cs0 = mu0[imu0];
		ix0 = imu0 * nmu;
//
//		Up: mu < 0 [1, p101, Eq.(4.22)]
		for (imu = 0; imu < nup; imu++) {
			ixl = (ix0 + imu) * nlr;
			cs = mu[imu];
			c0 = cs0 / (cs0 - cs);
			u = 1.0 / cs - 1.0 / cs0;
			for (ilr = 0; ilr < nlr; ilr++)
				e[ixl + ilr] =
					c0 * (1.0 - exp(tau[ilr] * u)) * exp(ta[ilr] * u);
		} // for imu (up)
//
//		Down: mu > 0 [1, p101, Eq.(4.21)]
		for (imu = nup; imu < nmu; imu++) {
			ixl = (ix0 + imu) * nlr;
			cs = mu[imu];
			if (fabs(cs0 - cs) < TINY) // lim(mu -> mu0)
				for (ilr = 0; ilr < nlr; ilr++) {
					ti = tau[ilr];
					e[ixl + ilr] =
						exp(-(ti + ta[ilr]) / cs0 - tb[ilr] / cs) * ti / cs0;
				} // for ilr
			else { // abs(cs0 - cs) > TINY
				c0 = cs0 / ( cs0 - cs);
				for (ilr = 0; ilr < nlr; ilr++) {
					ti = tau[ilr];
					tai = ta[ilr];
					tbi = tb[ilr];
					e[ixl + ilr] = c0 * ( exp( -(ti + tai) / cs0 - tbi / cs ) -
										  exp( -(ti + tbi) / cs - tai / cs0 ) );
				} // for ilr
			} // if |cs0 - cs| < TINY
		} // for imu (dn)
	} // for imu0
//
	ca1 = new double[nk];
	ca2 = new double[nk];
	cb1 = new double[nk];
	cb2 = new double[nk];
//
//  Precompute coefficients for polynomial recurrence [1, p207, Appendix B]
//  ca1, ca2, cb1, cb2 are not used for ik = 0, 1, 2 
	for (ik = 3; ik < nk; ik++) {
		ca = 1.0 / ik;
		ca1[ik] = 2.0 - ca;
		ca2[ik] = 1.0 - ca;
		cb = 1.0 * ik * ik - 4.0;
		cb1[ik] = (2.0 * ik - 1.0) / sqrt(cb);
		cb2[ik] = sqrt( (ik + 1.0) * (ik - 3.0) / cb ); 
	} // for ik
//
	for (iaz = 0; iaz < naz; iaz++) {
		ixa = iaz * nmu0 * nmu;
		ca = caz[iaz];
		sa = saz[iaz];
		for (imu0 = 0; imu0 < nmu0; imu0++) {
			ix0 = imu0 * nmu;
			cs0 = mu0[imu0];
			sn0 = smu0[imu0];
			for (imu = 0; imu < nmu; imu++) {
				ixl = (ix0 + imu) * nlr;
				cs = mu[imu];
				sn = smu[imu];
				x = cs * cs0 + sn * sn0 * ca; // [1, p70, Eq.(3.19)]
				x2 = x * x;
//
//				Accumulation over layers
				i1s = 0.0;
				q1s = 0.0;
				for (ilr = 0; ilr < nlr; ilr++) {
					ilk = ilr * nk;
					pa1 = x;
					pa2 = 1.5 * x2 - 0.5;
					a1 = a1k[ilk] + a1k[ilk + 1] * x + a1k[ilk + 2] * pa2; // [1, p53, Eq.(2.152)]
					pb1 = 0.0;
					pb2 = C22 * (1.0 - x2);
					b1 = -b1k[ilk + 2] * pb2; // [1, p53, Eq.(2.156)];  Pk20 = (-1)*Qk2
//
//					Compute [11] and [21] elements of the phase matrix at given x
					for (ik = 3; ik < nk; ik++) {
						pa3 = ca1[ik] * x * pa2 - ca2[ik] * pa1;
						pb3 = cb1[ik] * x * pb2 - cb2[ik] * pb1;
						a1 += a1k[ilk + ik] * pa3; // [1, p53, Eq.(2.152)]
						b1 -= b1k[ilk + ik] * pb3; // [1, p53, Eq.(2.156)];  Pk20 = (-1)*Qk2
//						Redefine for recurrence
						pa1 = pa2;
						pa2 = pa3;
						pb1 = pb2;
						pb2 = pb3;
					} // for ik
//
					et = e[ixl + ilr];
					i1s += a1 * et;
					q1s += b1 * et;
				} // for ilr
//
//				Apply rotation for the given scattering geometry
				isca = ixa + ix0 + imu;
				i1[isca] = i1s;
				rotator2(sn0, cs0, sn, cs, sa, ca, s2x, c2x);
				u1[isca] = q1s * s2x; // [1, p69, Eqs.(3.9-10)]
				q1[isca] = q1s * c2x; // [1, p69, Eqs.(3.9-10)]
			} // for imu
		} // for imu0
	} // for iaz
//
	delete[] e;
	delete[] ta;
	delete[] tb;
	delete[] ca1;
	delete[] ca2;
	delete[] cb1;
	delete[] cb2;
} // void sglscatk(...)
/*------------------------------------------------------------------------------
 18/04/28 - in rotator2: rot[] -> &s2x, &c2x; double CONST *array is now used.
 18/04/16 - SORD test 048 = IPRT B2: Rayleigh + absorption, 30 optical layers
			xi_i1, _q1, _u1 ~ e-10 (absolute deviation); all geometries.
 18/04/15 - First created and tested vs Wauben WMF, de Haan JF, Hovenier JW,
            1993: Low orders of scattering in plane-parallel homogeneous
			atmopshere, Astron.Astrophys., v.276, pp.589-602.
			See p.595, Section 6.1 for definition of the case.
			See p.596, Table 1, column "First" for the result.

			Moments: de Rooij WA and van der Stap CCAH, 1984: Expansion of Mie
			scattering matrices in generalized spherical functions,
			Astron.Astrophys., v.131, pp.237-248.
			See Table 4, p.243, model C.

			Tested for 1 layer and 4 layers of the dame tau0.
			aza = 0*, 90*, 180*
			mu0 = 1.0, 0.75, 0.6*, 0.4, 0.2
			mu = -1.0, -0.6*, -0.2, 0.1, 0.6*, 0.8, 1.0
			tau = [0.05, 0.03, 0.02, 0.1]
			
			*Benmchmark from Wauben et al:

            mu0   aza    mu    I               Q               U
			0.6   0.00  -0.6   1.2874000E-02  -1.0560000E-03   0.0000000E+00
			0.6   0.00   0.6   1.8164380E+00   0.0000000E+00   0.0000000E+00
			0.6  90.00  -0.6   4.6300000E-03   4.1600000E-04  -7.8000000E-04
			0.6  90.00   0.6   1.5914000E-02   4.9000000E-04  -9.1800000E-04
			0.6 180.00  -0.6   7.6640000E-03   0.0000000E+00   0.0000000E+00
			0.6 180.00   0.6   4.7870000E-03  -8.7500000E-04   0.0000000E+00
*/