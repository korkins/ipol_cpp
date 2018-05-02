void zrssg0(double const def, double *mu0, int const nsz, double *mug,
			int const ng2,  double *zr1, double *zr2) {
/*------------------------------------------------------------------------------
PURPOSE:
	To compute the 1st column of the Rayleigh phase matrix for the Fourier order
	m = 0 at Gauss nodes

INPUT:
	def   d(1)     Depolarization factor, def = [0 ... ~6/7]
	mu0   d(nsz)   cos(sza) > 0
	nsz   i(1)     Length of mu0
	mug   d(ng2)   Gauss nodes
	ng2   i(1)     Length of mug

OUTPUT:
	zr1, zr2   d(ng2*nsz)   Elements of the phase matrix, ig = 1:ng2 changes 1st

TREE:
	-

COMMENTS:
	Theoretical background for this subroutine is described in [1, p.93, Section
	3.4.4]. Note that polynomials used in [1] and in this subroutine differ by a
	factor of sqrt((k+m)!/(k-m)!). Nevertheless, the final result for the phase
	matrix, W [1, p.92, Eq.(3.128)], coincide exactly with the result of this
	subroutine.

	The 1st column of the phase matrix is

				| Qk0(mu) a1k Qk0(mu0)|
	ZR0M = sum( |-Rk0(mu) b1k Qk0(mu0)| ),                                   (1)
				|          0          |

	where summation is taken over k = 0, 1, 2 (3 values). The expansion moments
	are [1, p.56, Section 2.9]

	a1k = [1  0   DF/2],                                                     (2)
	b1k = [0  0  (sqrt(6)/2)DF],                                             (3)

	where

	DF = (1 - d)/(1 + d/2),                                                  (4)

	and d is the depolarization factor on input; d = 0 (D = 1) corresponds to
	pure Rayleigh scattering without depolarization.

	d = Il(SCAT_ANG = 90)/Ir(SCAT_ANG = 90)                                  (5)

	Range of the depolarization factor, d=def, is [0...6/7] [1, p.56]. It is
	often assumed spectrally independent [2, p.3494].

	Explicit values for polynomilas involved in computations are

	m = 0:                                                                   (6)
		Q00 = 1; Q10 = mu; Q20 = (3*mu*mu - 1)/2
		R00 = 0; R10 = 0; R20 = sqrt(3/8)*(1 - mu*mu)

REFERENCESS:
	1. HovenierJW et al. Transfer of polarized light, Kluwer 2004
	2. Bates DR, 1984, Planet. Space Sci., V32,pp.785–790

PROTOTYPE:
	void zrssg0(double const, double *, int const, double *, int const,
	            double *, double *)
------------------------------------------------------------------------------*/
//
//	LOCAL VARIABLES
	int
		ig, // Loop index over Gauss nodes
		is; // Loop index over Solar cosines
	double
		df, // Eq.(4)
		d1, // df/8
		d2, // -3*df/8
		y0; // 3*mu0^2 - 1
//
//	LOCAL ARRAYS
	double *x2; // x2[ng2] = mug*mug
//------------------------------------------------------------------------------
//
	df = (1.0 - def)/(1.0 + 0.5*def);
	d1 =  0.125*df;
	d2 = -0.375*df;
//
	x2  = new double [ng2];
	for (ig = 0; ig < ng2; ig++)			// over nodes
		x2[ig] = mug[ig]*mug[ig];
//
	for (is = 0; is < nsz; is++) {			// over mu0
		y0 = 3.0*mu0[is]*mu0[is] - 1.0;
		for (ig = 0; ig < ng2; ig++) {		// over nodes
			zr1[is*ng2+ig] = 1.0 + d1*(3.0*x2[ig] - 1.0)*y0;
			zr2[is*ng2+ig] = d2*(1.0 - x2[ig])*y0;
		} // ig
	} // is
//
	delete[] x2;
} // void zrssg0(...)
/*------------------------------------------------------------------------------
 18/03/17 - delete -> delete[]: http://www.cplusplus.com/doc/tutorial/dynamic/
 17/05/02 - First created and launched - compiled ok. Not yet tested.
*/