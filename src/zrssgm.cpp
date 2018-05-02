# include <math.h> /* sqrt */
//
void zrssgm(int const m, double const def, double *mu0, int const nsz,
			double *mug, int const ng2, double *zr1, double *zr2, double *zr3) {
/*------------------------------------------------------------------------------
PURPOSE:
	To compute the 1st column of the Rayleigh phase matrix for the Fourier order
	m > 0 at the Gauss nodes

INPUT:
	m     i(1)     The Fourier moment, m = 1 or 2
	def   d(1)     Depolarization factor, def = [0 ... ~6/7]
	mu0   d(nsz)   cos(sza) > 0; avoid mu0 == 1.0 (replace with 0.999...999)
	nsz   i(1)     Length of mu0
	mug   d(ng2)   Gauss nodes
	ng2   i(1)     Length of mug

OUTPUT:
	zr1, zr2, zr3   d(ng2*nsz)   Elements of the phase matrix; ig=1:ng2 runs 1st

TREE:
	-

COMMENTS:
	Theoretical background for this subroutine is described in [1, p.93, Section
	3.4.4]. Note that polynomials used in [1] and in this subroutine differ by a
	factor of sqrt((k+m)!/(k-m)!). Nevertheless, the final result for the phase
	matrix, W [1, p.92, Eq.(3.128)], coincide exactly with the result of this
	subroutine.

	The 1st column of the phase matrix is

				| Qkm(mu) a1k Qkm(mu0)|
	ZR0M = sum( |-Rkm(mu) b1k Qkm(mu0)| ),                                   (1)
				| Tkm(mu) b1k Qkm(mu0)          |

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

    m = 1:                                                                   (6)
        Q01 = 0; Q11 = sqrt(1 - mu*mu)/2; Q21 = 3*mu*sqrt((1-mu*mu)/6)
        R01 = 0; R11 = 0; R21 = -mu*sqrt(1 - mu*mu)/2
        T01 = 0; T11 = 0; T21 = -sqrt(1 - mu*mu)/2
    m = 2:                                                                   (7)
        Q02 = 0; Q12 = 0; Q22 = 3*(1-mu*mu)/(2*sqrt(6))
        R02 = 0; R12 = 0; R22 = (1 + mu*mu)/4
        T02 = 0; T12 = 0; T22 = mu/2

REFERENCESS:
	1. HovenierJW et al. Transfer of polarized light, Kluwer 2004
	2. Bates DR, 1984, Planet. Space Sci., V32,pp.785–790

PROTOTYPE:
	void zrssgm(int const, double const, double *, int const, double *,
	            int const, double *, double *, double *)
------------------------------------------------------------------------------*/
//
//	LOCAL VARIABLES
	int
		ig, // Loop index over Gauss nodes
		is; // Loop index over Solar cosines
	double
		df,	// Eq.(4)
		d1, // temp var
		d2, // temp var
		y0; // 1 - 3*mu0^2
//
//	LOCAL ARRAYS
	double
		*x2; // x2[ng2] = mug*mug
//------------------------------------------------------------------------------
//
	df = (1.0 - def)/(1.0 + 0.5*def);
//
	x2  = new double [ng2];
	for (ig = 0; ig < ng2; ig++)
		x2[ig] = mug[ig]*mug[ig];
//
	switch(m) {
		case 1: // [1, p.94, Eq.(3.138)]
			df *= 0.75;										// df*3/4
			for (is = 0; is < nsz; is++) {
				d1 = df*mu0[is];
				y0 = 1.0 - mu0[is]*mu0[is];
				for (ig = 0; ig < ng2; ig++) {
					zr3[is*ng2+ig] = d1*sqrt( (1.0 - x2[ig])*y0 );
					zr1[is*ng2+ig] = zr2[is*ng2+ig] = zr3[is*ng2+ig]*mug[ig];
					zr3[is*ng2+ig] = -zr3[is*ng2+ig];
				} // ig
			} // is
		case 2: // [1, p.94, Eq.(3.140)]
			df *= 0.1875;									// df*3/16 
			for (is = 0; is < nsz; is++) {
				d1 = df*(1.0 - mu0[is]*mu0[is]);
				d2 = 2.0*d1;
				for (ig = 0; ig < ng2; ig++) {
					zr1[is*ng2+ig] = d1*(1.0 - x2[ig]);
					zr2[is*ng2+ig] = -d1*(1.0 + x2[ig]);
					zr3[is*ng2+ig] = d2*mug[ig];
				} // ig
			} // is
	} // switch(m)
//
	delete[] x2;
} // void zrssgm(...)
/*------------------------------------------------------------------------------
 18/03/17 - delete -> delete[]: http://www.cplusplus.com/doc/tutorial/dynamic/
 17/05/02 - First created and launched - compiled ok. Not yet tested.
*/