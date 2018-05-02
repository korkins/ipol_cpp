# include <math.h>	/* sqrt */
//
void polqkm(int const m, int const k, double const *x, int const nx, double *q) {
/*------------------------------------------------------------------------------
PURPOSE:
	To compute the Qkm(x) plynomials for all k = m:K & the Fourier order m > 0.
	Qkm(x) = 0 is returned for all k < m.

INPUT:
	m    i[1]    The Fourier order, m > 0
	k    i[1]    Maximum order, k > 0
	x	 d[nx]   Abscissas, x = [-1:+1]
	nx   i[1]    Length of x

OUTPUT:
	q    d[(k+1)*nx]   Qkm(x) with k changing first

TREE:
	-

COMMENTS:
	Definition:

	Qkm(x) = sqrt[(k-m)!/(k+m)!]*Pkm,                                        (1)
	Pkm(x) = (1-x2)^(m/2)*(dPk(x)/dx)^m,                                     (2)

	where Pk(x) are the Legendre polynomials. Note, unlike in [2] (-1)^m is
	omitted in Qkm(x). Refer to [1-4] for details.

	For fast summation over k, this index changes first.

	Qkm(x) for a few initial values of m > 0 and k for testing:

	m = 1:
		Q01 = 0.0                                                // k = 0
		Q11 = sqrt( 0.5*(1.0 - x2) )                             // k = 1
		Q21 = 3.0*x[ix]*sqrt( (1.0 - x2)/6.0 )                   // k = 2
		Q31 = (3.0/4.0)*(5.0*x2 - 1.0)*sqrt( (1.0 - x2)/3.0 )    // k = 3

	m = 2:
		Q02 = 0.0                                                // k = 0
		Q12 = 0.0                                                // k = 1
		Q22 = 3.0/(2.0*sqrt(6.0))*(1.0 - x2);		             // k = 2
		Q32 = 15.0/sqrt(120.0)*x[ix]*(1.0 - x2);	             // k = 3
		Q42 = 15.0/(2.0*sqrt(360.0))*(7.0*x2 - 1.0)*(1.0 - x2)   // k = 4

	m = 3:
		Q03 = 0.0                                                // k = 0
		Q13 = 0.0												 // k = 1
		Q23 = 0.0                                                // k = 2
		Q33 = 15.0/sqrt(720.0)*(1.0 - x2)*sqrt(1.0 - x2);        // k = 3
		Q43 = 105.0/sqrt(5040.0)*(1.0 - x2)*x[ix]*sqrt(1.0 - x2) // k = 4

REFERENCESS:
	1. Gelfand IM et al., 1963: Representations of the rotation and Lorentz
       groups and their applications. Oxford: Pergamon Press.
    2. Hovenier JW et al., 2004: Transfer of Polarized Light in Planetary
       Atmosphere. Basic Concepts and Practical Methods, Dordrecht: Kluwer
       Academic Publishers.
    3. http://mathworld.wolfram.com/AssociatedLegendrePolynomial.html
	4. http://www.mathworks.com/help/matlab/ref/legendre.html

PROTOTYPE:
	void polqkm(int const, int const, double const *, int const, double *);
------------------------------------------------------------------------------*/
//
//	LOCAL VARIABLES
	int ik, im, ix, k1;
	double c0, m1, m2, m4, sqrx2;
//
//	LOCAL ARRAYS
	double *c1, *c2, *c3; // c1[k+1], c2[k+1], c3[k+1] - coefficients
//-----------------------------------------------------------------------------
//
//  For Qmm(x) - initial value for the recurrence relation
	c0 = 1.0;
	for (ik = 2; ik <= 2*m; ik+=2) c0 = c0 - c0/ik;
	c0 = sqrt(c0);
//
//	For the recurrence relation
	m2 = 1.0*m*m;
	m1 = m2 - 1.0;
	m4 = m2 - 4.0;
	k1 = k+1;
	c1 = new double [k1];
	c2 = new double [k1];
	c3 = new double [k1];
	for (ik = m+1; ik < k1; ik++) {
		c1[ik] = 2.0*ik - 1.0;
		c2[ik] = sqrt( (ik + 1.0)*(ik - 3.0) - m4 );
		c3[ik] = 1.0/sqrt( (ik + 1.0)*(ik - 1.0) - m1 );
	} // ik
//
//  Loop over abscissas, x
	for (ix = 0; ix < nx; ix++) {
//		Qkm(x) = 0 for k < m
		for (ik = 0; ik < m; ik++) q[ix*k1+ik] = 0.0;
//
//		Qmm(x)=c0*[sqrt(1-x2)]^m
		q[ix*k1+m] = c0;
		sqrx2 = sqrt( 1.0 - x[ix]*x[ix] );
		for (im = 1; im <= m; im++) q[ix*k1+m] *= sqrx2;
//
//		Q{k-1}m(x), Q{k-2}m(x) -> Qkm(x)
		for (ik = m+1; ik < k1; ik++)
			q[ix*k1+ik] = ( c1[ik]*x[ix]*q[ix*k1+ik-1] - 
							c2[ik]*q[ix*k1+ik-2] )*c3[ik];
	} // ix
//
	delete[] c1;
	delete[] c2;
	delete[] c3;
} // void polqkm(...)
/*------------------------------------------------------------------------------
 18/04/29 - double CONST *x
 18/03/17 - delete -> delete[]: http://www.cplusplus.com/doc/tutorial/dynamic/
 18/02/17 - First created and tested vs explicit expressions (see above) for
		    m = 1, 2, 3 and k = 0:4 for x= -1.0:0.001:1.0 and k_max = 1000.

		    k = 512 (in Fortrna 513), m = 256

		       x        POLQKM.f90               polqkm.cpp              Err,%
		    -1.00       0.000000000000000E+000  -0.00000000000000000	    -
		    -0.50 (!)  -2.601822304856592E-002  -0.026018223048565915    0.0
		     0.00       3.786666189291950E-002   0.037866661892919498    3E-13
			 0.25       9.592316443679009E-003   0.0095923164436790883  -8E-13
			 0.50 (!)  -2.601822304856592E-002  -0.026018223048565915    0.0
			 0.75      -2.785756308806302E-002  -0.027857563088062885    7E-13
			 1.00       0.000000000000000E+000   0.00000000000000000     -
*/