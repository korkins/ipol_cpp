void polleg(int const k, double const *x, int const nx, double *p) {
/*------------------------------------------------------------------------------
PURPOSE:
	To compute the Legendre polynomials, P{k}(x)
    
INPUT:
	k    i[1]    Maximum order, k > 0
	x	 d[nx]   Abscissas, x = [-1:+1]
	nx   i[1]    Length of x

OUTPUT:
	p    d[(k+1)*nx]   Legendre polynomials

TREE:
	-

COMMENT:
	Pk(x) = Qkm(x) for m = 0. The Bonnet recursion formula [1, 2]:

	(k+1)P{k+1}(x) = (2k+1)*P{k}(x) - k*P{k-1}(x),                           (1)

	where k = 0:K, P{0}(x) = 1.0, P{1}(x) = x.

	For fast summation over k, this index changes first.

REFERENCESS:
	1. https://en.wikipedia.org/wiki/Legendre_polynomials
	2. http://mathworld.wolfram.com/LegendrePolynomial.html

PROTOTYPE:
	void polleg(int const, double const *, int const, double *);
------------------------------------------------------------------------------*/
//
//	LOCAL VARIABLES
	int ik, ix, k1;
	double ki;
//
//	LOCAL ARRAYS
	double *c1, *c2; // c1[k+1], c2[k+1] - temporary arrays for coefficients
//-----------------------------------------------------------------------------
//
//	Precompute coefficients for the recursion: C1(0:1) & C2(0:1) are not used
	k1 = k+1;
	c1 = new double [k1];
	c2 = new double [k1];
	for (ik = 2; ik < k1; ik++) {
		ki = 1.0/ik;
		c1[ik] = 2.0 - ki;
		c2[ik] = 1.0 - ki;
	} // ik
//
//  Recursion
	for (ix = 0; ix < nx; ix++) {   // abscissa
		p[ix*k1] = 1.0;
		p[ix*k1+1] = x[ix];
		for (ik = 2; ik < k1; ik++) // order
			p[ix*k1+ik] = c1[ik]*x[ix]*p[ix*k1+ik-1] - c2[ik]*p[ix*k1+ik-2]; 
	} // ix
//
	delete[] c1;
	delete[] c2;
} // void polleg(...)
/*------------------------------------------------------------------------------
 18/04/28 - double CONST *x is now used
 18/03/17 - delete -> delete[]: http://www.cplusplus.com/doc/tutorial/dynamic/
 18/02/10 - First created and tested vs excplicit form of Pk(x) for k = 0, 1, 2,
		    3, 6, 7, 10, x = -1:0.001:1 (nx=2001) and k = 1000.
		    t_cpu = 0.01s for all P[1001*2001] values.
*/