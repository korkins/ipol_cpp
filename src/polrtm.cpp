# include <math.h> /* sqrt, pow */
# define max(a, b) ( ((a)>(b)) ? (a):(b) ) /* macros: max of 2 integers */
//
void polrtm(int const m, int const k, double const *x, int const nx,
														 double *r, double *t) {
/*------------------------------------------------------------------------------
PURPOSE:
	To compute the Rkm(x) & Tkm(x) plynomials for all k = m:K & the Fourier
	order m > 0. Rkm(x) = Tkm(x) = 0 is returned for all k < max(m, 2).

INPUT:
	m    i[1]    The Fourier order, m > 0 (as in theory)
	k    i[1]    Maximum order, k > 0
	x	 d[nx]   Abscissas, x = [-1:+1]
	nx   i[1]    Length of x

OUTPUT:
	r, t   d[(k+1)*nx]   Rkm(x) & Tkm(x) with k changing first

TREE:
	-

COMMENTS:
	The matrix polynomial is defined following [1]:

                      | Qkm(x)  0       0       0     |
                      | 0       Rkm(x) -Tkm(x)  0     |
             Pkm(x) = | 0      -Tkm(x)  Rkm(x)  0     |,
                      | 0       0       0       Qkm(x)|

	where

             Rkm(x) = -1/2*i^m*(Pkm{+2} + Pkm{-2}),
             Tkm(x) = -1/2*i^m*(Pkm{+2} - Pkm{-2}),

	are real functions of -1 <= x <= 1 and i = sqrt(-1).

             Qkm(x) = sqrt[(l-m)!/(l+m)!]*Pkm,
             Pkm(x) = (1-x^2)^(m/2)(d/dmu)^m{Pl(x)},

	where Pk(x) are the ordinary Legendre polynomials of an order k. Note
	that (-1)^m is not used in Qkm. Generalized Legendre polynomials are [2]

             Pkm{n} = Akm{n}*(1 - x)^(-(n - m)/2)*(1 + x)^(-(n + m)/2)*
                      [d/dmu]^(k-n)[(1 - x)^(k-m)*(1 + x)^(k+m)],

	where n = -2, -0, +0, + 2, and

	Akm{n} = [(-1)^(k-m)*i^(n-m)/2^l/(k-m)!]*[(k-m)!(k+n)!/(k+m)!/(k-n)!]^1/2.

	This definition coincides with that widely used in literature [3-5].

	THINK ME: Include M=1 in the general case and remove IF (M==1) check.
	THINK ME: Move case M > 1 in front of M=1 as more frequent

	To test the subroutine use:

	a. Symmetry relation

		Rkm(-x) = +(-1)^(m+k)*Rkm(x)  - note 'plus'
		Tkm(-x) = -(-1)^(m+k)*Tkm(x)  - note 'minus'

	b. Explicit expressions for m = k:

		Rkk(x) = Ck*fk(x)*(1 + x2)
		Tkk(x) = 2*Ck*fk(x)*x

	where
	
		Ck = (+1/2^k)*sqrt(2k!/(k-2)!/(k+2)!)
		fk(x) = (1 - x2)^( (k-2)/2 )

	c. Explicit expressions for m=1, k=6:

			r = (c/sx)*(99*x7 - 201*x5 + 121*x3 - 19*x)
			t = (c/sx)*(33*x6 -  51*x4 +  19*x2 - 1)
		
	   where c= sqrt(10)/32, sx = sqrt(1 - x2), xn = x^n.

	d. Explicit expressions for m=2, k=6 (*** TO BE SIMPLIFIED ***):
		
		p1 = (1/64)*(1 + x)*(1 + x)*(495*x4 - 660*x3 + 90*x2 + 108*x - 17.0)
		p2 = (1/64)*(1 - x)*(1 - x)*(495*x4 + 660*x3 + 90*x2 - 108*x - 17.0)

		r = 0.5*(p1 + p2)
		t = 0.5*(p1 - p2)

	e. Explicit expressions for m=3, k=6 (*** TO BE SIMPLIFIED ***):

		c = 3/32; sx = sqrt( (1.0 - x)/(1.0 + x[ix]) )

		p1 =    c*sx*(55*x6 + 110*x5 + 5*x4 - 92*x3 - 31*x2 + 14*x + 3)
		p2 = -(c/sx)*(55*x6 - 110*x5 + 5*x4 + 92*x3 - 31*x2 - 14*x + 3)

		r = 0.5*(p1 + p2)
		t = 0.5*(p1 - p2)

	f. Explicit expressions for m=3, k=6 (*** TO BE SIMPLIFIED ***):

		c = -sqrt(30)/64

		p1 = c*(x - 1)*(33*x5 + 77*x4 + 34*x3 - 30*x2 - 19*x + 1);
		p2 = c*(x + 1)*(33*x5 - 77*x4 + 34*x3 + 30*x2 - 19*x - 1);

		r = 0.5*(p1 + p2)
		t = 0.5*(p1 - p2).

REFERENCESS:
	1. Siewert CE, JQSRT (2000), 64, pp.227-254
	2. Gelfand IM et al., 1963: Representations of the rotation and Lorentz
	   groups and their applications. Oxford: Pergamon Press.
	3. Hovenier JW et al., 2004: Transfer of Polarized Light in Planetary
       Atmosphere. Basic Concepts and Practical Methods, Dordrecht: Kluwer
       Academic Publishers.
    4. Rozanov VV and Kokhanovsky AA, Atm.Res., 2006, 79, 241.
    5. Y.Ota et al., JQSRT, 2010, 111, 878.

PROTOTYPE:
	void polrtm(int const, int const, double const *, int const, double *, double *);
------------------------------------------------------------------------------*/
//
//	LOCAL VARIABLES
	int ik, ix, k1, mm2;
	double c0, cx, ik1, ik2, m2, x2;
//
//	LOCAL ARRAYS
	double *c1, *c2, *c3; // c1[k+1], c2[k+1], c3[k+1] - coefficients
//------------------------------------------------------------------------------
//
	k1 = k+1;
	c1 = new double [k1];
	c2 = new double [k1];
	c3 = new double [k1];
//
	m2 = 1.0*m*m;
	mm2 = max(m, 2);
	for (ik = mm2+1; ik < k1; ik++) {
		ik1 = ik - 1.0;
		ik2 = ik1*ik1;
		c0 = ik/sqrt( (ik - 2.0)*(ik + 2.0)*(ik+m)*(ik-m) );
		c1[ik] = (2.0*ik - 1.0)*c0;
		c2[ik] = 2.0*m*(2.0*ik - 1.0)*c0/(ik1*ik);
		c3[ik] = sqrt( (1.0 - m2/ik2)*(ik2 - 4.0) )*c0;
	} // ik
//
	if (m != 1) { // in most cases m > 1
		c0 = sqrt( 0.03125*m*(m-1) ); // 0.03125 = 1/32
		for (ik = 2; ik < m; ik++) c0 = 0.5*c0*sqrt( m/(ik + 1.0) + 1.0 );
		for (ix = 0; ix < nx; ix++) {
//			k < m: Rkm(x)=Tkm(x)=0
			for (ik = 0; ik < m; ik++) {
				r[ix*k1+ik] = 0.0;
				t[ix*k1+ik] = 0.0;
			} // ik
//			k = m: Rmm(x) & Tmm(x) - explicit expressions
		    x2 = x[ix]*x[ix];
			r[ix*k1+m] = c0*pow( sqrt( 1.0 - x2 ), m-2 ); // <<< Inefficient
			t[ix*k1+m] = 2.0*x[ix]*r[ix*k1+m];
			r[ix*k1+m] = (1.0 + x2)*r[ix*k1+m];
//			k > m: recurrence
			for (ik = m+1; ik < k1; ik++) {
				cx = c1[ik]*x[ix];
				r[ix*k1+ik] = cx*r[ix*k1+ik-1] - c2[ik]*t[ix*k1+ik-1] -
							  c3[ik]*r[ix*k1+ik-2];
				t[ix*k1+ik] = cx*t[ix*k1+ik-1] - c2[ik]*r[ix*k1+ik-1] -
							  c3[ik]*t[ix*k1+ik-2];
			} // ik
		} // ix
	} else { // m == 1
		for (ix = 0; ix < nx; ix++) {
//			k = 0, 1: Tkm(x) = Rkm(x) = 0 for m = 1
			r[ix*k1] = 0.0;
			r[ix*k1+1] = 0.0;
			t[ix*k1] = 0.0;
			t[ix*k1+1] = 0.0;
//			k = 2: explicit expressions
			t[ix*k1+2] = -0.5*sqrt( 1.0 - x[ix]*x[ix] );
			r[ix*k1+2] = t[ix*k1+2]*x[ix];
//			k > 2: recurrence
			for (ik = 3; ik < k1; ik++) {
				cx = c1[ik]*x[ix];
				r[ix*k1+ik] = cx*r[ix*k1+ik-1] - c2[ik]*t[ix*k1+ik-1] -
							  c3[ik]*r[ix*k1+ik-2];
				t[ix*k1+ik] = cx*t[ix*k1+ik-1] - c2[ik]*r[ix*k1+ik-1] -
							  c3[ik]*t[ix*k1+ik-2];
			} // ik
		} // ix
	} // if m != 1
//
	delete[] c1;
	delete[] c2;
	delete[] c3;
} // void polrtm(...)
/*------------------------------------------------------------------------------
 18/04/29 - double CONST *x
 18/03/17 - delete -> delete[]: http://www.cplusplus.com/doc/tutorial/dynamic/
 17/04/20 - Symmetry relations checked for x = 0.001:0.001:0.999999;
		    m = 1:1000; k = 0:1000.

 17/04/19 -  Tested against all the cases (c-f). Tested for m=999,k=1000; m=256,
             k=512. For m=128, k=255:

		       x     R.f90				  R.f90-R	T.f90                  T.f90-T
		     -0.95  -1.05268034363584E-17   -2.0E-31   1.02235651027377E-17   0.0E+00
		     -0.90  -9.21141805327410E-05    1.0E-18   7.39638912759868E-05   1.4E-18
		     -0.80  -9.61876451104626E-05	 -5.7E-16   7.98048018891205E-02  -4.0E-16
		     -0.70   2.38266484684450E-04	 -2.5E-16   6.94824164795929E-02  -1.6E-15
		     -0.60   2.05934563634704E-04	  3.2E-16   6.17205086475275E-02   3.0E-16
		     -0.50  -1.18563285251866E-02	  6.0E-16   4.44739866048293E-02   5.0E-16
		     -0.40  -1.98699738813941E-02	 -4.0E-16  -2.61984592363673E-02   9.0E-16
		     -0.30   2.13322312553310E-02	 -5.0E-16  -2.52825563533065E-02  -1.1E-15
		     -0.20  -1.09985749531716E-02	  1.1E-15   4.32761863993168E-02   3.1E-16
		     -0.10   3.82484302902520E-03	 -9.7E-16  -4.63484667213980E-02  -3.0E-16
		      0.00   0.00000000000000E+00	  8.2E-16   4.65292818279906E-02   9.7E-17
		      0.10  -3.82484302902586E-03	 -1.3E-15  -4.63484667213979E-02   0.0E+00
			  0.20   1.09985749531716E-02	  3.0E-16   4.32761863993168E-02  -1.0E-15
			  0.30  -2.13322312553310E-02	 -6.0E-16  -2.52825563533065E-02   1.3E-15
			  0.40   1.98699738813946E-02	  0.0E+00  -2.61984592363664E-02   1.0E-16
			  0.50   1.18563285251866E-02	  3.0E-16   4.44739866048293E-02  -1.4E-15
			  0.60  -2.05934563634300E-04	  6.8E-16   6.17205086475270E-02   9.7E-17
			  0.70  -2.38266484684450E-04	  3.3E-16   6.94824164795929E-02  -1.4E-15
			  0.80   9.61876451104626E-05	 -1.5E-15   7.98048018891205E-02  -3.4E-15
			  0.90   9.21141805327449E-05	 -5.8E-18   7.39638912759902E-05  -3.3E-18
			  0.95   1.05268034363584E-17	  2.0E-31   1.02235651027377E-17   0.0E+00

 17/04/17 - Tested against explicit Rmm(x) and Tmm(x) for m = k = 2:20.
		    x = -1.0:0.1:+1 (nx=21).
 17/04/14 - First created, compiled. Not yet tested.
*/