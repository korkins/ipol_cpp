void polleg(int const nk, double const x, double *p) {
/*------------------------------------------------------------------------------
PURPOSE:
	To compute the Legendre polynomials, P{k}(x), for k = 0...nk-1.
    
INPUT:
	nk   i[1]    Number of polynomials, nk > 1
	x	 d[1]    Abscissa, x = [-1, +1]

OUTPUT:
	p    d[nk]   Legendre polynomials

TREE:
	-

COMMENT:
	Pk(x) = Qkm(x) for m = 0. The Bonnet recursion formula [1, 2]:

	(k+1)P{k+1}(x) = (2k+1)*P{k}(x) - k*P{k-1}(x),                           (1)

	where k = 0:K, P{0}(x) = 1.0, P{1}(x) = x.

	In Fortran and preveious version of this subroutine, an array of x was
	assumed on input with (2k+1)/(k+1) & k/(k+1) precomputed. This might save
	time but it is nothing compared to singular value decomposition. If x is not
	an array, the subroutine is shorter and not much slower:
	2000 values of x, 1000 values of k, win32-release, 0.03(new) vs 0.02(old)

REFERENCESS:
	1. https://en.wikipedia.org/wiki/Legendre_polynomials
	2. http://mathworld.wolfram.com/LegendrePolynomial.html

PROTOTYPE:
	void polleg(int const, double const, double *);
------------------------------------------------------------------------------*/
//
//	LOCAL VARIABLES
	int ik;
//-----------------------------------------------------------------------------
//
	p[0] = 1.0;
	p[1] = x;
	for (ik = 2; ik < nk; ik++)
		p[ik] = (2.0 - 1.0 / ik) * x * p[ik-1] - (1.0 - 1.0 / ik) * p[ik-2];
//
} // void polleg(...)
/*------------------------------------------------------------------------------
 16Mar17 - On input x is now a variable; int const nx is removed;
 10Feb17 - First created and tested vs excplicit form of Pk(x) for k = 0, 1, 2,
		   3, 6, 7, 10, x = -1:0.001:1 (nx=2001) and k = 1000.
		   t_cpu = 0.01s for all P[1001*2001] values.
*/