#include <cmath>     /* cos, abs */
#define YEPS 3.0e-14 // to compare doubles
#define PI   3.1415926535897938
//
void gauszw(double const x1, double const x2, int const n,
	        double *x, double *w) {
/*------------------------------------------------------------------------------
PURPOSE:
	To compute the abscissas, x in [x1:x2] and weights, w, for the n-th order
	Gaussian-Legendre quadrature.

INPUT:
	x1, x2   d[1]   interval
	n	     i[1]   order of the quadrature

OUTPUT:
	x, w   d[n]   zeros and weights

COMMENTS:
	This subroutine was cretaed from [1], idicies modified: 0- vs 1-based.
    Note that computation time for N ~ 1000-10000 depends on how small the
    small number YEPS is. Example of the output data for N = 8 and [X1, X2]=
    [-1,+1] is (note the symmetry relation):
                         xi            wi
              i=1:  -0.96028986    0.10122854
              i=2:  -0.79666648    0.22238103
              i=3:  -0.52553241    0.31370665
              i=4:  -0.18343464    0.36268378
              i=5:   0.18343464    0.36268378
              i=6:   0.52553241    0.31370665
              i=7:   0.79666648    0.22238103
              i=8:   0.96028986    0.10122854

	Refer to [2, 3] for details. Gaussian numerical integration is

              int(f(x), -1..1) ~= sum(wj*f(xj), j = 1..n).

    Double-Gaussian [4] scheme for an even n is

              int(f(x), -1..1) = int(f(x), -1..0) + int(f(x), 0..1) =
    ~= sum(w2j*f(x2j), j=1..n/2, x2j<0) + sum(w2j*f(x2j), j=1..N/2, x2j>0).

    The relations between Gaussian and Double-Gaussian [4] weights and zeros are

              w2(+j) = 0.5wj, x2(+j) = 0.5(xj + 1) > 0,
              w2(-j) = w2(j), x2(-j) = 0.5(xj - 1) < 0,

    where xj and wj are computed using Gauss scheme (i.e. this subroutine) and
    N/2 knots.

	On input, x = [-1.0, +1.0] - yields ordinary qudrature, x = [ 0.0, +1.0] -
	yields double-Gaussian qudrature.

REFERENCES:
    1. GAUSZW.f90
	2. https://en.wikipedia.org/wiki/Gaussian_quadrature
	3. http://mathworld.wolfram.com/GaussianQuadrature.html
	4. Sykes J. Approximate integration of the equation of transfer,
       Mon.Not.Roy.Astron.Soc, 1951, 11, 377.

PROTOTYPE:
	void gauszw(double const, double const, int const, double *, double *);
------------------------------------------------------------------------------*/
//
//	LOCAL VARIABLES
	int m, i, j;
	double yxm, yxl, yz, yp1, yp2, yp3, ypp, yz1;
//-----------------------------------------------------------------------------
//
	m = (n + 1)/2;
	yxm = 0.5 * (x2 + x1);
	yxl = 0.5 * (x2 - x1);
	for (i = 0; i < m; i++) {
//      yz - initial approximation for i-th root
		yz = cos(PI * (i + 0.75) / (n + 0.5)); // i_Fortran vs i_C; see j below
		do {
			yp1 = 1.0; // Legendre polynomial at z
			yp2 = 0.0; // Legender polynomial of the previous order
			for (j = 0; j < n; j++) {
				yp3 = yp2;
				yp2 = yp1;
				yp1 = ( (2.0 * j + 1.0) * yz * yp2 - j * yp3 ) / (j + 1);
			} // for j
			ypp = n * (yz * yp1 - yp2) / (yz * yz - 1.0); // dLegPol(x)/dx
			yz1 = yz;
			yz = yz1 - yp1 / ypp; // Find the root using Newton method
			if (abs(yz - yz1) < YEPS) break; // from do...while loop
		} while (1); // 1 = true
		x[i] = yxm - yz * yxl;
		x[n - 1 - i] = yxm + yxl * yz;
		w[i] = 2.0 * yxl / ((1.0 - yz * yz) * ypp * ypp);
		w[n - 1 - i] = w[i];
	} // i
} // void gauszw(...)
/*------------------------------------------------------------------------------
 18/02/18 - Tested vs GAUSZW.f90 for n=128 & x=[-1:1] and n=13 & x=[0:pi].
            Print out 8 digits - ok.
*/