# include <math.h>	/* sqrt */
# define C22 0.61237243569579452454932101867647 // 3/(2*sqrt(6)) for Q22(x)
//
void sumka1b1(double *x, int const nx, double *a1k, double *b1k,
	          int const nk, double *a1, double *b1) {
/*------------------------------------------------------------------------------
PURPOSE:
	To compute two elements of the scattering matrix, a1 & b1, at x.
 
INPUT:
	x     d[nx]   abscissae, [-1:+1], - cos(scattering_angle)
	nx    i[1]    length of x
	a1k	  d[nk]   a1 expansion moments, 2k+1 included
	b1k	  d[nk]   b1 expansion moments, 2k+1 included; b1k[0:1] = 0.0
	nk    i[1]    number of expansion moments, nk > 2

OUTPUT:
	a1, b1   d[nx]   a1(x) and b1(x)

TREE:
	-

COMMENTS:
	This subroutine was created from [1]. On input, a1k and b1k are scaled by a
	factor of (2k+1), k=0...nk-1. The scattering matrix [2, p.49, Eq.(2.135)]
	and the matrix of moments [similar bot not exactly equal to 2, p.91,
	Eq.(3.122)] are defined as:

               |a1 b1  0  0 |
               |b1 a2  0  0 |
           X = |0  0   a3 b2|                                               (1)
               |0  0  -b2 a4|

                | a1k -b1k  0    0  |
                |-b1k  a2k  0    0  |
           Xk = | 0    0    a3k -b2k|                                       (2)
                | 0    0    b2k  a4k|

	Note that Xk, Eq.(2), contains the expansion moments of Eq.(1) over matrix
	polynomials PI(mu) [2, p.91, Eq.(3.125)]. On input, Eq.(2) b1k is defined
	WITHOUT "-".

	Eqs.(2.152), (2.153), (2.156)-(2.159) [2, pp.53, 54] are used. Note that
	the [(l-2)!/(l+2)!]^1/2 is included in the Q polynomial. Using the
	following relation [2, Eq.(B.20), p.209]

           Pk20 = (-1)*Qk2                                                  (3)

	we write the following expressions that are used in this subroutine

           a1(x) =  sum(a1k*Qk0(x), k = 0..K)								(4)
           b1(x) = -sum(b1k*Qk2(x), k = 2..K)								(5)

REFERENCES:
	1. SORD_IP/SUMKA1B1.f90
	2. Hovenier JW, van der Mee C, Domke H., 2004: Transfer of polarized light in
       planetary atmospheres. Basic concepts and practical methods. Kluwer Publ.

PROTOTYPE:
	void sumka1b1(double *, int const, double *, double *, int const,
	              double *, double *)
------------------------------------------------------------------------------*/
//
//	LOCAL VARIABLES
	int ik, ix;
	double ki,
		   pa1, pa2, pa3,   // Legendre polynomials, Pk(x), for a1
		   pb1, pb2, pb3,   // Qk2(x) polynomials for b1
		   x2;              // x^2
//
//	LOCAL ARRAYS
	double *ca1, *ca2, *cb1, *cb2, *cb3; // coefficients for Pk(x) and Qk2(x)
//-----------------------------------------------------------------------------
//
	ca1 = new double[nk];
	ca2 = new double[nk];
	cb1 = new double[nk];
	cb2 = new double[nk];
	cb3 = new double[nk];
/*
//  Not used - skip initialization
	ca1[0] = 0.0; ca2[0] = 0.0; cb1[0] = 0.0; cb2[0] = 0.0; cb3[0] = 0.0;
	ca1[1] = 0.0; ca2[1] = 0.0; cb1[1] = 0.0; cb2[1] = 0.0; cb3[1] = 0.0;
	ca1[2] = 0.0; ca2[2] = 0.0; cb1[2] = 0.0; cb2[2] = 0.0; cb3[2] = 0.0;
*/
//  Precompute coefficients (note: sqrt)
	for (ik = 3; ik < nk; ik++) {
		ki = 1.0 / ik;
		ca1[ik] = 2.0 - ki;
		ca2[ik] = 1.0 - ki;
		cb1[ik] = 2.0 * ik - 1.0;
		cb2[ik] = sqrt((ik + 1.0) * (ik - 3.0));
		cb3[ik] = 1.0 / sqrt(1.0 * ik * ik - 4.0);
	} // for ik
//
	for (ix = 0; ix < nx; ix++) {
		x2 = x[ix] * x[ix];
//
//		Initialize a1 as sum(a1k*Pk0(x), k = 0, 1, 2):
		pa1 = x[ix];
		pa2 = 1.5 * x2 - 0.5;
		a1[ix] = a1k[0] + a1k[1] * x[ix] + a1k[2] * pa2; // a1k[0] = 1 or ssa
//
//		Initialize b1 as 0 + 0 + b1k*Q22(x)
		pb1 = 0.0;
		pb2 = C22 * (1.0 - x2);
		b1[ix] = -b1k[2] * pb2; // Pk20 = (-1)*Qk2
//		
//		Summation over k
		for (ik = 3; ik < nk; ik++) {
			pa3 = ca1[ik] * x[ix] * pa2 - ca2[ik] * pa1;
			pb3 = (cb1[ik] * x[ix] * pb2 - cb2[ik] * pb1) * cb3[ik];
			a1[ix] += a1k[ik] * pa3;
			b1[ix] -= b1k[ik] * pb3; // Pk20 = (-1)*Qk2
//
//			Redefine for recurrence
			pa1 = pa2;
			pa2 = pa3;
			pb1 = pb2;
			pb2 = pb3;
		} // for ik
	} // for ix
//
	delete[] ca1;
	delete[] ca2;
	delete[] cb1;
	delete[] cb2;
	delete[] cb3;
} // void sumka1b1(...)
/*------------------------------------------------------------------------------
 18/03/17 - delete -> delete[]: http://www.cplusplus.com/doc/tutorial/dynamic/
 18/03/16 - Tested for coarse aerosol from Kokhanovsky, JQSRT 2010 - ok:
		    da1 = [-0.02:0.02]%; db1 = [-0.025:0.01]% (round off, printing)
		    nmu = 1900; nk = 932.
 18/02/19 - First created and tested for Rayleigh - ok (ik = 0, 1, 2 only)
*/