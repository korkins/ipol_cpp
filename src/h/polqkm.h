void polqkm(int const m, int const k, double *x, int const nx, double *q);
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
------------------------------------------------------------------------------*/