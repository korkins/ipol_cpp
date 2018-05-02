void polrtm(int const m, int const k, double *x, int const nx,
														 double *r, double *t);
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
------------------------------------------------------------------------------*/