void polleg(int const k, double *x, int const nx, double *p);
/*------------------------------------------------------------------------------
PURPOSE:
	To compute the Legendre polynomials, P{k}(x)
    
INPUT:
	k    i[1]    Maximum order, k > 0
	x	 d[nx]   Abscissas, x = [-1:+1]
	nx   i[1]    Length of x

OUTPUT:
	p    d[(k+1)*nx]   Legendre polynomials
------------------------------------------------------------------------------*/