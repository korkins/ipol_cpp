void zrssg0(double const def, double *mu0, int const nsz, double *mug,
			int const ng2,  double *zr1, double *zr2);
/*------------------------------------------------------------------------------
PURPOSE:
	To compute the 1st column of the Rayleigh phase matrix for the Fourier order
	m = 0 at Gauss nodes

INPUT:
	def   d(1)     Depolarization factor, def = [0 ... ~6/7]
	mu0   d(nsz)   cos(sza) > 0
	nsz   i(1)     Length of mu0
	mug   d(ng2)   Gauss nodes
	ng2   i(1)     Length of mug

OUTPUT:
	zr1, zr2   d(ng2*nsz)   Elements of the phase matrix, ig = 1:ng2 changes 1st
------------------------------------------------------------------------------*/