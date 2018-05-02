void zrssgm(int const m, double const def, double *mu0, int const nsz,
			double *mug, int const ng2, double *zr1, double *zr2, double *zr3);
/*------------------------------------------------------------------------------
PURPOSE:
	To compute the 1st column of the Rayleigh phase matrix for the Fourier order
	m > 0 at the Gauss nodes

INPUT:
	m     i(1)     The Fourier moment, m = 1 or 2
	def   d(1)     Depolarization factor, def = [0 ... ~6/7]
	mu0   d(nsz)   cos(sza) > 0; avoid mu0 == 1.0 (replace with 0.999...999)
	nsz   i(1)     Length of mu0
	mug   d(ng2)   Gauss nodes
	ng2   i(1)     Length of mug

OUTPUT:
	zr1, zr2, zr3   d(ng2*nsz)   Elements of the phase matrix; ig=1:ng2 runs 1st
------------------------------------------------------------------------------*/