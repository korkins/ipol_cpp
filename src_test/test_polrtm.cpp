#include <stdio.h>				/* printf */
#include <time.h>
#include <iostream>
#include <math.h>
//using namespace std;

void polrtm(int const, int const, double *, int const, double *, double *);
double factorial(int);

void main()
{
	double *x, *rkm, *tkm, *rkk, *tkk, *xp, *xn, *rp, *rn, *tp, *tn;
	char cvar;
	int ix, ik, k, k1, nx, kk, kk1, im;
	double dx, a, b, c, d, ck, fk, x2, x3, x4, x5, x6, x7, drkk, dtkk, 
		   dr, dt, r, t, sx, p1, p2;
	
	FILE *fp;
//
	nx = 21;
	k = 10;
	k1 = k+1;
	dx = 2.0/(nx-1);
	x  = new double [nx];
	rkm = new double [(k+1)*nx];
	tkm = new double [(k+1)*nx];
//
	x[0] = -1.0;
	for (ix = 1; ix < nx; ix++)
		x[ix] = x[ix-1] + dx;
    x[nx-1] = 1.0;
//  to avoid sx = sqrt(1 - x2) = 0:
	x[0] = -1.0 + dx/2.0;
	x[nx-1] = 1.0 - dx/2.0;
//
	double start = (double)clock() /(double) CLOCKS_PER_SEC;
	polrtm(1, k, x, nx, rkm, tkm);
	double end = (double)clock() / (double) CLOCKS_PER_SEC;

	for (ix = 0; ix < nx; ix++)
		for (ik = 0; ik < k1; ik++)
			printf("%6.2f %6i %18.6e %18.6e \n", x[ix], ik, rkm[ix*k1+ik], tkm[ix*k1+ik]);

	printf("\n");
	printf("cpu time %fs\n", (end - start));

	delete(rkm);
	delete(tkm);

//-------------------------------------------------------------------------
//
	printf("\n");
	printf("\n");
	printf("Test vs Rkk & Tkk \n");
	printf("\n");
	for (kk = 2; kk < 21; kk++) {
		a = 1.0/pow(2.0, kk);
		b = factorial(2*kk);
		c = factorial(kk-2);
		d = factorial(kk+2);

		ck = a*sqrt(b/c/d);

		rkk = new double [nx];
		tkk = new double [nx];
		
		for (ix = 0; ix < nx; ix++) {
			x2 = x[ix]*x[ix];
			fk = pow( (1.0 - x2), (kk-2.0)/2.0 );
			rkk[ix] = ck*fk*(1.0 + x2);
			tkk[ix] = 2.0*ck*fk*x[ix];
		}
	//
		kk1 = kk+1;
		rkm = new double [(kk1)*nx];
		tkm = new double [(kk1)*nx];
		polrtm(kk, kk, x, nx, rkm, tkm);
		drkk = 0.0;
		dtkk = 0.0;
		for (ix = 0; ix < nx; ix++) {
		/*
			printf("%6.2f %6i %10.1e %10.1e %10.1e %10.1e %10.1e %10.1e \n", x[ix], kk, 
				rkm[ix*kk1+kk], rkk[ix], rkm[ix*kk1+kk] - rkk[ix], 
					tkm[ix*kk1+kk], tkk[ix], tkm[ix*kk1+kk] - tkk[ix]); */
			drkk += abs( rkm[ix*kk1+kk] - rkk[ix] );
			dtkk += abs( tkm[ix*kk1+kk] - tkk[ix] );
		}
//
		printf("%6i %10.1e %10.1e \n", kk, drkk/nx, dtkk/nx);
//		
		delete(rkk);
		delete(tkk);
		delete(rkm);
		delete(tkm);
	} // for (kk = 3; kk < 5; kk++)
//
//-------------------------------------------------------------------------
	printf("\n");
	printf("\n");
	printf("theory: m=1, k=6 \n");
	kk1 = 6+1;
	rkm = new double [(kk1)*nx];
	tkm = new double [(kk1)*nx];
	polrtm(1, 6, x, nx, rkm, tkm);
/*
	ik = 4;
	for (ix = 0; ix < nx; ix++)
		printf("%6.2f %6i %18.6e %18.6e \n", x[ix], ik, rkm[ix*7+ik], tkm[ix*7+ik]);
*/	

	dr = 0.0;
	dt = 0.0;

	c = sqrt(10.0)/32.0;
	for (ix = 0; ix < nx; ix++) {
	    x2 = x[ix]*x[ix];
		sx = sqrt(1.0 - x2);
		x3 = x[ix]*x2;
		x4 = x[ix]*x3;
		x5 = x[ix]*x4;
		x6 = x[ix]*x5;
		x7 = x[ix]*x6;
		r = c*(99.0*x7 - 201.0*x5 + 121.0*x3 - 19.0*x[ix])/sx;
		t = c*(33.0*x6 -  51.0*x4 +  19.0*x2 -  1.0)/sx;
		dr += abs( rkm[ix*kk1+6] - r );
		dt += abs( tkm[ix*kk1+6] - t );
	}

	printf("%6i %10.1e %10.1e \n", 6, dr/nx, dt/nx);

	delete(rkm);
	delete(tkm);

//-------------------------------------------------------------------------
	printf("\n");
	printf("\n");
	printf("theory: m=2, k=6 \n");
	kk1 = 6+1;
	rkm = new double [(kk1)*nx];
	tkm = new double [(kk1)*nx];
	polrtm(2, 6, x, nx, rkm, tkm);

/*
	ik = 2;
	for (ix = 0; ix < nx; ix++)
		printf("%6.2f %6i %18.6e %18.6e \n", x[ix], ik, rkm[ix*7+ik], tkm[ix*7+ik]);
*/	

	dr = 0.0;
	dt = 0.0;

	c = 1.0/64.0;
	for (ix = 0; ix < nx; ix++) {
		x2 = x[ix]*x[ix];
		x3 = x[ix]*x2;
		x4 = x[ix]*x3;
		x5 = x[ix]*x4;
		x6 = x[ix]*x5;
		x7 = x[ix]*x6;

		p1 = c*(1.0 + x[ix])*(1.0 + x[ix])*(495.0*x4 - 660.0*x3 + 90.0*x2 + 108.0*x[ix] - 17.0);
		p2 = c*(1.0 - x[ix])*(1.0 - x[ix])*(495.0*x4 + 660.0*x3 + 90.0*x2 - 108.0*x[ix] - 17.0);

		r = 0.5*(p1 + p2);
		t = 0.5*(p1 - p2);
		dr += abs( rkm[ix*kk1+6] - r );
		dt += abs( tkm[ix*kk1+6] - t );
	}

	printf("%6i %10.1e %10.1e \n", 6, dr/nx, dt/nx);

	delete(rkm);
	delete(tkm);

//-------------------------------------------------------------------------

	printf("\n");
	printf("\n");
	printf("theory: m=3, k=6 \n");
	kk1 = 6+1;
	rkm = new double [(kk1)*nx];
	tkm = new double [(kk1)*nx];
	polrtm(3, 6, x, nx, rkm, tkm);

/*
	ik = 2;
	for (ix = 0; ix < nx; ix++)
		printf("%6.2f %6i %18.6e %18.6e \n", x[ix], ik, rkm[ix*7+ik], tkm[ix*7+ik]);
*/	

	dr = 0.0;
	dt = 0.0;

	c = 3.0/32.0;
	for (ix = 0; ix < nx; ix++) { 
		x2 = x[ix]*x[ix];
		sx = sqrt( (1.0 - x[ix])/(1.0 + x[ix]) );
		x3 = x[ix]*x2;
		x4 = x[ix]*x3;
		x5 = x[ix]*x4;
		x6 = x[ix]*x5;
		x7 = x[ix]*x6;

		p1 = c*sx*(55.0*x6 + 110.0*x5 + 5.0*x4 - 92.0*x3 - 31.0*x2 + 14.0*x[ix] + 3.0);
		p2 = -(c/sx)*(55.0*x6 - 110.0*x5 + 5.0*x4 + 92.0*x3 - 31.0*x2 - 14.0*x[ix] + 3.0);

		r = 0.5*(p1 + p2);
		t = 0.5*(p1 - p2);
		dr += abs( rkm[ix*kk1+6] - r );
		dt += abs( tkm[ix*kk1+6] - t );
	}

	printf("%6i %10.1e %10.1e \n", 6, dr/nx, dt/nx);

	delete(rkm);
	delete(tkm);

//-------------------------------------------------------------------------

	printf("\n");
	printf("\n");
	printf("theory: m=4, k=6 \n");
	kk1 = 6+1;
	rkm = new double [(kk1)*nx];
	tkm = new double [(kk1)*nx];
	polrtm(4, 6, x, nx, rkm, tkm);

/*
	ik = 2;
	for (ix = 0; ix < nx; ix++)
		printf("%6.2f %6i %18.6e %18.6e \n", x[ix], ik, rkm[ix*7+ik], tkm[ix*7+ik]);
*/	

	dr = 0.0;
	dt = 0.0;

	c = -sqrt(30.0)/64.0;
	for (ix = 0; ix < nx; ix++) { 
		x2 = x[ix]*x[ix];
		sx = sqrt( (1.0 - x[ix])/(1.0 + x[ix]) );
		x3 = x[ix]*x2;
		x4 = x[ix]*x3;
		x5 = x[ix]*x4;
		x6 = x[ix]*x5;
		x7 = x[ix]*x6;

		p1 = c*(x[ix] - 1.0)*(33.0*x5 + 77.0*x4 + 34.0*x3 - 30.0*x2 - 19.0*x[ix] + 1.0);
		p2 = c*(x[ix] + 1.0)*(33.0*x5 - 77.0*x4 + 34.0*x3 + 30.0*x2 - 19.0*x[ix] - 1.0);

		r = 0.5*(p1 + p2);
		t = 0.5*(p1 - p2);
		dr += abs( rkm[ix*kk1+6] - r );
		dt += abs( tkm[ix*kk1+6] - t );
	}

	printf("%6i %10.1e %10.1e \n", 6, dr/nx, dt/nx);

	delete(rkm);
	delete(tkm);

//-------------------------------------------------------------------------

	printf("\n");
	printf("\n");
	printf("theory: m=128, k=0:1000, print out k = 255 \n");
	kk1 = 1000;
	rkm = new double [(kk1)*nx];
	tkm = new double [(kk1)*nx];
	polrtm(128, kk1-1, x, nx, rkm, tkm);

	ik = 255;
	for (ix = 0; ix < nx; ix++)
		printf("%6.2f %6i %18.6e %18.6e \n", x[ix], ik, rkm[ix*kk1+ik], tkm[ix*kk1+ik]);

	fp = fopen("RTcpp.txt", "w");
	ik = 255;
	for (ix = 0; ix < nx; ix++)
		fprintf(fp, "%6.2f %6i %24.16e %24.16e \n", x[ix], ik, rkm[ix*kk1+ik], tkm[ix*kk1+ik]);
	
	fclose(fp);
	delete(rkm);
	delete(tkm);

//-------------------------------------------------------------------------
	printf("\n");
	printf("\n");
	printf("SYMMETRY \n");
	nx = 1000;
	k = 1000;
	kk1 = k+1;
	xp = new double [nx];
	rp = new double [(kk1)*nx];
	tp = new double [(kk1)*nx];
	xn = new double [nx];
	rn = new double [(kk1)*nx];
	tn = new double [(kk1)*nx];

	xp[0] = 1.0/nx;
	xn[0] = -xp[0];
	for (ix = 1; ix < nx; ix++) { 
		xp[ix] = xp[ix-1] + 1.0/nx;
		xn[ix] = -xp[ix];
	}
	dr = xp[nx-1];
	dt = xn[nx-1];
	xp[nx-1] =  0.999999;
    xn[nx-1] = -0.999999;

	dr = 0.0;
	dt = 0.0;
	for (im = 1; im < k; im++) {
		printf("%i \n", im);
		polrtm(im, k, xp, nx, rp, tp);
		polrtm(im, k, xn, nx, rn, tn);
		for (ix = 0; ix < nx; ix++) {
			for (ik = im; ik < kk1; ik++){
				dr = dr + abs( rp[ix*kk1+ik] - pow( -1.0, (im+ik) )*rn[ix*kk1+ik] );
				//printf("%6.2f %6i %16.6e %16.6e \n", xp[ix], ik, rp[ix*kk1+ik],
											//pow( -1.0, (im+ik) )*rn[ix*kk1+ik]);
				dt = dt + abs( tp[ix*kk1+ik] + pow( -1.0, (im+ik) )*tn[ix*kk1+ik] );
			}
		}
	}
	//dr = dr; ///(kk1*k*nx);

	printf("%18.6e %18.6e \n", dr, dt);

	printf("Done! \n");
	std::cin >> cvar;
//
}



double factorial(int n){
	double f;
	f = 1.0;
	for (int i = 1; i <= n; i++)
		f *= i; 
	return f;
}