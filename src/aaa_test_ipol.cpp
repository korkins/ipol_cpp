# include <math.h>
# include <stdio.h>
//
void main() {
/*------------------------------------------------------------------------------
PURPOSE:
	To test ipol using published benchmarks - for the list of tests see SORD

COMMENTS:
	1) new/delete vs malloc/free [1]:

 Feature                  | new/delete                     | malloc/free                   
--------------------------+--------------------------------+-------------------------------
 Memory allocated from    | 'Free Store'                   | 'Heap'                        
 Returns                  | Fully typed pointer            | void*                         
 On failure               | Throws (never returns NULL)    | Returns NULL                  
 Required size            | Calculated by compiler         | Must be specified in bytes    
 Handling arrays          | Has an explicit version        | Requires manual calculations  
 Reallocating             | Not handled intuitively        | Simple (no copy constructor)  
 Call of reverse          | Implementation defined         | No                            
 Low memory cases         | Can add a new memory allocator | Not handled by user code      
 Overridable              | Yes                            | No                            
 Use of (con-)/destructor | Yes                            | No  

	2) Example with arrays: M[slow_index][fast_index]

		M[0][0] = 1;
		M[0][1] = 2;
		M[0][2] = 3;

		M[1][0] = 4;
		M[1][1] = 5;
		M[1][2] = 6;

		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 3; j++) printf("M[%i][%i]=%i ", i, j, M[i][j]);
			printf("\n");
		}

REFERENCES:
	1. http://stackoverflow.com/questions/240212/what-is-the-difference-between-new-delete-and-malloc-free
------------------------------------------------------------------------------*/
//
// PARAMETERS
//
//	Max sizes of arrays
	const int
		NAZ_MAX = 181,       // Azimuth: view
		NGA_MAX = 360,       // Azimuth: Gauss nodes
		NVZ_MAX = 180,       // Zenith: view (NUD - user defined - in SORD)
		NSZ_MAX = 90,        // Number of solar zeniths
		NG1_MAX = 1024,      // Zenith: Gauss nodes per hemisphere
		NG2_MAX = 2*NG1_MAX, // Zenith: Gauss nodes per sphere
		NK_MAX  = NG2_MAX,   // Moments of the phase matrix
		NLR_MAX = 100,       // Number of optical layers
		NM_MAX  = NG1_MAX,   // Fourier moments
		NTST_MAX = 100;      // Number of tests
	const double
		PI = 3.1415926535897932,
		PI2 = 2.0*PI,
		D2R = PI/180.0,
		R2D = 180.0/PI,
		TINY = 1.0e-12;
//
//	LOCAL VARIABLES
//
	int
		itest, irun, ia, iz, isrf, nks, ng1, ng2, nga, nup, ndn, nmu, nl, nlr,
		nml, nm, naz, nvz, nsz, nk, nx, ny, ibr, nz, nkm, ivz, isz, iaz, nup, ndn;
	double
		cpu_t1, cpu_t2, time, depf, tau0, epsi;
//
//	LOCAL ARRAYS
//
	int
		tstnum[NTST_MAX];
	double
		sza[NSZ_MAX], mu0[NSZ_MAX], vza[NVZ_MAX], mu[NVZ_MAX], aza[NVZ_MAX],
		azi[NVZ_MAX], tau[NLR_MAX], ssa[NLR_MAX];

	int
		M[2][3];
		
		//*c1, *c2, *c3; // c1[k+1], c2[k+1], c3[k+1] - coefficients
//-----------------------------------------------------------------------------
//
	printf("Try formats: %s %i\n %8.2f\n %12.3e\n", "NG2_MAX = ", NG2_MAX, 1.0*NG2_MAX, 1.0*NG2_MAX);
//
//  Bu default: all tests ON
    for (itest = 0; itest < NTST_MAX; itest++) tstnum[itest] = 1;
//
	/*
	I.
	for (itest = 0; itest < NTST_MAX; itest++)
		printf("itest = %i; tstnum[itst] = %i\n", itest, tstnum[itest]);
	II.
		M[0][0] = 1;
		M[0][1] = 2;
		M[0][2] = 3;

		M[1][0] = 4;
		M[1][1] = 5;
		M[1][2] = 6;

		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 3; j++) printf("M[%i][%i]=%i ", i, j, M[i][j]);
			printf("\n");
		}
	*/
//
	/*
	************** OPEN FILES HERE *****************
	*/
//	Loop over runs
for (irun = 0; irun < 1; irun++)
{
/*
*******************************************************************************
	001. GOAL: to test single scattering approximation in a homogeneous layer
			   with aerosol scattering over black surface
*******************************************************************************
*/
	itest = 0;
	if (tstnum[itest] == 1)
	{
// 1		
//		ACCURACY
// 1
//		Number of Gauss nodes per sphere (N/U)
		ng1 = 2;
//		Total number of the Fourier moments (N/U)
		nm = 0;
//		Convergence accuracy parameter
		epsi = -999.0;
// 1
//		GEOMETRY
// 1
//		Solar zenith angle, SZA, in degrees or mu0=cos(SZA)
		mu0[0] = 0.6;
		sza[0] = acos(mu0[0])*R2D;
//		View zenith angles, VZA, in degrees or mu = -cos(VZA)
		nvz = 2;
		nup = 1;
		mu[0] = -mu0[0]; // up
		ndn = nvz-nup;
		mu[1] =  mu0[0]; // down
		for (ivz = 0; ivz < nvz; ivz++) vza[ivz] = 180.0 - acos(mu[ivz])*R2D;
	} // if tstnum[0] == 1
/*
*******************************************************************************
	002. GOAL: to test double scattering approximation in a homogeneous layer
			   with aerosol scattering over black surface
*******************************************************************************
*/
	itest += 1;
	if (tstnum[itest] == 1)
	{

	} // if tstnum[0] == 1
	
} // for irun
//
} // void main