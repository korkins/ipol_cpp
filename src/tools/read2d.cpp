#include <stdio.h> /* fprintf */
//
void read2d(char *fn, int const ny, int const nx, double **a) {
/*------------------------------------------------------------------------------
PURPOSE:
	To read a 2d array of doubles from a txt file

INPUT:
	fn   a[x]   file name
	ny   i[1]   1st dimension of a
	nx   i[1]   2nd dimension of a

OUTPUT:
	a    d[ny][nx]   array of doubles

TREE:
	-

COMMENT:
	The txt file will be located in the solution folder.
	User defined path did not work for some reason...

REFERENCESS:
	-

PROTOTYPE:
	void read2d(char *, int const, int const, double **);
------------------------------------------------------------------------------*/
//
//	LOCAL VARIABLES
	FILE *pFile;
	int ix, iy;
//-----------------------------------------------------------------------------
//
	pFile = fopen(fn, "r");
	for (iy = 0; iy < ny; iy++)
		for (ix = 0; ix < nx; ix++)
			fscanf(pFile, "%lf", &a[iy][ix]); // %lf if a[][] is DOUBLE
	fclose(pFile);
//
} // void read2d
/*------------------------------------------------------------------------------
 16Mar18 - First created and tested for a 4x3 array in ain.txt:
				 1.0 -2.0e+1  3.0
				 4.0 -4.0  5.0
				-0.1  0.2e-1  0.3e+3
				 0.0 -1.0 -100.0e-1

		   print2d gives in aOUT.txt - reads ok
				 1.00  -20.00    3.00
				 4.00   -4.00    5.00
				 -0.10    0.02  300.00
				 0.00   -1.00  -10.00
*/