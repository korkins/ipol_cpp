#include <stdio.h> /* fprintf */
//
void print2d(double const **a, int const ny, int const nx, char *fn) {
/*------------------------------------------------------------------------------
PURPOSE:
	To print a 2d array of doubles into a txt file

INPUT:
	a    d[ny][nx]   array of doubles
	ny   i[1]        1st dimension of a
	nx   i[1]        2nd dimension of a
	fn   a[nf]       file name

OUTPUT:
	-

TREE:
	-

COMMENT:
	The txt file will be located in the solution folder.
	User defined path did not work for some reason...

REFERENCESS:
	-

PROTOTYPE:
	void print2d(double **, int const, int const, char *);
------------------------------------------------------------------------------*/
//
//	LOCAL VARIABLES
	FILE *pFile;
	int ix, iy;
//-----------------------------------------------------------------------------
//
	pFile = fopen(fn, "w");
	for (iy = 0; iy < ny; iy++) {
		for (ix = 0; ix < nx; ix++)
			if (ix < 3)
				fprintf(pFile, "%8.2f", a[iy][ix]); // angles
			else
				fprintf(pFile, "%12.8f", a[iy][ix]);
				//fprintf(pFile, "%16.8e", a[iy][ix]);
		fprintf(pFile, "\n");
	}
	fclose(pFile);
//
} // void print2d
/*------------------------------------------------------------------------------
 18/04/28 - double CONST **a is now used. 
 18/01/13 - First created and tested for a 2x3 array.
*/