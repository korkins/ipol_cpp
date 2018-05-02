#include <stdio.h> /* fprintf */
//
void print1d(double *x, int const nx, char *fn) {
/*------------------------------------------------------------------------------
PURPOSE:
	To print a 1d array of doubles into a txt file

INPUT:
	x    d[nx]   array of doubles
	nx   i[1]    2nd dimension of a
	fn   a[nf]   file name

OUTPUT:
	-

TREE:
	-

COMMENT:
	The txt file will be located in the solution folder.
	User defined path did not work for some reason...
	The input vector, x, is printed as a column.

REFERENCESS:
	-

PROTOTYPE:
	void print1d(double *, int const, char *);
------------------------------------------------------------------------------*/
//
//	LOCAL VARIABLES
	FILE *pFile;
	int ix;
//-----------------------------------------------------------------------------
//
	pFile = fopen(fn, "w");
	for (ix = 0; ix < nx; ix++)
		fprintf(pFile, "%12.8f \n", x[ix]);
	fclose(pFile);
//
} // void print1d
/*------------------------------------------------------------------------------
 18/01/14 - First created and tested for a vector
*/