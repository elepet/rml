#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

const int size = 4;

#include "rml.h"

int main(){
	float mat[size][size];
	rmlFill(0.0, mat);
	mat[1][1] = 3.0;
	mat[2][2] = 3.5;
	rmlPrint(mat);
	float vec[] = {2,3,4};
	rmlScaleMat(vec, mat);
	rmlPrint(mat);
	rmlRotateMat(vec, M_PI/2, mat);
	rmlPrint(mat);
	rmlTranslateMat(vec, mat);
	rmlPrint(mat);
	rmlProjectMat(800, 600, M_PI/2, 0.1, 40, mat);
	rmlPrint(mat);
	return 0;
}
