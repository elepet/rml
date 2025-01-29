#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "rml.h"

int main(){
	//Usage examples.
	rmlVector* vec = rmlAllocateVector(3);
	vec->val[0] = 1.0;
	rmlMatrix* mat = rmlRotationMatrix(vec, M_PI / 2);
	rmlPrintMatrix(mat);
	rmlVector* vec3 = rmlAllocateVector(4);
	vec3->val[0] = 1.0;
	vec3->val[1] = 2.0;
	vec3->val[2] = 3.0;
	vec3->val[3] = 4.0;
	rmlVector* vec2 = rmlDotMatrixVector(mat, vec3);
	rmlPrintVector(vec2);
	rmlFreeVector(vec2);
	rmlFreeMatrix(mat);
	rmlFreeVector(vec);
	return 0;
}
