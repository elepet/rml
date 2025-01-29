#include <stdio.h>
#include <stdlib.h>

#include "rml.h"

int main(){
	rmlVector* vec = rmlAllocateVector(3);
	vec->val[0] = 2.0;
	vec->val[1] = 3.0;
	vec->val[2] = 4.0;
	rmlMatrix* mat = rmlScalingMatrix(vec);
	rmlPrintMatrix(mat);
	rmlFreeMatrix(mat);
	rmlFreeVector(vec);

	return 0;
}
