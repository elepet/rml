#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include "rml.h"

int main(){
	//Usage examples.
	rmlVector* vec = rmlAllocateVector(3);
	rmlVector* vec2 = rmlAllocateVector(3);
	vec->val[0] = 1.0;
	vec->val[1] = 1.0;
	vec->val[2] = 0,0;
	vec2->val[0] = 1.3;
	vec2->val[1] = 0.8;
	vec2->val[2] = 3.9;
	rmlVector* vec3 = rmlAddVector(vec, 0, vec2);
	rmlPrintVector(vec3);
	printf("%.2f\n",rmlLengthOfVector(vec));
	rmlFreeVector(vec);
	rmlFreeVector(vec2);
	rmlFreeVector(vec3);
	return 0;
}
