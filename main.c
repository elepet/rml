#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include "rml.h"

int main(){
	rmlMatrix* mat = rmlAllocateMatrix(3,4);
	rmlFillMatrix(1.0, mat);
	rmlPrintMatrix(mat);
	rmlFreeMatrix(mat);
	return 0;
}
