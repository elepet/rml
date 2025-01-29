//Renderer Matrix Library
//Operations modify, calls for specific mat/vec return pointers.
//TODO: maybe add garbage collector?

#include <stdio.h>
#include <stdlib.h>

#include "rml.h"

//Returns pointer to row*col matrix initialised with zeros.
//Dynamically allocates. User must free with rmlFreeMatrix.
rmlMatrix* rmlAllocateMatrix(unsigned int row, unsigned int col) {
	if (row == 0 || col == 0) {
		printf("RML_ERROR_ALLOCATEMATRIX: row == 0 || col == 0.\n");
		return NULL;
	}
	rmlMatrix* mat = (rmlMatrix*)malloc(sizeof(rmlMatrix));
	if (mat == NULL) {
		printf("RML_ERROR_ALLOCATEMATRIX: allocation failed.\n");
		return NULL;
	}
	float** val = (float**)calloc(row, sizeof(float*));
	if (val == NULL) {
		free(mat);
		printf("RML_ERROR_ALLOCATEMATRIX: allocation failed.\n");
		return NULL;
	}
	for (int i = 0; i < row; i++) {
		val[i] = (float*)calloc(col, sizeof(float));
		if (val[i] == NULL) {
			for (int j = 0; j < i; j++) {
				free(val[j]);
			}
			free(val);
			free(mat);
			printf("RML_ERROR_ALLOCATEMATRIX: allocation failed.\n");
			return NULL;
		}
	}
	mat->row = row;
	mat->col = col;
	mat->val = val;
	return mat;
}

void rmlFreeMatrix(rmlMatrix* mat) {
	if (mat != NULL) {
		for (int i = 0; i < mat->row; i++) {
			free(mat->val[i]);
		}
		free(mat->val);
		free(mat);
	} else printf("RML_ERROR_FREEMATRIX: passed NULL.\n");
}

void rmlPrintMatrix(rmlMatrix* mat) {
	if (mat != NULL) {
		for (int i = 0; i < mat->row; i++) {
			for (int j = 0; j < mat->col; j++) {
				printf("%.2f ", mat->val[i][j]);
			}
		printf("\n");
		}
	} else printf("RML_ERROR_PRINTMATRIX: passed NULL.\n");
}

rmlVector* rmlAllocateVector(unsigned int size) {
	if (size == 0) {
		printf("RML_ERROR_ALLOCATEVECTOR: size == 0.\n");
		return NULL;
	}
	rmlVector* vec = (rmlVector*)malloc(sizeof(rmlVector));
	if (vec == NULL) {
		printf("RML_ERROR_ALLOCATEVECTOR: allocation failed.\n");
		return NULL;
	}
	float* val = (float*)calloc(size, sizeof(float));
	if (val == NULL) {
		free(vec);
		printf("RML_ERROR_ALLOCATEVECTOR: allocation failed.\n");
		return NULL;
	}
	vec->size = size;
	vec->val = val;
	return vec;
}

void rmlFreeVector(rmlVector* vec) {
	if (vec != NULL) {
		free(vec->val);
		free(vec);
	} else printf("RML_ERROR_FREEVECTOR: passed NULL.\n");
}

void rmlPrintVector(rmlVector* vec) {
	if (vec != NULL) {
		for (int i = 0; i < vec->size; i++) {
			printf("%.2f\n", vec->val[i]);
		}
	} else printf("RML_ERROR_PRINTVECTOR: passed NULL.\n");
}

//Modifies matrix pointed to by 'out' to dot product of two matrices (or vectors).
//Commutative for vectors but not for matrices.
void rmlDotProduct(rmlMatrix* in1, rmlMatrix* in2, rmlMatrix* out) {
	if (in1->col != in2->row) {
		printf("RML_ERROR_DOTPRODUCT: in1->col != in2->row.\n");
		return;
	}
	if (out->row != in1->row || out->col != in2->col) {
		printf("RML_ERROR_DOTPRODUCT: out->row != in1->row || out->col != in2->col.\n");
		return;
	}
	for (int i = 0; i < in1->row; i++) {
		for (int j = 0; j < in2->col; j++) {
			for (int k = 0; k < in1->col; k++) {
				out->val[i][j] += in1->val[i][k] * in2->val[k][j];
			}
		}
	}
}; 

void rmlScalarProduct(float scalar, rmlMatrix* mat) {
	if (mat != NULL) {
		for (int i = 0; i < mat->row; i++) {
			for (int j = 0; j < mat->col; j++) {
				mat->val[i][j] = scalar * mat->val[i][j];
			}
		}
	} else printf("RML_ERROR_SCALARPRODUCT: passed NULL.\n");
}

//TODO: Fix!!
//Modifies matrix pointed to by 'mat' from row-major to column-major or vice versa.
void rmlSwitchMajor(rmlMatrix* mat) {
}

//Returns pointer to size*size identity matrix.
//Dynamically allocates. User must free with rmlFreeMatrix.
rmlMatrix* rmlIdentityMatrix(unsigned int size) {
	rmlMatrix* mat = rmlAllocateMatrix(size, size);
	for (int i = 0; i < mat->row; i++) {
		for (int j = 0; j < mat->col; j++) {
			if (i == j) {
				mat->val[i][j] = 1.0;
			} else mat->val[i][j] = 0.0;
		}
	}
	return mat;
}

//Returns pointer to 3*3 or 4*4 scaling matrix specified by 'vec'.
//Dynamically allocates. User must free with rmlFreeMatrix.
rmlMatrix* rmlScalingMatrix(rmlVector* vec) {	
	if (!(vec->size == 2 || vec->size  == 3)) {
		printf("RML_ERROR_SCALINGMATRIX: !(vec->size == 2 || vec->size  == 3).\n");
		return NULL;
	}
	rmlMatrix* mat = rmlAllocateMatrix(vec->size + 1, vec->size + 1);
	for (int i = 0; i < mat->row; i++) {
		for (int j = 0; j < mat->col; j++) {
			if (i == j && i != vec->size) {
				mat->val[i][j] = vec->val[i];
			} else if (i ==j && i == vec->size) mat->val[i][j] = 1.0;
		}
	}
	return mat;
}

//Returns pointer to 3*3 or 4*4 translation matrix specified by 'vec'.
//Dynamically allocates. User must free with rmlFreeMatrix.
rmlMatrix* rmlTranslationMatrix(rmlVector* vec) {
	if (!(vec->size == 2 || vec->size  == 3)) {
		printf("RML_ERROR_TRANSLATIONMATRIX: !(vec->size == 2 || vec->size  == 3).\n");
		return NULL;
	}
	rmlMatrix* mat = rmlIdentityMatrix(vec->size + 1);
	mat->val[0][mat->col - 1] = vec->val[0];
	mat->val[1][mat->col - 1] = vec->val[1];
	if (vec->size == 3) mat->val[2][mat->col - 1] = vec->val[2];
	return mat;
}

//Local space -> world space.
//Import model in local space coords.
void rmlModelMatrix(float **mat, unsigned int row, unsigned int col, float x, float y, float z) {
	
}

//World space -> view space.
void rmlViewMatrix(float **mat, unsigned int row, unsigned int col, float x, float y, float z) {

}

//View space -> clip space.
//Run viewport transform after to get screen space.
void rmlProjectionMatrix(unsigned int row, unsigned int col) {
	//Apply matrix or return matrix?
}

