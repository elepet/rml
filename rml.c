//Renderer Matrix Library
//TODO: garbage collector, quaternions, transpose/major, camera

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
		printf("\n");
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
		printf("\n");
	} else printf("RML_ERROR_PRINTVECTOR: passed NULL.\n");
}

//Returns pointer to matrix that is dot product of input matrices.
//If input is m*n.n*p, output is m*p.
//Not commutative.
//Dynamically allocates. User must free with rmlFreeMatrix.
rmlMatrix* rmlDotMatrix(rmlMatrix* in1, rmlMatrix* in2) {
	if (in1->col != in2->row) {
		printf("RML_ERROR_DOTMATRIX: in1->col != in2->row.\n");
		return NULL;
	}
	rmlMatrix* out = rmlAllocateMatrix(in1->row, in2->col);
	for (int i = 0; i < in1->row; i++) {
		for (int j = 0; j < in2->col; j++) {
			for (int k = 0; k < in1->col; k++) {
				out->val[i][j] += in1->val[i][k] * in2->val[k][j];
			}
		}
	}
	return out;
}; 

//Returns dot product of two vectors.
//Commutative.
float rmlDotVector(rmlVector* in1, rmlVector* in2) {
	if (in1->size != in2->size) {
		printf("RML_ERROR_DOTVECTOR: in1->size != in2->size.\n");
		return NAN;
	}
	float out = 0.0;
	for (int i = 0; i < in1->size; i++) {
		out += in1->val[i] * in2->val[i];
	}
	return out;
}; 

//Returns pointer to vector that is dot product of input matrix and vector.
//Dynamically allocates. User must free with rmlFreeMatrix.
rmlVector* rmlDotMatrixVector(rmlMatrix* mat, rmlVector* vec) {
	if (mat->col != vec->size) {
		printf("RML_ERROR_DOTMATRIXVECTOR: mat->col != vec->size.\n");
		return NULL;
	}
	rmlVector* out = rmlAllocateVector(mat->row);
	for (int i = 0; i < mat->row; i++) {
		for (int j = 0; j < mat->col; j++) {
			out->val[i] += mat->val[i][j] * vec->val[j];
		}
	}
	return out;
}; 

//Multiplies every member of matrix by scalar.
//Uniform scale.
void rmlScalarByMatrix(float scalar, rmlMatrix* mat) {
	if (mat != NULL) {
		for (int i = 0; i < mat->row; i++) {
			for (int j = 0; j < mat->col; j++) {
				mat->val[i][j] = scalar * mat->val[i][j];
			}
		}
	} else printf("RML_ERROR_SCALARBYMATRIX: passed NULL.\n");
}

//Multiplies every member of vector by scalar.
//Uniform scale.
void rmlScalarByVector(float scalar, rmlVector* vec) {
	if (vec != NULL) {
		for (int i = 0; i < vec->size; i++) {
			vec->val[i] = scalar * vec->val[i];
		}
	} else printf("RML_ERROR_SCALARBYVECTOR: passed NULL.\n");
}

//Returns pointer to vector that is cross product of input vectors.
//Only defined for size = 3.
//Anticommutative, so a x b = - b x a.
//Dynamically allocates. User must free with rmlFreeMatrix.
rmlVector* rmlCrossVector(rmlVector* in1, rmlVector* in2) {
	if (in1->size != in2->size || in1->size != 3) {
		printf("RML_ERROR_CROSSVECTOR: in1->size != in2->size || in1->size != 3.\n");
		return NULL;
	}
	rmlVector* out = rmlAllocateVector(3);
	out->val[0] = in1->val[1] * in2->val[2] - in1->val[2] * in2->val[1];		
	out->val[1] = in1->val[2] * in2->val[0] - in1->val[0] * in2->val[2];		
	out->val[2] = in1->val[0] * in2->val[1] - in1->val[1] * in2->val[0];		
	return out;
}

void rmlNormaliseVector(rmlVector* vec) {
	if (vec != NULL) {
		float mag = sqrt(vec->val[0]*vec->val[0]+vec->val[1]*vec->val[1]+vec->val[2]*vec->val[2]);
		if (mag == 0.0) {
			printf("RML_ERROR_NORMALVECTOR: mag == 0.0.\n");
			return;
		}
		vec->val[0] = vec->val[0] / mag;
		vec->val[1] = vec->val[1] / mag;
		vec->val[2] = vec->val[2] / mag;
	} else printf("RML_ERROR_NORMALVECTOR: passed NULL.\n");
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
//Not uniform scale.
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
			} else if (i == j && i == vec->size) mat->val[i][j] = 1.0;
		}
	}
	return mat;
}

//Returns pointer to 3*3 or 4*4 translation matrix specified by 'vec'.
//Dynamically allocates. User must free with rmlFreeMatrix.
rmlMatrix* rmlTranslationMatrix(rmlVector* vec) {
	if (!(vec->size == 2 || vec->size == 3)) {
		printf("RML_ERROR_TRANSLATIONMATRIX: !(vec->size == 2 || vec->size == 3).\n");
		return NULL;
	}
	rmlMatrix* mat = rmlIdentityMatrix(vec->size + 1);
	mat->val[0][mat->col - 1] = vec->val[0];
	mat->val[1][mat->col - 1] = vec->val[1];
	if (vec->size == 3) mat->val[2][mat->col - 1] = vec->val[2];
	return mat;
}

//Returns pointer to 3*3 or 4*4 rotation matrix specified by rotation axis 'vec' and angle 'rad' (ccw +ve radians).
//Can result in gimbal lock.
//Dynamically allocates. User must free with rmlFreeMatrix.
rmlMatrix* rmlRotationMatrix(rmlVector* vec, double rad) {
	if (!(vec->size == 2 || vec->size == 3)) {
		printf("RML_ERROR_ROTATIONMATRIX: !(vec->size == 2 || vec->size == 3).\n");
		return NULL;
	}
	rmlMatrix* mat = rmlAllocateMatrix(vec->size + 1, vec->size + 1);
	mat->val[mat->row - 1][mat->col - 1] = 1.0;
	if (vec->size == 2) {
		mat->val[0][0] = cos(rad);	
		mat->val[0][1] = sin(rad);	
		mat->val[1][0] = -sin(rad);	
		mat->val[1][1] = cos(rad);	
	} else {
		mat->val[0][0] = cos(rad) + vec->val[0] * vec->val[0] * (1 - cos(rad));
		mat->val[0][1] = vec->val[0] * vec->val[1] * (1 - cos(rad)) - vec->val[2] * sin(rad);	
		mat->val[0][2] = vec->val[0] * vec->val[2] * (1 - cos(rad)) + vec->val[1] * sin(rad);
		mat->val[1][0] = vec->val[1] * vec->val[0] * (1 - cos(rad)) + vec->val[2] * sin(rad);
		mat->val[1][1] = cos(rad) + vec->val[1] * vec->val[1] * (1 - cos(rad));
		mat->val[1][2] = vec->val[1] * vec->val[2] * (1 - cos(rad)) - vec->val[0] * sin(rad);
		mat->val[2][0] = vec->val[2] * vec->val[0] * (1 - cos(rad)) - vec->val[1] * sin(rad);
		mat->val[2][1] = vec->val[2] * vec->val[1] * (1 - cos(rad)) + vec->val[0] * sin(rad);
		mat->val[2][2] = cos(rad) + vec->val[2] * vec->val[2] * (1 - cos(rad));
	}	
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

//TODO: Fix!!
//Modifies matrix pointed to by 'mat' from row-major to column-major or vice versa.
void rmlSwitchMajor(rmlMatrix* mat) {
}
