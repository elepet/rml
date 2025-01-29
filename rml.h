//Renderer Matrix Library

typedef struct {
	unsigned int row;
	unsigned int col;
	float** val; // Row-major matrix values.
} rmlMatrix;

typedef struct {
	unsigned int size;
	float* val;
} rmlVector;

rmlMatrix* rmlAllocateMatrix(unsigned int row, unsigned int col);
void rmlFreeMatrix(rmlMatrix* mat);
void rmlPrintMatrix(rmlMatrix* mat);
rmlVector* rmlAllocateVector(unsigned int size);
void rmlFreeVector(rmlVector* vec);
void rmlPrintVector(rmlVector* vec);
rmlMatrix* rmlDotMatrix(rmlMatrix* in1, rmlMatrix* in2);
float rmlDotVector(rmlVector* in1, rmlVector* in2);
rmlVector* rmlDotMatrixVector(rmlMatrix* mat, rmlVector* vec);
void rmlScalarByMatrix(float scalar, rmlMatrix* mat);
void rmlScalarByVector(float scalar, rmlVector* vec);
rmlVector* rmlCrossVector(rmlVector* in1, rmlVector* in2);
void rmlNormaliseVector(rmlVector* vec);
rmlMatrix* rmlIdentityMatrix(unsigned int size);
rmlMatrix* rmlScalingMatrix(rmlVector* vec);
rmlMatrix* rmlTranslationMatrix(rmlVector* vec);
rmlMatrix* rmlRotationMatrix(rmlVector* vec, double rad);
