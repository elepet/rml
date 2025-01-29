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
void rmlDotProduct(rmlMatrix* in1, rmlMatrix* in2, rmlMatrix* out);
rmlMatrix* rmlIdentityMatrix(unsigned int size);
rmlMatrix* rmlScalingMatrix(rmlVector* vec);
rmlMatrix* rmlTranslationMatrix(rmlVector* vec);
