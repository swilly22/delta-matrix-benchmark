#include "rg_matrix.h"

// free multi-edge arrays GraphBLAS unary operation
static GrB_UnaryOp free_multi_edge_op = NULL;

// free RG_Matrix's internal matrices:
// M, delta-plus, delta-minus and transpose
void RG_Matrix_free
(
	RG_Matrix *C
) {
	RG_Matrix M = *C;

	GrB_Info info;

	if(RG_MATRIX_MAINTAIN_TRANSPOSE(M)) RG_Matrix_free(&M->transposed);

	GrB_Matrix m  = RG_MATRIX_M(M);
	GrB_Matrix dp = RG_MATRIX_DELTA_PLUS(M);

	info = GrB_Matrix_free(&M->matrix);
	info = GrB_Matrix_free(&M->delta_plus);
	info = GrB_Matrix_free(&M->delta_minus);

	pthread_mutex_destroy(&M->mutex);

	free(M);
	
	*C = NULL;
}

