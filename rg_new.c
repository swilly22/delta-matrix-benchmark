#include "rg_matrix.h"

static GrB_Info _RG_Matrix_init
(
	RG_Matrix A,
	GrB_Type type,
	GrB_Index nrows,
	GrB_Index ncols
) {
	GrB_Info info;
	A->dirty = false;

	//--------------------------------------------------------------------------
	// create m, delta-plus and delta-minus
	//--------------------------------------------------------------------------

	//--------------------------------------------------------------------------
	// m, can be either hypersparse or sparse
	//--------------------------------------------------------------------------
	info = GrB_Matrix_new(&A->matrix, type, nrows, ncols);
	info = GxB_set(A->matrix, GxB_SPARSITY_CONTROL, GxB_SPARSE | GxB_HYPERSPARSE);

	//--------------------------------------------------------------------------
	// delta-plus, always hypersparse
	//--------------------------------------------------------------------------
	info = GrB_Matrix_new(&A->delta_plus, type, nrows, ncols);
	info = GxB_set(A->delta_plus, GxB_SPARSITY_CONTROL, GxB_HYPERSPARSE);
	info = GxB_set(A->delta_plus, GxB_HYPER_SWITCH, GxB_ALWAYS_HYPER);

	//--------------------------------------------------------------------------
	// delta-minus, always hypersparse
	//--------------------------------------------------------------------------
	info = GrB_Matrix_new(&A->delta_minus, GrB_BOOL, nrows, ncols);
	info = GxB_set(A->delta_minus, GxB_SPARSITY_CONTROL, GxB_HYPERSPARSE);
	info = GxB_set(A->delta_minus, GxB_HYPER_SWITCH, GxB_ALWAYS_HYPER);

	return info;
}

// creates a new matrix
GrB_Info RG_Matrix_new
(
	RG_Matrix *A,
	GrB_Type type,
	GrB_Index nrows,
	GrB_Index ncols
) {
	GrB_Info info;
	RG_Matrix matrix = calloc(1, sizeof(_RG_Matrix));

	//--------------------------------------------------------------------------
	// input validations
	//--------------------------------------------------------------------------

	// supported types: boolean and uint64
	info = _RG_Matrix_init(matrix, type, nrows, ncols);

	//--------------------------------------------------------------------------
	// create transpose matrix if required
	//--------------------------------------------------------------------------

	if(type == GrB_UINT64) {
		matrix->transposed = calloc(1, sizeof(_RG_Matrix));
		info = _RG_Matrix_init(matrix->transposed, GrB_BOOL, ncols, nrows);
	}

	int mutex_res = pthread_mutex_init(&matrix->mutex, NULL);

	*A = matrix;
	return info;
}

