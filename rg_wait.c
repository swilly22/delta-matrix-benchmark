#include "rg_matrix.h"

static inline void _SetUndirty
(
	RG_Matrix C
) {
	C->dirty = false;

	if(RG_MATRIX_MAINTAIN_TRANSPOSE(C)) {
		C->transposed->dirty = false;
	}
}

static GrB_Info RG_Matrix_sync
(
	RG_Matrix C
) {
	GrB_Matrix      m    = RG_MATRIX_M(C);
	GrB_Matrix      dp   = RG_MATRIX_DELTA_PLUS(C);
	GrB_Matrix      dm   = RG_MATRIX_DELTA_MINUS(C);
	GrB_Descriptor  desc = GrB_NULL;
	GrB_Matrix      mask = GrB_NULL;

	GrB_Info info;
	GrB_Index dp_nvals;
	GrB_Index dm_nvals;

	//--------------------------------------------------------------------------
	// determin change set
	//--------------------------------------------------------------------------

	GrB_Matrix_nvals(&dp_nvals, dp);
	GrB_Matrix_nvals(&dm_nvals, dm);

	bool additions = dp_nvals > 0;
	bool deletions = dm_nvals > 0;

	//--------------------------------------------------------------------------
	// perform deletions
	//--------------------------------------------------------------------------

	if(deletions) {
		info = GrB_transpose(m, dm, GrB_NULL, m, GrB_DESC_RSCT0);

		// clear delta minus
		info = GrB_Matrix_clear(dm);
	}

	//--------------------------------------------------------------------------
	// perform additions
	//--------------------------------------------------------------------------

	if(additions) {
		GrB_Type t;
		GrB_Semiring s;
		info = GxB_Matrix_type(&t, m);

		s = (t == GrB_BOOL) ? GxB_ANY_PAIR_BOOL : GxB_ANY_PAIR_UINT64;
		info = GrB_Matrix_eWiseAdd_Semiring(m, NULL, NULL, s, m, dp, NULL);

		// clear delta plus
		info = GrB_Matrix_clear(dp);
	}

	//--------------------------------------------------------------------------
	// validate that both delta-plus and delta-minus are cleared
	//--------------------------------------------------------------------------

	GrB_Index nvals;
	GrB_Matrix_nvals(&nvals, dp);
	GrB_Matrix_nvals(&nvals, dm);

	// wait on all 3 matrices
	info = GrB_wait(m,  GrB_MATERIALIZE);
	info = GrB_wait(dm, GrB_MATERIALIZE);
	info = GrB_wait(dp, GrB_MATERIALIZE);

	return info;
}

GrB_Info RG_Matrix_wait
(
	RG_Matrix A,
	bool force_sync
) {
	if(RG_MATRIX_MAINTAIN_TRANSPOSE(A)) {
		RG_Matrix_wait(A->transposed, force_sync);
	}
	
	GrB_Info   info        = GrB_SUCCESS;
	GrB_Matrix m           = RG_MATRIX_M(A);
	GrB_Matrix delta_plus  = RG_MATRIX_DELTA_PLUS(A);
	GrB_Matrix delta_minus = RG_MATRIX_DELTA_MINUS(A);

	// check if merge is required
	GrB_Index delta_plus_nvals;
	GrB_Index delta_minus_nvals;
	GrB_Matrix_nvals(&delta_plus_nvals, delta_plus);
	GrB_Matrix_nvals(&delta_minus_nvals, delta_minus);

	uint64_t delta_max_pending_changes = 10000;

	if(force_sync ||
	   delta_plus_nvals + delta_minus_nvals >= delta_max_pending_changes) {
		info = RG_Matrix_sync(A);
	} else {
		// wait on 'm', in most cases 'm' won't contain any pending work
		// but it might need to build its internal hyper-hash
		info = GrB_wait(m, GrB_MATERIALIZE);
		info = GrB_wait(delta_plus, GrB_MATERIALIZE);
		info = GrB_wait(delta_minus, GrB_MATERIALIZE);
	}

	_SetUndirty(A);

	return info;
}

