#include "rg_matrix.h"

void RG_Matrix_setDirty
(
	RG_Matrix C
) {
	C->dirty = true;
	if(RG_MATRIX_MAINTAIN_TRANSPOSE(C)) C->transposed->dirty = true;
}

RG_Matrix RG_Matrix_getTranspose
(
	const RG_Matrix C
) {
	return C->transposed;
}

bool RG_Matrix_isDirty
(
	const RG_Matrix C
) {

	if(C->dirty) {
		return true;
	}

	bool pending_M;
	bool pending_DP;
	bool pending_DM;

	GxB_Matrix_Pending(RG_MATRIX_M(C), &pending_M);
	GxB_Matrix_Pending(RG_MATRIX_DELTA_PLUS(C), &pending_DP);
	GxB_Matrix_Pending(RG_MATRIX_DELTA_MINUS(C), &pending_DM);

	return (pending_M | pending_DM | pending_DP);
}

// checks if C is fully synced
// a synced delta matrix does not contains any entries in
// either its delta-plus and delta-minus internal matrices
bool RG_Matrix_Synced
(
	const RG_Matrix C  // matrix to inquery
) {

	// quick indication, if the matrix is marked as dirty that means
	// entires exists in either DP or DM
	if(C->dirty) {
		return false;
	}

	GrB_Index dp_nvals;
	GrB_Index dm_nvals;
	GrB_Matrix_nvals(&dp_nvals, RG_MATRIX_DELTA_PLUS(C));
	GrB_Matrix_nvals(&dm_nvals, RG_MATRIX_DELTA_MINUS(C));

	return ((dp_nvals + dm_nvals) == 0);
}

// locks the matrix
void RG_Matrix_Lock
(
	RG_Matrix C
) {
	pthread_mutex_lock(&C->mutex);
}

// unlocks the matrix
void RG_Matrix_Unlock
(
	RG_Matrix C
) {
	pthread_mutex_unlock(&C->mutex);
}

GrB_Info RG_Matrix_nrows
(
	GrB_Index *nrows,
	const RG_Matrix C
) {
	GrB_Matrix m = RG_MATRIX_M(C);
	return GrB_Matrix_nrows(nrows, m);
}

GrB_Info RG_Matrix_ncols
(
	GrB_Index *ncols,
	const RG_Matrix C
) {
	GrB_Matrix m = RG_MATRIX_M(C);
	return GrB_Matrix_ncols(ncols, m);
}

GrB_Info RG_Matrix_nvals    // get the number of entries in a matrix
(
    GrB_Index *nvals,       // matrix has nvals entries
    const RG_Matrix A       // matrix to query
) {
	GrB_Matrix  m;
	GrB_Matrix  dp;
	GrB_Matrix  dm;
	GrB_Info    info;

	GrB_Index m_nvals  = 0;
	GrB_Index dp_nvals = 0;
	GrB_Index dm_nvals = 0;

	// nvals = nvals(M) + nvals(DP) - nvals(DM)

	m  = RG_MATRIX_M(A);
	dp = RG_MATRIX_DELTA_PLUS(A);
	dm = RG_MATRIX_DELTA_MINUS(A);

	info = GrB_Matrix_nvals(&m_nvals, m);
	info = GrB_Matrix_nvals(&dp_nvals, dp);
	info = GrB_Matrix_nvals(&dm_nvals, dm);

	*nvals = m_nvals + dp_nvals - dm_nvals;
	return info;
}

GrB_Info RG_Matrix_clear
(
    RG_Matrix A
) {
	GrB_Matrix  m           = RG_MATRIX_M(A);
	GrB_Info    info        = GrB_SUCCESS;
	GrB_Matrix  delta_plus  = RG_MATRIX_DELTA_PLUS(A);
	GrB_Matrix  delta_minus = RG_MATRIX_DELTA_MINUS(A);

	info = GrB_Matrix_clear(m);
	info = GrB_Matrix_clear(m);
	info = GrB_Matrix_clear(m);

	A->dirty = false;
	if(RG_MATRIX_MAINTAIN_TRANSPOSE(A)) A->transposed->dirty = false;

	return info;
}

GrB_Info RG_Matrix_type
(
	GrB_Type *type,
	RG_Matrix A
) {
	GrB_Matrix M = RG_MATRIX_M(A);
	GrB_Info info = GxB_Matrix_type(type, M);
	return info;
}

