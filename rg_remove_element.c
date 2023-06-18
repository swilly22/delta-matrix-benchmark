#include "rg_matrix.h"

GrB_Info RG_Matrix_removeElement_BOOL
(
	RG_Matrix C,                    // matrix to remove entry from
	GrB_Index i,                    // row index
	GrB_Index j                     // column index
) {
	bool        m_x;
	bool        dm_x;
	bool        dp_x;
	GrB_Info    info;
	GrB_Type    type;
	bool        in_m  = false;
	bool        in_dp = false;
	bool        in_dm = false;
	GrB_Matrix  m     = RG_MATRIX_M(C);
	GrB_Matrix  dp    = RG_MATRIX_DELTA_PLUS(C);
	GrB_Matrix  dm    = RG_MATRIX_DELTA_MINUS(C);

	if(RG_MATRIX_MAINTAIN_TRANSPOSE(C)) {
		info = RG_Matrix_removeElement_BOOL(C->transposed, j, i);
		if(info != GrB_SUCCESS) {
			return info;
		}
	}

	//--------------------------------------------------------------------------
	// entry exists in 'M'
	//--------------------------------------------------------------------------

	info = GrB_Matrix_extractElement(&m_x, m, i, j);
	in_m = (info == GrB_SUCCESS);

	if(in_m) {
		// mark deletion in delta minus
		info = GrB_Matrix_setElement(dm, true, i, j);
		RG_Matrix_setDirty(C);
		return info;
	}

	//--------------------------------------------------------------------------
	// entry exists in 'delta-plus'
	//--------------------------------------------------------------------------

	// remove entry from 'dp'
	info = GrB_Matrix_removeElement(dp, i, j);
	RG_Matrix_setDirty(C);
	return info;
}

GrB_Info RG_Matrix_removeElement_UINT64
(
    RG_Matrix C,                    // matrix to remove entry from
    GrB_Index i,                    // row index
    GrB_Index j                     // column index
) {
	uint64_t    m_x;
	uint64_t    dm_x;
	uint64_t    dp_x;
	GrB_Info    info;
	GrB_Type    type;
	bool        in_m  = false;
	bool        in_dp = false;
	bool        in_dm = false;
	GrB_Matrix  m     = RG_MATRIX_M(C);
	GrB_Matrix  dp    = RG_MATRIX_DELTA_PLUS(C);
	GrB_Matrix  dm    = RG_MATRIX_DELTA_MINUS(C);

	if(RG_MATRIX_MAINTAIN_TRANSPOSE(C)) {
		info = RG_Matrix_removeElement_BOOL(C->transposed, j, i);
		if(info != GrB_SUCCESS) {
			return info;
		}
	}

	info = GrB_Matrix_extractElement(&m_x, m, i, j);
	in_m = (info == GrB_SUCCESS);

	info = GrB_Matrix_extractElement(&dp_x, dp, i, j);
	in_dp = (info == GrB_SUCCESS);

	info = GrB_Matrix_extractElement(&dm_x, dm, i, j);
	in_dm = (info == GrB_SUCCESS);

	// mask 'in_m' incase it is marked for deletion
	in_m = in_m && !(in_dm);

	// entry missing from both 'm' and 'dp'
	if(!(in_m || in_dp)) {
		return GrB_NO_VALUE;
	}

	// entry can't exists in both 'm' and 'dp'

	//--------------------------------------------------------------------------
	// entry exists in 'M'
	//--------------------------------------------------------------------------

	if(in_m) {
		// mark deletion in delta minus
		info = GrB_Matrix_setElement(dm, true, i, j);
	}

	//--------------------------------------------------------------------------
	// entry exists in 'delta-plus'
	//--------------------------------------------------------------------------

	if(in_dp) {
		// remove entry from 'dp'
		info = GrB_Matrix_removeElement(dp, i, j);
	}

	RG_Matrix_setDirty(C);
	return info;
}

