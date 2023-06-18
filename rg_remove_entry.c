#include "rg_matrix.h"

GrB_Info RG_Matrix_removeEntry_UINT64
(
	RG_Matrix C,                    // matrix to remove entry from
	GrB_Index i,                    // row index
	GrB_Index j,                    // column index
	uint64_t  v,                    // value to remove
	bool     *entry_deleted         // is entry deleted
) {
	uint64_t    m_x;
	uint64_t    dp_x;
	GrB_Info    info;
	GrB_Type    type;
	bool        in_m        =  false;
	GrB_Matrix  m           =  RG_MATRIX_M(C);
	GrB_Matrix  dp          =  RG_MATRIX_DELTA_PLUS(C);
	GrB_Matrix  dm          =  RG_MATRIX_DELTA_MINUS(C);

	*entry_deleted = false;
	
	// entry should exists in either delta-plus or main
	// locate entry
	info = GrB_Matrix_extractElement(&m_x, m, i, j);
	in_m = (info == GrB_SUCCESS);

	//--------------------------------------------------------------------------
	// entry exists in 'M'
	//--------------------------------------------------------------------------

	if(in_m) {
		*entry_deleted = true;
		// mark deletion in delta minus
		info = GrB_Matrix_setElement(dm, true, i, j);
		info = RG_Matrix_removeElement_BOOL(C->transposed, j, i);
		RG_Matrix_setDirty(C);
		return info;
	}

	//--------------------------------------------------------------------------
	// entry exists in 'delta-plus'
	//--------------------------------------------------------------------------

	info = GrB_Matrix_extractElement(&dp_x, dp, i, j);
	if(info != GrB_SUCCESS) return info;

	*entry_deleted = true;
	info = GrB_Matrix_removeElement(dp, i, j);
	info = RG_Matrix_removeElement_BOOL(C->transposed, j, i);
	RG_Matrix_setDirty(C);
	return info;
}

