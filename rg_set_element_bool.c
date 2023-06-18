#include "rg_matrix.h"

GrB_Info RG_Matrix_setElement_BOOL      // C (i,j) = x
(
    RG_Matrix C,                        // matrix to modify
    GrB_Index i,                        // row index
    GrB_Index j                         // column index
) {
	bool v;
	GrB_Info info;

	GrB_Matrix m  = RG_MATRIX_M(C);
	GrB_Matrix dp = RG_MATRIX_DELTA_PLUS(C);
	GrB_Matrix dm = RG_MATRIX_DELTA_MINUS(C);

	bool already_allocated   = false;  // M[i,j] exists
	bool marked_for_deletion = false;  // dm[i,j] exists

	if(RG_MATRIX_MAINTAIN_TRANSPOSE(C)) {
		info = RG_Matrix_setElement_BOOL(C->transposed, j, i);
	}

	info = GrB_Matrix_extractElement(&v, dm, i, j);
	marked_for_deletion = (info == GrB_SUCCESS);

	if(marked_for_deletion) {
		// unset delta-minus m already assign to true
		info = GrB_Matrix_removeElement(dm, i, j);
	} else {
		info = GrB_Matrix_extractElement(&v, m, i, j);
		already_allocated = (info == GrB_SUCCESS);

		if(!already_allocated) {
			// update entry to dp[i, j]
			info = GrB_Matrix_setElement_BOOL(dp, true, i, j);
		}
	}

	RG_Matrix_setDirty(C);

	return info;
}

