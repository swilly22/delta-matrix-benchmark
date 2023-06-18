#include "rg_matrix.h"
#include <stdlib.h>

uint64_t iters = 5000; // number of operations to perform

// create a random matrix with nvals entries
RG_Matrix random_delta_matrix
(
	GrB_Index nrows,  // number of rows
	GrB_Index ncols,  // number of columns
	GrB_Index nvals   // number of entries
) {
	RG_Matrix A;
	RG_Matrix_new(&A, GrB_BOOL, nrows, ncols);

	// set nvals random entries
	for(GrB_Index i = 0; i < nvals; i++) {
		GrB_Index row = rand() % nrows;
		GrB_Index col = rand() % ncols;
		RG_Matrix_setElement_BOOL(A, row, col);
	}

	return A;
}

float bench_delta
(
	RG_Matrix A,
	GrB_Index *set_indices,
	GrB_Index *del_indices,
	uint64_t nset,
	uint64_t ndel
) {
	//--------------------------------------------------------------------------
	// perform random set/delete operations
	//--------------------------------------------------------------------------

	GrB_Index nrows;
	GrB_Index ncols;
	RG_Matrix_nrows(&nrows, A);
	RG_Matrix_ncols(&ncols, A);
	uint64_t set_idx = 0;
	uint64_t del_idx = 0;

	// start timer
	time_t t = clock();

	for(uint64_t i = 0; i < (nset + ndel) / 2; i++) {
		GrB_Index row;
		GrB_Index col;
		// set operation
		row = set_indices[set_idx++];
		col = set_indices[set_idx++];
		RG_Matrix_setElement_BOOL(A, row, col);
		RG_Matrix_wait(A, false);
		
		// delete operation
		row = del_indices[del_idx++];
		col = del_indices[del_idx++];
		RG_Matrix_removeElement_BOOL(A, row, col);
		RG_Matrix_wait(A, false);
		
		// read operation
		bool x;
		row = rand() % nrows;
		col = rand() % ncols;
		RG_Matrix_extractElement_BOOL(&x, A, row, col);
	}

	// return elapsed time in milliseconds
	float elapsed = ((float)(clock() - t))/CLOCKS_PER_SEC*1000;
	return elapsed;
}

float bench_sparse
(
	GrB_Matrix A,
	GrB_Index *set_indices,
	GrB_Index *del_indices,
	uint64_t nset,
	uint64_t ndel
) {
	//--------------------------------------------------------------------------
	// perform random set/delete operations
	//--------------------------------------------------------------------------

	GrB_Index nrows;
	GrB_Index ncols;
	GrB_Matrix_nrows(&nrows, A);
	GrB_Matrix_ncols(&ncols, A);
	uint64_t set_idx = 0;
	uint64_t del_idx = 0;

	// start timer
	time_t t = clock();

	for(uint64_t i = 0; i < iters; i++) {
		GrB_Index row;
		GrB_Index col;
		int op = i % 3;
		if(op == 0) {
			// set operation
			row = set_indices[set_idx++];
			col = set_indices[set_idx++];
			GrB_Matrix_setElement_BOOL(A, true, row, col);
			GrB_Matrix_wait(A, GrB_MATERIALIZE);
		} else if(op == 1) {
			// delete operation
			row = del_indices[del_idx++];
			col = del_indices[del_idx++];
			GrB_Matrix_removeElement(A, row, col);
			GrB_Matrix_wait(A, GrB_MATERIALIZE);
		} else {
			// read operation
			bool x;
			row = rand() % nrows;
			col = rand() % ncols;
			GrB_Matrix_extractElement_BOOL(&x, A, row, col);
		}
	}

	// return elapsed time in milliseconds
	float elapsed = ((float)(clock() - t))/CLOCKS_PER_SEC*1000;
	return elapsed;
}

// compute random sets for both set and del operations
void compute_change_set
(
	GrB_Matrix A,             // matrix compute change set for
	GrB_Index **set_indices,  // indices of set operations
	GrB_Index **del_indices,  // indices of delete operations
	uint64_t *nset,           // number of set operations
	uint64_t *ndel,           // number of delete operations
	uint64_t changes          // requested number of changes
) {
	// create a copy of the matrix
	GrB_Matrix _A;
	GrB_Matrix_dup(&_A, A);

	GrB_Index nrows;
	GrB_Index ncols;
	GrB_Matrix_nrows(&nrows, _A);
	GrB_Matrix_ncols(&ncols, _A);

	uint64_t _nset = changes / 3;  // 1/3 of changes are set operations
	uint64_t _ndel = changes / 3;  // 1/3 of changes are delete operations
	
	// compute set indices
	GrB_Index *_set_indices = malloc(sizeof(GrB_Index) * 2 * _nset);
	for(uint64_t i = 0; i < _nset; i++) {
		GrB_Index row = rand() % nrows;
		GrB_Index col = rand() % ncols;

		// save set indices
		_set_indices[i*2]   = row;
		_set_indices[i*2+1] = col;
		GrB_Matrix_setElement_BOOL(_A, true, row, col);
	}

	GrB_Index nvals;
	GrB_Matrix_nvals(&nvals, _A);

	// compute delete indices
	GrB_Index *_del_indices = malloc(sizeof(GrB_Index) * 2 * _ndel);

	GrB_Index *I = malloc(sizeof(GrB_Index) * nvals);
	GrB_Index *J = malloc(sizeof(GrB_Index) * nvals);
	GrB_Info info = GrB_Matrix_extractTuples_BOOL(I, J, NULL, &nvals, _A);

	for(uint64_t i = 0; i < _ndel; i++) {
		uint64_t idx = rand() % nvals;

		// save del indices
		_del_indices[i*2]   = I[idx];
		_del_indices[i*2+1] = J[idx];
	}

	free(I);
	free(J);

	*nset = _nset;
	*ndel = _ndel;
	*set_indices = _set_indices;
	*del_indices = _del_indices;
}

int main(int argv, char **argc) {
	// seed the random number generator
	srand(time(NULL));

	// initialize GraphBLAS
	GrB_init(GrB_NONBLOCKING);

	// matrix dimensions
	int n = 5;
	GrB_Index dimensions[15] = {
		10000,    10000,    500000,    // 10K  X 10K,  500K
		100000,   100000,   4000000,   // 100K X 100K, 4M
		1000000,  1000000,  30000000,  // 1M   X 1M,   30M
		10000000, 10000000, 50000000,  // 10M  X 10M,  50M
		80000000, 80000000, 300000000  // 80M  X 80M,  300M
	};

	// elapsed time for delta and sparse matrices
	float elapsed_delta[5]  = {0};
	float elapsed_sparse[5] = {0};

	RG_Matrix  delta;
	GrB_Matrix sparse;
	float      elapsed;

	for(int i = 0; i < n; i++) {
		GrB_Index nrows = dimensions[i * 3];
		GrB_Index ncols = dimensions[i * 3 + 1];
		GrB_Index nvals = dimensions[i * 3 + 2];

		printf("benchmarking %ld x %ld matrix\n", nrows, ncols);

		// create random sparse and delta matrix
		delta = random_delta_matrix(nrows, ncols, nvals);
		RG_Matrix_wait(delta, true);

		// create a sparse matrix from the delta matrix
		GrB_Matrix_dup(&sparse, RG_MATRIX_M(delta));

		RG_Matrix_nvals(&nvals, delta);
		printf("delta matrix nvals = %ld\n", nvals);

		GrB_Matrix_nvals(&nvals, sparse);
		printf("sparse matrix nvals = %ld\n", nvals);

		uint64_t nset;
		uint64_t ndel;
		GrB_Index *set_indices;
		GrB_Index *del_indices;
		uint64_t changes = iters;
		compute_change_set(sparse, &set_indices, &del_indices, &nset, &ndel,
				changes);

		elapsed = bench_delta(delta, set_indices, del_indices, nset, ndel);
		printf("delta matrix time = %f ms\n", elapsed);
		elapsed_delta[i] = elapsed;

		elapsed = bench_sparse(sparse, set_indices, del_indices, nset, ndel);
		printf("sparse matrix time = %f ms\n", elapsed);
		elapsed_sparse[i] = elapsed;

		// clean up
		free(set_indices);
		free(del_indices);
		GrB_free(&sparse);
		RG_Matrix_free(&delta);
	}

	// print results
	for(int i = 0; i < n; i++) {
		GrB_Index nrows = dimensions[i * 2];
		GrB_Index ncols = dimensions[i * 2 + 1];
		printf("delta: %ld x %ld, %f\n", nrows, ncols, elapsed_delta[i]);
		printf("sparse: %ld x %ld, %f\n", nrows, ncols, elapsed_sparse[i]);
		printf("\n");
	}

	GrB_finalize();
}

