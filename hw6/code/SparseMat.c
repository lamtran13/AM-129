/* File: SparseMat.c
 * Author: Lam Tran
 * Purpose: Define a simplified CSR sparse matrix format
 */

#include "SparseMat.h"

/* Construct sparse matrix data structure in CSR format
 * with unknown non-zero layout 
 */
void CreateSparseMat(int m, int n, SparseMat *a_A)
{
  SparseMat A = malloc(sizeof(*A));
  /* Set overall dimensions */
  A->m = m;
  A->n = n;
  A->totnnz = 0;
  
  /* These are all arrays of length m, and can be allocated now */
  A->nnz = malloc(A->m*sizeof(*A->nnz));
  A->cols = malloc(A->m*sizeof(*A->cols));
  A->rowVals = malloc(A->m*sizeof(*A->rowVals));
  
  /* The length of these arrays won't be known until later, set NULL for now */
  A->layout = NULL;
  A->vals = NULL;

  /* Pass the created struct back out */
  *a_A = A;
}

/* Deallocate everything inside struct, then struct itself */
void DestroySparseMat(SparseMat *a_A)
{
  if(*a_A) {
    SparseMat A = *a_A;
    if(A->nnz) { free(A->nnz); }
    if(A->cols) { free(A->cols); }
    if(A->rowVals) { free(A->rowVals); }
    if(A->layout) { free(A->layout); }
    if(A->vals) { free(A->vals); }
    free(*a_A);
  }
}

/* Specify number of non-zeros per row, and set up internal links */
void SparseMatSetNNZ(SparseMat A, int *nnz)
{
  /* Copy in nnz, and track number of non-zeros in all rows */
  A->totnnz = 0;
  for(int i=0; i<A->m; i++) {
    A->nnz[i] = nnz[i];
    A->totnnz += nnz[i];
  }
  
  /* The length of these arrays is now known, allocate them now */
  A->layout = malloc(A->totnnz*sizeof(*A->layout));
  A->vals = calloc(A->totnnz, sizeof(*A->vals));

  /* Create aliases into the long arrays above */
  int ctr = 0;
  for(int i=0; i<A->m; i++) {
    A->cols[i] = &A->layout[ctr];
    A->rowVals[i] = &A->vals[ctr];
    ctr += A->nnz[i];
  }
}

/* Specify values and locations for a given row.
 * Must be compatible with corresponding entry in nnz given
 * during call to SparseMatSetNNZ */
void SparseMatSetRow(SparseMat A, int row, int *cols, double *vals)
{
  for(int j=0; j<A->nnz[row]; j++) {
    A->cols[row][j] = cols[j];
    A->rowVals[row][j] = vals[j];
  }
}

/* Perform sparse matrix vector product Ax and store into y.
 * Dimensions not checked, and need to be compatible
 */
void SparseMatMult(SparseMat A, double *x, double *y)
{
  for(int i=0; i<A->m; i++){
      y[i]=0;
    for(int k=0; k<A->nnz[i]-1; k++){
      y[i] += A->rowVals[i][k] * x[A->cols[i][k]];
  }
  }
}

/* Supporting functions */
void SparseMatPrintInfo(SparseMat A, const char *name)
{
  if(name) {
    printf("Matrix %s has:\n", name);
  } else {
    printf("Matrix has:\n");
  }
  printf("  %d rows\n", A->m);
  printf("  %d columns\n", A->n);
  printf("  %d total non-zeros\n", A->totnnz);
}

void SparseMatPrintAll(SparseMat A, const char *name)
{
  SparseMatPrintInfo(A, name);
  for(int i=0; i<A->m; i++) {
    if(A->nnz[i]) {
      printf(" Row %d has %d nz cols: (%d", i, A->nnz[i], A->cols[i][0]);
    } else {
      printf(" Row %d cols: (", i);
    }
    for(int c=1; c<A->nnz[i]; c++) {
      printf(",%d", A->cols[i][c]);
    }
    if(A->nnz[i]) {
      printf(") vals: (%e", A->rowVals[i][0]);
    } else {
      printf(") vals: (");
    }
    for(int c=1; c<A->nnz[i]; c++) {
      printf(",%e", A->rowVals[i][c]);
    }
    printf(")\n");
  }
}
