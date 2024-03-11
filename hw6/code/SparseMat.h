/* File: SparseMat.h
 * Author: Lam Tran
 * Purpose: Define a simplified CSR sparse matrix format
 */

#include <stdlib.h>
#include <stdio.h>

#ifndef SPARSEMAT_H
#define SPARSEMAT_H

typedef struct _p_SparseMat *SparseMat;
struct _p_SparseMat {
  int m, n, totnnz;
  int *nnz, *layout;
  int **cols;
  double *vals;
  double **rowVals;
};

/* Construction - destruction */
void CreateSparseMat(int,int,SparseMat*);
void DestroySparseMat(SparseMat*);

/* Filling the matrix */
void SparseMatSetNNZ(SparseMat,int*);
void SparseMatSetRow(SparseMat,int,int*,double*);

/* Matrix vector products */
void SparseMatMult(SparseMat,double*,double*);

/* Supporting functions */
void SparseMatPrintInfo(SparseMat,const char*);
void SparseMatPrintAll(SparseMat,const char*);

#endif
