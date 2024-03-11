/* File: Driver.c
 * Author: Lam Tran
 * Purpose: This program computes the stationary distribution for a particle
 *          on a one-dimensional lattice. The probability to stay stationary
 *          depends on position, and their is a small bias in which in the
 *          probability of hopping left vs. right.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "SparseMat.h"

#define ITER_MAX 100000
#define THRESH 1.e-6
#define BIAS 1.e-3

/* Probability to stay on a lattice point */
double Stay(double x);

/* Probability to hop left */
double HopL(double x);

/* Probability to hop right */
double HopR(double x);

/* Normalize a vector by it's l1 norm */
void normalize_l1(int m, double *x);

/* Normalize a vector by it's l2 norm */
void normalize_l2(int m, double *x);

/* Compute the l1 norm of the difference of 2 vectors */
double normDiff(int m, double *u, double *v);

int main(int argc, char **argv)
{
  /* Set number of grid points */
  int m = argc>1 ? atoi(argv[1]) : 50;
  
  /* Create domain, (0,1) exclusive of boundaries */
  double dx = 1./(m+1);
  double x[m];
  for(int i=0; i<m; i++) {
    x[i] = (i+1)*dx;
  }

  /* Create vectors, intialize u to all ones */
  double u[m];
  double v[m];
  for(int i=0; i<m; i++) {
    u[i] = 1.;
  }
  normalize_l2(m,u);

  /* Create empty sparse matrix */
  SparseMat mat = NULL;
  CreateSparseMat(m,m,&mat);

  /* Set number of non-zeros */
  /* ------------ Add nnz array here! --------- */
  int nnz[m];
  nnz[0] = 2;
  nnz[m-1] = 2;
  for(int i=1; i<m-1; i++) {
    nnz[i] = 3;
  }
  SparseMatSetNNZ(mat,nnz);

  /* Fill first row */
  int colsFR[] = {0,1};
  double valsFR[] = {1.-HopR(x[0]),HopL(x[1])};
  SparseMatSetRow(mat,0,colsFR,valsFR);

  /* Fill in last two rows */
  int colsLR[] = {m-2,m-1};
  double valsLR[] = {HopR(x[m-2]),1.-HopL(x[m-1])};
  SparseMatSetRow(mat,m-1,colsLR,valsLR);
  
  /* Fill interior rows */

  for(int i=1; i<m-1; i++) {
    int cols[] = {i-1, i, i+1};
    double vals[] = {HopR(x[i-1]),Stay(x[i]), HopL(x[i+1])};
    SparseMatSetRow(mat,i,cols,vals);
   }
  /* Nothing below here needs to be changed */
  
  /* Run power iteration */
  for(int nI=0; nI<ITER_MAX; nI++) {
    if(nI%2==0) {
      SparseMatMult(mat,u,v);
      normalize_l2(m,v);
    } else {
      SparseMatMult(mat,v,u);
      normalize_l2(m,u);
    }
    double diff = normDiff(m,u,v);
    if(diff<THRESH) {
      printf("Converged in %d iterations\n",nI+1);
      break;
    } else if(nI==ITER_MAX-1) {
      printf("Failed to converge, diff = %e\n",diff);
    }
  }
  /* Normalize instead by l1 to make this a PMF */
  normalize_l1(m,v);

  /* Save grid and vector to file */
  FILE *fp = fopen("dist.dat","w");
  for(int i=0; i<m; i++) {
    fprintf(fp, "%e %e\n",x[i],v[i]);
  }
  fclose(fp);

  /* Clean up */
  DestroySparseMat(&mat);
  
  return 0;
}

/* Probability to stay on a lattice point */
double Stay(double x)
{
  return 0.75*exp(-fabs(x-0.25));
}

/* Probability to hop left */
double HopL(double x)
{
  return (0.5 - BIAS)*(1.-Stay(x));
}

/* Probability to hop right */
double HopR(double x)
{
  return (0.5 + BIAS)*(1.-Stay(x));
}

/* Normalize a vector by it's l2 norm */
void normalize_l2(int m, double *x)
{
  double norm = 0;
  for(int i=0; i<m; i++) {
    norm += x[i]*x[i];
  }
  norm = sqrt(norm);
  for(int i=0; i<m; i++) {
    x[i] /= norm;
  }
}

/* Normalize a vector by it's l1 norm */
void normalize_l1(int m, double *x)
{
  double norm = 0;
  for(int i=0; i<m; i++) {
    norm += fabs(x[i]);
  }
  for(int i=0; i<m; i++) {
    x[i] /= norm;
  }
}

/* Compute the l2 norm of the difference of 2 vectors */
double normDiff(int m, double *u, double *v) 
{
  double norm = 0;
  for(int i=0; i<m; i++) {
    norm += pow(u[i]-v[i],2);
  }
  return sqrt(norm);
}
