#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include "Utils.h"

void rzSN(double *absz, double *M, double *RV, int *n, int *p, int *N, double *yVec, double *xiVec, double *psiVec, double *GVec){
 gsl_matrix *y = gsl_matrix_calloc(*n, *p);
 gsl_matrix *xi = gsl_matrix_calloc(*p, 1);
 gsl_matrix *psi = gsl_matrix_calloc(*p, 1);
 gsl_matrix *G = gsl_matrix_calloc(*p, *p);
 gsl_matrix *invG = gsl_matrix_alloc (*p, *p);
 gsl_permutation *perm = gsl_permutation_alloc(*p);
 gsl_matrix *tpsiinvG = gsl_matrix_calloc(1, *p);	// psi' * G^(-1)
 gsl_matrix *tpsiinvGpsi = gsl_matrix_calloc(1, 1);	// psi' * G^(-1) * psi
 gsl_matrix *ymxi = gsl_matrix_calloc(*p, 1);		// y[i,] - xi[iN,]
 gsl_matrix *tpsiinvGymxi = gsl_matrix_calloc(1, 1);	// psi' * G^(-1) * (y - xi)

 double Vi[1];
 double tpsiinvGpsiDouble[1];
 double tpsiinvGymxiDouble[1];
 int s;
 size_t i, j, jj, iN;
 // Define the data matrix y
 for(i = 0; i < *n; i++){
  for (j = 0; j < *p; j++){
   gsl_matrix_set(y, i, j, yVec[i * *p + j]);
  }
 }
 // Find M and RV
 for(iN = 0; iN < *N; iN++){	// iN = 0;
  for (j = 0; j < *p; j++){
   gsl_matrix_set(xi, j, 0, xiVec[iN * *p + j]);
   gsl_matrix_set(psi, j, 0, psiVec[iN * *p + j]);
   for (jj = 0; jj < *p; jj++){
    gsl_matrix_set(G, j, jj, GVec[iN * *p * *p + jj * *p + j]);
   }
  }
  gsl_linalg_LU_decomp (G, perm, &s);
  gsl_linalg_LU_invert(G, perm, invG);
  gsl_matrix_set_zero(tpsiinvG);
  gsl_matrix_set_zero(tpsiinvGpsi);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, psi, invG, 1.0, tpsiinvG);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tpsiinvG, psi, 1.0, tpsiinvGpsi);
  tpsiinvGpsiDouble[0] = gsl_matrix_get(tpsiinvGpsi, 0, 0);
  Vi[0] = 1 / (1 + tpsiinvGpsiDouble[0]);
  for(i = 0; i < *n; i++){
   RV[iN * *n + i] = Vi[0];
   gsl_matrix_set_zero(tpsiinvGymxi);
   for (j = 0; j < *p; j++){
    gsl_matrix_set(ymxi, j, 0, gsl_matrix_get(y, i, j) - gsl_matrix_get(xi, j, 0));
   }
   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tpsiinvG, ymxi, 1.0, tpsiinvGymxi);
   tpsiinvGymxiDouble[0] = gsl_matrix_get(tpsiinvGymxi, 0, 0);
   M[iN * *n + i] = Vi[0] * tpsiinvGymxiDouble[0];
   leftTruncNorm(&(M[iN * *n + i]), &(RV[iN * *n + i]), &(absz[iN * *n + i]));
  }
 }
 gsl_matrix_free(y);
 gsl_matrix_free(xi);
 gsl_matrix_free(psi);
 gsl_matrix_free(G);
 gsl_matrix_free(invG);
 gsl_matrix_free(ymxi);
 gsl_permutation_free(perm);
 gsl_matrix_free(tpsiinvG);
 gsl_matrix_free(tpsiinvGpsi);
 gsl_matrix_free(tpsiinvGymxi);
}

