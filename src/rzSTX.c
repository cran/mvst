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

void rzSTX(double *absz, double *M, double *RV, int *n, int *p, int *k, int *N, double *yVec, double *XVec, double *BVec, double *psiVec, double *GVec, double *vVec){
 gsl_matrix *y = gsl_matrix_calloc(*n, *p);
 gsl_matrix *psi = gsl_matrix_calloc(*p, 1);
 gsl_matrix *G = gsl_matrix_calloc(*p, *p);
 gsl_matrix *invG = gsl_matrix_alloc (*p, *p);
 gsl_permutation *perm = gsl_permutation_alloc(*p);
 gsl_matrix *tpsiinvG = gsl_matrix_calloc(1, *p);	// psi' * G^(-1)
 gsl_matrix *tpsiinvGpsi = gsl_matrix_calloc(1, 1);	// psi' * G^(-1) * psi
 gsl_matrix *ymBX = gsl_matrix_calloc(*p, 1);		// y - BX
 gsl_matrix *tpsiinvGymBX = gsl_matrix_calloc(1, 1);	// psi' * G^(-1) * (y - BX)

 double Vi[1];
 double tpsiinvGpsiDouble[1];
 double tpsiinvGymBXDouble[1];
 int s;
 size_t i, j, l, jj, iN;
 double BXj;
 // Define the data matrix y
 for(i = 0; i < *n; i++){
  for(j = 0; j < *p; j++){
   gsl_matrix_set(y, i, j, yVec[i * *p + j]);
  }
 }
 // Find M and RV
 for(iN = 0; iN < *N; iN++){	// iN = 0;
  for (j = 0; j < *p; j++){
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
   gsl_matrix_set_zero(tpsiinvGymBX);
   for (j = 0; j < *p; j++){
    BXj = 0;
    for (l = 0; l < *k; l++){
     BXj += BVec[iN * (*k * *p) + j * *k + l] * XVec[i * *k + l];
    }
    gsl_matrix_set(ymBX, j, 0, gsl_matrix_get(y, i, j) - BXj);
   }
   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tpsiinvG, ymBX, 1.0, tpsiinvGymBX);
   tpsiinvGymBXDouble[0] = gsl_matrix_get(tpsiinvGymBX, 0, 0);
   M[iN * *n + i] = Vi[0] * sqrt(vVec[iN * *n + i]) * tpsiinvGymBXDouble[0];
   leftTruncNorm(&(M[iN * *n + i]), &(RV[iN * *n + i]), &(absz[iN * *n + i]));
  }
 }
 gsl_matrix_free(y);
 gsl_matrix_free(psi);
 gsl_matrix_free(G);
 gsl_matrix_free(invG);
 gsl_permutation_free(perm);
 gsl_matrix_free(tpsiinvG);
 gsl_matrix_free(tpsiinvGpsi);
 gsl_matrix_free(ymBX);
 gsl_matrix_free(tpsiinvGymBX);
}

