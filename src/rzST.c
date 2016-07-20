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

void rzST(double *absz, double *M, double *RV, int *n, int *p, int *N, double *yVec, double *xiVec, double *psiVec, double *GVec, double *vVec){
 gsl_matrix *y = gsl_matrix_calloc(*n, *p);
 gsl_matrix *xi = gsl_matrix_calloc(*p, 1);
 gsl_matrix *psi = gsl_matrix_calloc(*p, 1);
 gsl_matrix *G = gsl_matrix_calloc(*p, *p);
 gsl_matrix *invG = gsl_matrix_alloc (*p, *p);
 gsl_permutation *perm = gsl_permutation_alloc(*p);
 gsl_matrix *tpsiinvG = gsl_matrix_calloc(1, *p);	// psi' * G^(-1)
 gsl_matrix *tpsiinvGpsi = gsl_matrix_calloc(1, 1);	// psi' * G^(-1) * psi
 gsl_matrix *ymxi = gsl_matrix_calloc(*p, 1);		// y[i,] - xi[k,]
 gsl_matrix *tpsiinvGymxi = gsl_matrix_calloc(1, 1);	// psi' * G^(-1) * (y - xi)

 double Vi;
 int s;
 size_t i, j, jj, k;
 // Define the data matrix y
 for(i = 0; i < *n; i++){
  for (j = 0; j < *p; j++){
   gsl_matrix_set(y, i, j, yVec[i * *p + j]);
  }
 }
 // Find M and RV
 for(k = 0; k < *N; k++){	// k = 0;
  for (j = 0; j < *p; j++){
   gsl_matrix_set(xi, j, 0, xiVec[k * *p + j]);
   gsl_matrix_set(psi, j, 0, psiVec[k * *p + j]);
   for (jj = 0; jj < *p; jj++){
    gsl_matrix_set(G, j, jj, GVec[k * *p * *p + jj * *p + j]);
   }
  }
//		double Gel = gsl_matrix_get(G, 0, 0);
//		*temp = Gel;
  gsl_linalg_LU_decomp (G, perm, &s);
  gsl_linalg_LU_invert(G, perm, invG);
/*  for (j = 0; j < *p; j++){
   for (jj = 0; jj < *p; jj++){
    double invGelement = gsl_matrix_get(invG, j, jj);
//    *temp = invGelement;
   }
  }*/
  gsl_matrix_set_zero(tpsiinvG);
  gsl_matrix_set_zero(tpsiinvGpsi);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, psi, invG, 1.0, tpsiinvG);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tpsiinvG, psi, 1.0, tpsiinvGpsi);
  double tpsiinvGpsiDouble = gsl_matrix_get(tpsiinvGpsi, 0, 0);
//  	double tpsiinvGDouble = gsl_matrix_get(tpsiinvG, 0, 0);
// 	*temp = tpsiinvGDouble;
  Vi = 1 / (1 + tpsiinvGpsiDouble);
  for(i = 0; i < *n; i++){
   RV[k * *n + i] = Vi;
   gsl_matrix_set_zero(tpsiinvGymxi);
   for (j = 0; j < *p; j++){
    gsl_matrix_set(ymxi, j, 0, gsl_matrix_get(y, i, j) - gsl_matrix_get(xi, j, 0));
   }
   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tpsiinvG, ymxi, 1.0, tpsiinvGymxi);
   double tpsiinvGymxiDouble = gsl_matrix_get(tpsiinvGymxi, 0, 0);
   M[k * *n + i] = Vi * sqrt(vVec[k * *n + i]) * tpsiinvGymxiDouble;
   leftTruncNorm(&(M[k * *n + i]), &(RV[k * *n + i]), &(absz[k * *n + i]));
  }
 }
}
