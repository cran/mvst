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

void rvT(double *values, double *propdens, int *n, int *p, int *N, double *yVec, double *nu, double *xiVec, double *GVec){
 gsl_matrix *y = gsl_matrix_calloc(*n, *p);
 gsl_matrix *xi = gsl_matrix_calloc(*p, 1);
 gsl_matrix *G = gsl_matrix_calloc(*p, *p);
 gsl_matrix *invG = gsl_matrix_calloc(*p, *p);
 gsl_permutation *perm = gsl_permutation_alloc(*p);
 gsl_matrix *ymxi = gsl_matrix_calloc(*p, 1);		// y[i,] - xi[,j]
 gsl_matrix *Product1 = gsl_matrix_calloc(1, *p);
 gsl_matrix *Product2 = gsl_matrix_calloc(1, 1);
 double Ai[1];
 double Alphav[1];

 size_t i, j;
 for (j = 0; j < *p; j++){
  for(i = 0; i < *n; i++){
   gsl_matrix_set(y, i, j, yVec[i * *p + j]);
  }
 }

 size_t jj, iN;
 int s;
 int counter;
 counter = 0;
 for(iN = 0; iN < *N; iN++){	// iN = 0;
  Alphav[0] = 0.5 * (nu[iN] + *p);
  for (j = 0; j < *p; j++){
   gsl_matrix_set(xi, j, 0, xiVec[iN * *p + j]);
//   gsl_matrix_set(psi, j, 0, psiVec[iN * *p + j]);
   for (jj = 0; jj < *p; jj++){
//  gsl_matrix_set(G, j, jj, GVec[jj * *p + j]);
    gsl_matrix_set(G, j, jj, GVec[iN * *p * *p + jj * *p + j]);
   }
  }
  gsl_linalg_LU_decomp (G, perm, &s);
  gsl_linalg_LU_invert(G, perm, invG);
  for(i = 0; i < *n; i++){
//   zi = zVec[iN * *n + i];
   for (j = 0; j < *p; j++){
    gsl_matrix_set(ymxi, j, 0, gsl_matrix_get(y, i, j) - gsl_matrix_get(xi, j, 0));
   }
   gsl_matrix_set_zero(Product1);
   gsl_matrix_set_zero(Product2);
//   gsl_matrix_set_zero(Product3);
   gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, ymxi, invG, 1.0, Product1);
   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Product1, ymxi, 1.0, Product2);
//   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Product1, psi, 1.0, Product3);
   double productA = gsl_matrix_get(Product2, 0, 0);
//   double productB = gsl_matrix_get(Product3, 0, 0);
   Ai[0] = 0.5 * (nu[iN] + productA);
/*   Bi[0] = - productB * fabsl(zi);
   g[0] = - 0.5 * Bi[0] * expl(gsl_sf_lngamma(Alphav[0] + 0.5) - gsl_sf_lngamma(Alphav[0]));
   sqrtBetav[0] = (- g[0] + sqrt(g[0] * g[0] + 4 * Ai[0] * Alphav[0] * Alphav[0])) / (2 * Alphav[0]);
   Betav[0] = sqrtBetav[0] * sqrtBetav[0]; */
   GetRNGstate();
   values[counter] = rgamma(Alphav[0], 1/Ai[0]);
   PutRNGstate();
   propdens[iN] = propdens[iN] + dgamma(values[counter], Alphav[0], 1/Ai[0], 1);
   counter++;
  }
 }

 gsl_matrix_free(y);
 gsl_matrix_free(xi);
 gsl_matrix_free(G);
 gsl_matrix_free(invG);
 gsl_matrix_free(ymxi);
 gsl_permutation_free(perm);
 gsl_matrix_free(Product1);
 gsl_matrix_free(Product2);
}


void rvTX(double *values, double *propdens, int *n, int *p, int *k, int *N, double *yVec, double *XVec, double *nu, double *BVec, double *GVec){
 gsl_matrix *y = gsl_matrix_calloc(*n, *p);
 gsl_matrix *G = gsl_matrix_calloc(*p, *p);
 gsl_matrix *invG = gsl_matrix_calloc(*p, *p);
 gsl_permutation *perm = gsl_permutation_alloc(*p);
 gsl_matrix *ymBX = gsl_matrix_calloc(*p, 1);		// y - BX
 gsl_matrix *Product1 = gsl_matrix_calloc(1, *p);
 gsl_matrix *Product2 = gsl_matrix_calloc(1, 1);
 double Ai[1];
 double Alphav[1];

 size_t i, j, l;
 for (j = 0; j < *p; j++){
  for(i = 0; i < *n; i++){
   gsl_matrix_set(y, i, j, yVec[i * *p + j]);
  }
 }

 size_t jj, iN;
 int s;
// double zi;
 double BXj;
 int counter;
 counter = 0;
 for(iN = 0; iN < *N; iN++){	// iN = 0;
  Alphav[0] = 0.5 * (nu[iN] + *p);
  for (j = 0; j < *p; j++){
//   gsl_matrix_set(xi, j, 0, xiVec[iN * *p + j]);
//   gsl_matrix_set(psi, j, 0, psiVec[iN * *p + j]);
   for (jj = 0; jj < *p; jj++){
//  gsl_matrix_set(G, j, jj, GVec[jj * *p + j]);
    gsl_matrix_set(G, j, jj, GVec[iN * *p * *p + jj * *p + j]);
   }
  }
  gsl_linalg_LU_decomp(G, perm, &s);
  gsl_linalg_LU_invert(G, perm, invG);
  for(i = 0; i < *n; i++){
//   zi = zVec[iN * *n + i];
   for(j = 0; j < *p; j++){
    BXj = 0;
    for (l = 0; l < *k; l++){
     BXj += BVec[iN * (*k * *p) + j * *k + l] * XVec[i * *k + l];
    }
    gsl_matrix_set(ymBX, j, 0, gsl_matrix_get(y, i, j) - BXj);
//    gsl_matrix_set(ymxi, j, 0, gsl_matrix_get(y, i, j) - gsl_matrix_get(xi, j, 0));
   }
   gsl_matrix_set_zero(Product1);
   gsl_matrix_set_zero(Product2);
//   gsl_matrix_set_zero(Product3);
   gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, ymBX, invG, 1.0, Product1);
   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Product1, ymBX, 1.0, Product2);
//   gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, ymxi, invG, 1.0, Product1);
//   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Product1, ymxi, 1.0, Product2);
//   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Product1, psi, 1.0, Product3);
   double productA = gsl_matrix_get(Product2, 0, 0);
//   double productB = gsl_matrix_get(Product3, 0, 0);
   Ai[0] = 0.5 * (nu[iN] + productA);
/*   Bi[0] = - productB * fabsl(zi);
   g[0] = - 0.5 * Bi[0] * expl(gsl_sf_lngamma(Alphav[0] + 0.5) - gsl_sf_lngamma(Alphav[0]));
   sqrtBetav[0] = (- g[0] + sqrt(g[0] * g[0] + 4 * Ai[0] * Alphav[0] * Alphav[0])) / (2 * Alphav[0]);
   Betav[0] = sqrtBetav[0] * sqrtBetav[0]; */
   GetRNGstate();
   values[counter] = rgamma(Alphav[0], 1/Ai[0]);
   PutRNGstate();
   propdens[iN] = propdens[iN] + dgamma(values[counter], Alphav[0], 1/Ai[0], 1);
   counter++;
  }
 }

 gsl_matrix_free(y);
 gsl_matrix_free(G);
 gsl_matrix_free(invG);
 gsl_matrix_free(ymBX);
 gsl_permutation_free(perm);
 gsl_matrix_free(Product1);
 gsl_matrix_free(Product2);
}

