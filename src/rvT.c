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

//int rvT(double *values, double *propdens, int *n, int *p, int *N, double *yVec, double *nu, double *xiVec, double *psiVec, double *GVec, double *zVec){
void rvT(double *values, double *propdens, int *n, int *p, int *N, double *yVec, double *nu, double *xiVec, double *GVec){
 gsl_matrix *y = gsl_matrix_calloc(*n, *p);
 gsl_matrix *xi = gsl_matrix_calloc(*p, 1);
// gsl_matrix *psi = gsl_matrix_calloc(*p, 1);
 gsl_matrix *G = gsl_matrix_calloc(*p, *p);
 gsl_matrix *invG = gsl_matrix_calloc(*p, *p);
 gsl_permutation *perm = gsl_permutation_alloc(*p);
 gsl_matrix *ymxi = gsl_matrix_calloc(*p, 1);		// y[i,] - xi[,j]
 gsl_matrix *Product1 = gsl_matrix_calloc(1, *p);
 gsl_matrix *Product2 = gsl_matrix_calloc(1, 1);
// gsl_matrix *Product3 = gsl_matrix_calloc(1, 1);
 double Ai[1];
/* double Bi[1];
 double g[1];*/
 double Alphav[1];
/* double sqrtBetav[1];
 double Betav[1];*/
// double vi[1];

// int i, j, k;
 size_t i, j;
 for (j = 0; j < *p; j++){
  for(i = 0; i < *n; i++){
   gsl_matrix_set(y, i, j, yVec[i * *p + j]);
  }
 }

// int jj, k;
 size_t jj, k;
 int s;
 double zi;
 int counter;
 counter = 0;
 for(k = 0; k < *N; k++){	// k = 0;
  Alphav[0] = 0.5 * (nu[k] + *p);
  for (j = 0; j < *p; j++){
   gsl_matrix_set(xi, j, 0, xiVec[k * *p + j]);
//   gsl_matrix_set(psi, j, 0, psiVec[k * *p + j]);
   for (jj = 0; jj < *p; jj++){
//  gsl_matrix_set(G, j, jj, GVec[jj * *p + j]);
    gsl_matrix_set(G, j, jj, GVec[k * *p * *p + jj * *p + j]);
   }
  }
  gsl_linalg_LU_decomp (G, perm, &s);
  gsl_linalg_LU_invert(G, perm, invG);
  for(i = 0; i < *n; i++){
//   zi = zVec[k * *n + i];
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
   Ai[0] = 0.5 * (nu[k] + productA);
/*   Bi[0] = - productB * fabsl(zi);
   g[0] = - 0.5 * Bi[0] * expl(gsl_sf_lngamma(Alphav[0] + 0.5) - gsl_sf_lngamma(Alphav[0]));
   sqrtBetav[0] = (- g[0] + sqrt(g[0] * g[0] + 4 * Ai[0] * Alphav[0] * Alphav[0])) / (2 * Alphav[0]);
   Betav[0] = sqrtBetav[0] * sqrtBetav[0]; */
   GetRNGstate();
   values[counter] = rgamma(Alphav[0], 1/Ai[0]);
   PutRNGstate();
   propdens[k] = propdens[k] + dgamma(values[counter], Alphav[0], 1/Ai[0], 1);
   counter++;
  }
 }

 gsl_matrix_free(y);
 gsl_matrix_free(xi);
// gsl_matrix_free(psi);
 gsl_matrix_free(G);
 gsl_matrix_free(invG);
 gsl_matrix_free(ymxi);
// gsl_matrix_free(perm);
 gsl_matrix_free(Product1);
 gsl_matrix_free(Product2);
// gsl_matrix_free(Product3);
// return 0;
}

