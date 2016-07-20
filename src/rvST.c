#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_integration.h>

struct fcv_params {double Ai; double Bi; double CC;};

double fcv_density (double x, void * p) {
 struct fcv_params * params = (struct fcv_params *)p;
 double Ai = (params->Ai);
 double Bi = (params->Bi);
 double CC = (params->CC);
 double logDens;

 logDens = (CC-1) * logl(x) - Ai * x - Bi * pow(x, 0.5);
 return expl(logDens);
}

double fcv_integrate (double *pars) {
 gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (1000);

 double lower_limit = 0.0;
 double abs_error = 1.0e-8;
 double rel_error = 1.0e-8;
 double result;
 double error;

 gsl_function fcv;
 struct fcv_params params;
 params.Ai = pars[0];
 params.Bi = pars[1];
 params.CC = pars[2];

 fcv.function = &fcv_density;
 fcv.params = &params;

 gsl_integration_qagiu (&fcv, lower_limit,
  abs_error, rel_error, 1000, work_ptr, &result,
  &error);
// printf ("result = % .18f\n", result);

 return result;
}

int rvST(double *values, double *propdens, int *n, int *p, int *N, double *yVec, double *nu, double *xiVec, double *psiVec, double *GVec, double *zVec){
 gsl_matrix *y = gsl_matrix_calloc(*n, *p);
 gsl_matrix *xi = gsl_matrix_calloc(*p, 1);
 gsl_matrix *psi = gsl_matrix_calloc(*p, 1);
 gsl_matrix *G = gsl_matrix_calloc(*p, *p);
 gsl_matrix *invG = gsl_matrix_calloc(*p, *p);
 gsl_permutation *perm = gsl_permutation_alloc(*p);
 gsl_matrix *ymxi = gsl_matrix_calloc(*p, 1);		// y[i,] - xi[,j]
 gsl_matrix *Product1 = gsl_matrix_calloc(1, *p);
 gsl_matrix *Product2 = gsl_matrix_calloc(1, 1);
 gsl_matrix *Product3 = gsl_matrix_calloc(1, 1);
 double Ai[1];
 double Bi[1];
 double CC[1];
 double DD[1];
// double normConstFC[1];
 double kpi[1];
// double logNormConstProp[1];
// double g[1];
 double Alphav[1];
// double sqrtBetav[1];
 double Betav[1];
 double logM[1];
 double logm[1];
 double viStar[1];
// double vi[1];
 double sqrtv[1];
 double prop[1];
 double u[1];
 double accept[1];
 
 void my_error_handler(const char *reason, const char *file, int line, int err){
  if(err > 0){
   Rprintf("\nAn error occurred in the evaluation of a 1F1 function - fixed.");
  }
 }

// int i, j, k;
 size_t i, j;
 for(j = 0; j < *p; j++){
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
  for(j = 0; j < *p; j++){
   gsl_matrix_set(xi, j, 0, xiVec[k * *p + j]);
   gsl_matrix_set(psi, j, 0, psiVec[k * *p + j]);
   for(jj = 0; jj < *p; jj++){
//  gsl_matrix_set(G, j, jj, GVec[jj * *p + j]);
    gsl_matrix_set(G, j, jj, GVec[k * *p * *p + jj * *p + j]);
   }
  }
  gsl_linalg_LU_decomp(G, perm, &s);
  gsl_linalg_LU_invert(G, perm, invG);
  for(i = 0; i < *n; i++){
   zi = zVec[k * *n + i];
//  gsl_matrix_set_zero(ymxi);
//  ymxi = yi;
//  gsl_matrix_sub(ymxi, xi);
   for(j = 0; j < *p; j++){
    gsl_matrix_set(ymxi, j, 0, gsl_matrix_get(y, i, j) - gsl_matrix_get(xi, j, 0));
   }
   gsl_matrix_set_zero(Product1);
   gsl_matrix_set_zero(Product2);
   gsl_matrix_set_zero(Product3);
   gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, ymxi, invG, 0.0, Product1);
   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Product1, ymxi, 0.0, Product2);
   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Product1, psi, 0.0, Product3);
   double productA = gsl_matrix_get(Product2, 0, 0);
   double productB = gsl_matrix_get(Product3, 0, 0);
   Ai[0] = 0.5 * (nu[k] + productA);
   Bi[0] = - productB * fabsl(zi);
   CC[0] = 0.5 * (nu[k] + *p);
   gsl_error_handler_t *old_handler = gsl_set_error_handler(&my_error_handler);
   DD[0] = pow(2.0, -CC[0]) * expl(Bi[0] * Bi[0] / (8.0 * Ai[0])) * ((sqrt(M_PI) / gsl_sf_gamma(0.5 * (1 + 2 * CC[0]))) * gsl_sf_hyperg_1F1(CC[0], 0.5, Bi[0] * Bi[0]/(4.0 * Ai[0])) - (sqrt(2.0 * M_PI) * Bi[0] / (sqrt(2.0 * Ai[0]) * gsl_sf_gamma(CC[0]))) * gsl_sf_hyperg_1F1(0.5 * (1.0 + 2.0 * CC[0]), 1.5, Bi[0] * Bi[0] / (4.0 * Ai[0])));
   gsl_set_error_handler(old_handler);
   if(isnan(DD[0]) == 0){
    kpi[0] = 2.0 * pow(2.0 * Ai[0], -CC[0]) * gsl_sf_gamma(2.0 * CC[0]) * expl(-Bi[0] * Bi[0]/(8.0 * Ai[0])) * DD[0];
   } else {
    Rprintf("\nNumerical integration of the 1F1 function for particle %d unit %d", k, i); 
    double params[3];
    params[0] = Ai[0];
    params[1] = Bi[0];
    params[2] = CC[0];
    kpi[0] = fcv_integrate(params);
   }
   if(kpi[0] <= 0.0){
/* values[counter] = 100;
    propdens[k] = R_NegInf;
    printf("\n RejSamp.C: A negative value for the normalizing constant has been replaced with a 100\n");
*/
    Rprintf("\nNumerical integration for particle %d unit %d", k, i); 
    double params[3];
    params[0] = Ai[0];
    params[1] = Bi[0];
    params[2] = CC[0];
    kpi[0] = fcv_integrate(params);
   }
   Alphav[0] = nu[k] + *p;
   Betav[0] = 0.5 * (Bi[0] + sqrt(Bi[0] * Bi[0] + 8 * Ai[0] * (nu[k] + *p + 1)));
   viStar[0] = expl(2.0 * (logl(Betav[0] - Bi[0]) - logl(2.0) - logl(Ai[0])));
   logM[0] = logl(2.0) + gsl_sf_lngamma(Alphav[0]) - Alphav[0] * logl(Betav[0]) - logl(kpi[0]) - Ai[0] * viStar[0] - (Bi[0] - Betav[0]) * pow(viStar[0], 0.5);

//                printf ("logm =%f, logM =%f, prop =%f\n", logm[0], logM[0], prop);
//                printf ("logM = %f, viStar = %f, DD = %f, kpi = %f\n", logM[0], viStar[0], DD[0], kpi[0]);

   accept[0] = 0.0;
   GetRNGstate();
   while(accept[0] == 0.0){
    sqrtv[0] = rgamma(Alphav[0], 1.0 / Betav[0]);
    prop[0] = sqrtv[0] * sqrtv[0];
    u[0] = runif(0.0, 1.0);
//     logm[0] = Alphav[0] * logl(Betav[0]) - logl(2.0) - gsl_sf_lngamma(Alphav[0]) - logl(kpi[0]) - Ai[0] * prop[0] - (Bi[0] - Betav[0]) * pow(prop[0], 0.5);
    logm[0] = logl(2.0) + gsl_sf_lngamma(Alphav[0]) - Alphav[0] * logl(Betav[0]) - logl(kpi[0]) - Ai[0] * prop[0] - (Bi[0] - Betav[0]) * pow(prop[0], 0.5);
    if(u[0] < expl(logm[0] - logM[0])) accept[0] = 1.0;
   }
   values[counter] = prop[0];
   PutRNGstate();
   propdens[k] = propdens[k] + (CC[0]-1) * logl(prop[0]) - Ai[0] * prop[0] - Bi[0] * pow(prop[0], 0.5) - logl(kpi[0]);
   counter++;
  }
 }
//                printf ("logM = %f, viStar = %f, DD = %f, kpi = %f, alphav = %f, betav = %f\n", logM[0], viStar[0], DD[0], kpi[0], Alphav[0], Betav[0]);

 gsl_matrix_free(y);
 gsl_matrix_free(xi);
 gsl_matrix_free(psi);
 gsl_matrix_free(G);
 gsl_matrix_free(invG);
 gsl_matrix_free(ymxi);
// gsl_matrix_free(perm);
 gsl_matrix_free(Product1);
 gsl_matrix_free(Product2);
 gsl_matrix_free(Product3);
 return 0;
}

