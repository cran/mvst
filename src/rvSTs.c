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

void my_error_handler(const char *reason, const char *file, int line, int err){
 if(err > 0){
  Rprintf("\nAn error occurred in the evaluation of a 1F1 function - fixed.");
 }
}

double fcv_density (double x, void * p) {
 struct fcv_params * params = (struct fcv_params *)p;
 double Ai = (params->Ai);
 double Bi = (params->Bi);
 double CC = (params->CC);
 double logDens;

 logDens = (CC-1) * log(x) - Ai * x - Bi * pow(x, 0.5);
 return exp(logDens);
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

 gsl_integration_workspace_free(work_ptr);

 return result;
}

void rvST(double *values, double *propdens, int *n, int *p, int *N, double *yVec, double *nu, double *xiVec, double *psiVec, double *GVec, double *zVec){
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
 double kpi[1];
 double Alphav[1];
 double Betav[1];
 double logM[1];
 double logm[1];
 double viStar[1];
 double sqrtv[1];
 double prop[1];
 double u[1];
 double accept[1];
 double params[3];
 double z[1];                  // argument of the parabolic cylinder function
 double halfz2[1];
 double productA[1];
 double productB[1];

// int i, j, iN;
 size_t i, j;
 for(j = 0; j < *p; j++){
  for(i = 0; i < *n; i++){
   gsl_matrix_set(y, i, j, yVec[i * *p + j]);
  }
 }

// int jj, iN;
 size_t jj, iN;
 int s;
 double zi;
 int counter;
 counter = 0;
 for(iN = 0; iN < *N; iN++){	// iN = 0;
  for(j = 0; j < *p; j++){
   gsl_matrix_set(xi, j, 0, xiVec[iN * *p + j]);
   gsl_matrix_set(psi, j, 0, psiVec[iN * *p + j]);
   for(jj = 0; jj < *p; jj++){
    gsl_matrix_set(G, j, jj, GVec[iN * *p * *p + jj * *p + j]);
   }
  }
  gsl_linalg_LU_decomp(G, perm, &s);
  gsl_linalg_LU_invert(G, perm, invG);
  for(i = 0; i < *n; i++){
   zi = zVec[iN * *n + i];
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
   productA[0] = gsl_matrix_get(Product2, 0, 0);
   productB[0] = gsl_matrix_get(Product3, 0, 0);
   Ai[0] = 0.5 * (nu[iN] + productA[0]);
   Bi[0] = - productB[0] * fabs(zi);
   CC[0] = 0.5 * (nu[iN] + *p);
   z[0] = Bi[0] / sqrt(2.0 * Ai[0]);
   halfz2[0] = Bi[0] * Bi[0] / (4.0 * Ai[0]);
/*
   gsl_error_handler_t *old_handler = gsl_set_error_handler(&my_error_handler);
   DD[0] = pow(2.0, -CC[0]) * exp(Bi[0] * Bi[0] / (8.0 * Ai[0])) * ((sqrt(M_PI) / gsl_sf_gamma(0.5 * (1 + 2 * CC[0]))) * gsl_sf_hyperg_1F1(CC[0], 0.5, Bi[0] * Bi[0]/(4.0 * Ai[0])) - (sqrt(2.0 * M_PI) * Bi[0] / (sqrt(2.0 * Ai[0]) * gsl_sf_gamma(CC[0]))) * gsl_sf_hyperg_1F1(0.5 * (1.0 + 2.0 * CC[0]), 1.5, Bi[0] * Bi[0] / (4.0 * Ai[0])));
   gsl_set_error_handler(old_handler);
*/
   if(Bi[0] < 0.0){
// First version
//    kpi[0] = exp((-2*CC[0]+1) * M_LN2 - CC[0] * log(Ai[0]) + gsl_sf_lngamma(2.0 * CC[0])) * (M_SQRTPI * gsl_sf_hyperg_1F1(CC[0], 0.5, halfz2[0]) / gsl_sf_gamma(CC[0] + 0.5) - M_SQRT2 * M_SQRTPI * z[0] * gsl_sf_hyperg_1F1(CC[0] + 0.5, 1.5, halfz2[0]) / gsl_sf_gamma(CC[0]));
// Revised version
    kpi[0] = exp((-2*CC[0]+1) * M_LN2 + (- CC[0] + 0.5) * log(Ai[0]) + gsl_sf_lngamma(2.0 * CC[0] - 1.0)) * (M_SQRTPI * gsl_sf_hyperg_1F1(CC[0] - 0.5, 0.5, halfz2[0]) / gsl_sf_gamma(CC[0]) - M_SQRT2 * M_SQRTPI * z[0] * gsl_sf_hyperg_1F1(CC[0], 1.5, halfz2[0]) / gsl_sf_gamma(CC[0] - 0.5));
   } else {
    // Rprintf("\nNumerical integration for particle %d unit %d", iN, i); 
    params[0] = Ai[0];
    params[1] = Bi[0];
    params[2] = CC[0];
    kpi[0] = fcv_integrate(params);
   }
//           Rprintf("\n particle %d unit %d Ai %f Bi %f CC %f DD %f halfz2 %f", iN, i, Ai[0], Bi[0], CC[0], DD[0], halfz2[0]); 
   Alphav[0] = nu[iN] + *p;
   Betav[0] = 0.5 * (Bi[0] + sqrt(Bi[0] * Bi[0] + 8 * Ai[0] * (nu[iN] + *p + 1)));
   viStar[0] = exp(2.0 * (log(Betav[0] - Bi[0]) - M_LN2 - log(Ai[0])));
   logM[0] = M_LN2 + gsl_sf_lngamma(Alphav[0]) - Alphav[0] * log(Betav[0]) - log(kpi[0]) - Ai[0] * viStar[0] - (Bi[0] - Betav[0]) * pow(viStar[0], 0.5);
//                printf ("logm =%f, logM =%f, prop =%f\n", logm[0], logM[0], prop);
//                printf ("logM = %f, viStar = %f, DD = %f, kpi = %f\n", logM[0], viStar[0], DD[0], kpi[0]);
   accept[0] = 0.0;
   GetRNGstate();
   while(accept[0] == 0.0){
    sqrtv[0] = rgamma(Alphav[0], 1.0 / Betav[0]);
    prop[0] = sqrtv[0] * sqrtv[0];
    u[0] = runif(0.0, 1.0);
//      logm[0] = Alphav[0] * log(Betav[0]) - M_LN2 - gsl_sf_lngamma(Alphav[0]) - log(kpi[0]) - Ai[0] * prop[0] - (Bi[0] - Betav[0]) * pow(prop[0], 0.5);
    logm[0] = M_LN2 + gsl_sf_lngamma(Alphav[0]) - Alphav[0] * log(Betav[0]) - log(kpi[0]) - Ai[0] * prop[0] - (Bi[0] - Betav[0]) * pow(prop[0], 0.5);
    if(u[0] < exp(logm[0] - logM[0])) accept[0] = 1.0;
   }
   values[counter] = prop[0];
   PutRNGstate();
   propdens[iN] = propdens[iN] + (CC[0]-1) * log(prop[0]) - Ai[0] * prop[0] - Bi[0] * pow(prop[0], 0.5) - log(kpi[0]);
/*
    values[counter] = 1.0;
    propdens[iN] = R_PosInf;
*/
/*
   if(kpi[0] <= 0.0){
    Rprintf("\nNumerical integration for particle %d unit %d", iN, i); 
    params[0] = Ai[0];
    params[1] = Bi[0];
    params[2] = CC[0];
    kpi[0] = fcv_integrate(params);
   }
*/
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
 gsl_permutation_free(perm);
 gsl_matrix_free(Product1);
 gsl_matrix_free(Product2);
 gsl_matrix_free(Product3);
}

void rvSTX(double *values, double *propdens, int *n, int *p, int *k, int *N, double *yVec, double *XVec, double *nu, double *BVec, double *psiVec, double *GVec, double *zVec){
 gsl_matrix *y = gsl_matrix_calloc(*n, *p);
 gsl_matrix *psi = gsl_matrix_calloc(*p, 1);
 gsl_matrix *G = gsl_matrix_calloc(*p, *p);
 gsl_matrix *invG = gsl_matrix_calloc(*p, *p);
 gsl_permutation *perm = gsl_permutation_alloc(*p);
 gsl_matrix *ymBX = gsl_matrix_calloc(*p, 1);		// y - BX
 gsl_matrix *Product1 = gsl_matrix_calloc(1, *p);
 gsl_matrix *Product2 = gsl_matrix_calloc(1, 1);
 gsl_matrix *Product3 = gsl_matrix_calloc(1, 1);
 double Ai[1];
 double Bi[1];
 double CC[1];
 double kpi[1];
 double Alphav[1];
 double Betav[1];
 double logM[1];
 double logm[1];
 double viStar[1];
 double sqrtv[1];
 double prop[1];
 double u[1];
 double accept[1];
 double params[3];
 double z[1];                  // argument of the parabolic cylinder function
 double halfz2[1];
 double productA[1];
 double productB[1];

 size_t i, j, l;
 for(j = 0; j < *p; j++){
  for(i = 0; i < *n; i++){
   gsl_matrix_set(y, i, j, yVec[i * *p + j]);
  }
 }

 size_t jj, iN;
 int s;
 double zi;
 double BXj;
 int counter;
 counter = 0;
 for(iN = 0; iN < *N; iN++){	// iN = 0;
  for(j = 0; j < *p; j++){
   gsl_matrix_set(psi, j, 0, psiVec[iN * *p + j]);
   for(jj = 0; jj < *p; jj++){
    gsl_matrix_set(G, j, jj, GVec[iN * *p * *p + jj * *p + j]);
   }
  }
  gsl_linalg_LU_decomp(G, perm, &s);
  gsl_linalg_LU_invert(G, perm, invG);
  for(i = 0; i < *n; i++){
   zi = zVec[iN * *n + i];
   for(j = 0; j < *p; j++){
    BXj = 0;
    for (l = 0; l < *k; l++){
     BXj += BVec[iN * (*k * *p) + j * *k + l] * XVec[i * *k + l];
    }
    gsl_matrix_set(ymBX, j, 0, gsl_matrix_get(y, i, j) - BXj);
   }
   gsl_matrix_set_zero(Product1);
   gsl_matrix_set_zero(Product2);
   gsl_matrix_set_zero(Product3);
   gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, ymBX, invG, 0.0, Product1);
   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Product1, ymBX, 0.0, Product2);
   gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Product1, psi, 0.0, Product3);
   productA[0] = gsl_matrix_get(Product2, 0, 0);
   productB[0] = gsl_matrix_get(Product3, 0, 0);
   Ai[0] = 0.5 * (nu[iN] + productA[0]);
   Bi[0] = - productB[0] * fabs(zi);
   CC[0] = 0.5 * (nu[iN] + *p);
   z[0] = Bi[0] / sqrt(2.0 * Ai[0]);
   halfz2[0] = Bi[0] * Bi[0] / (4.0 * Ai[0]);
/*
   gsl_error_handler_t *old_handler = gsl_set_error_handler(&my_error_handler);
   DD[0] = pow(2.0, -CC[0]) * exp(Bi[0] * Bi[0] / (8.0 * Ai[0])) * ((sqrt(M_PI) / gsl_sf_gamma(0.5 * (1 + 2 * CC[0]))) * gsl_sf_hyperg_1F1(CC[0], 0.5, Bi[0] * Bi[0]/(4.0 * Ai[0])) - (sqrt(2.0 * M_PI) * Bi[0] / (sqrt(2.0 * Ai[0]) * gsl_sf_gamma(CC[0]))) * gsl_sf_hyperg_1F1(0.5 * (1.0 + 2.0 * CC[0]), 1.5, Bi[0] * Bi[0] / (4.0 * Ai[0])));
   gsl_set_error_handler(old_handler);
*/
   if(Bi[0] < 0.0){
// First version
    kpi[0] = exp((-2*CC[0]+1) * M_LN2 - CC[0] * log(Ai[0]) + gsl_sf_lngamma(2.0 * CC[0])) * (M_SQRTPI * gsl_sf_hyperg_1F1(CC[0], 0.5, halfz2[0]) / gsl_sf_gamma(CC[0] + 0.5) - M_SQRT2 * M_SQRTPI * z[0] * gsl_sf_hyperg_1F1(CC[0] + 0.5, 1.5, halfz2[0]) / gsl_sf_gamma(CC[0]));
// Revised version
//    kpi[0] = pow(Ai[0], -CC[0]) * (gsl_sf_gamma(CC[0]) * gsl_sf_hyperg_1F1(CC[0], 0.5, halfz2[0]) - gsl_sf_gamma(CC[0] + 0.5) * Bi[0] * gsl_sf_hyperg_1F1(CC[0] + 0.5, 1.5, halfz2[0]) / sqrt(Ai[0]));
//    kpi[0] = pow(Ai[0], -CC[0]) * (gsl_sf_gamma(CC[0]) * gsl_sf_hyperg_1F1(CC[0], 0.5, halfz2[0]) - gsl_sf_gamma(CC[0] + 0.5) * Bi[0] * gsl_sf_hyperg_1F1(CC[0] + 0.5, 1.5, halfz2[0]) / sqrt(Ai[0]));
   } else {
    // Rprintf("\nNumerical integration for particle %d unit %d", iN, i); 
    params[0] = Ai[0];
    params[1] = Bi[0];
    params[2] = CC[0];
    kpi[0] = fcv_integrate(params);
   }
//           Rprintf("\n particle %d unit %d Ai %f Bi %f CC %f DD %f HGFarg %f", iN, i, Ai[0], Bi[0], CC[0], DD[0], HGFarg[0]); 
   /*
   if(HGFarg[0] > 500.0){       // If HGF > 100, delete the particle (vi = 1 for all i's, log.dq = 0 for the particle)
    Rprintf("\n particle %d unit %d Ai %f Bi %f CC %f DD %f HGFarg %f - skipped", iN, i, Ai[0], Bi[0], CC[0], DD[0], HGFarg[0]);
    for(i = 0; i < *n; i++){
     GetRNGstate();
     values[counter] = rgamma(CC[0], 1.0 / Ai[0]);
     PutRNGstate();
    }
    propdens[iN] = R_PosInf;
    counter++;
    break;
   } */
   /*
   if(((Bi[0] > 0) && (HGFarg[0] > 5.0)) || ((Bi[0] < 0) && (HGFarg[0] > 500.0))){  // In the problematic cases, use the numerical approximation
    Rprintf("\n particle %d unit %d Ai %f Bi %f CC %f DD %f HGFarg %f - numeric approx", iN, i, Ai[0], Bi[0], CC[0], DD[0], HGFarg[0]);
    params[0] = Ai[0];
    params[1] = Bi[0];
    params[2] = CC[0];
    kpi[0] = fcv_integrate(params);
   } else {       // If !((HGFarg > 5 & B > 0)), use analytic formulae
//             Rprintf("\n particle %d unit %d Ai %f Bi %f CC %f DD %f HGFarg %f", iN, i, Ai[0], Bi[0], CC[0], DD[0], HGFarg[0]);
    DD[0] = pow(2.0, -CC[0]) * exp(Bi[0] * Bi[0] / (8.0 * Ai[0])) * ((sqrt(M_PI) / gsl_sf_gamma(0.5 * (1 + 2 * CC[0]))) * gsl_sf_hyperg_1F1(CC[0], 0.5, Bi[0] * Bi[0]/(4.0 * Ai[0])) - (sqrt(2.0 * M_PI) * Bi[0] / (sqrt(2.0 * Ai[0]) * gsl_sf_gamma(CC[0]))) * gsl_sf_hyperg_1F1(0.5 * (1.0 + 2.0 * CC[0]), 1.5, Bi[0] * Bi[0] / (4.0 * Ai[0])));
    if(DD[0] > 0){
//         Rprintf("\n particle %d unit %d Ai %f Bi %f CC %f DD %f HGFarg %f", iN, i, Ai[0], Bi[0], CC[0], DD[0], HGFarg[0]);
     kpi[0] = 2.0 * pow(2.0 * Ai[0], -CC[0]) * gsl_sf_gamma(2.0 * CC[0]) * exp(-Bi[0] * Bi[0]/(8.0 * Ai[0])) * DD[0];
    } else {     // If the analytic value of D < 0, use the numerical approximation
         Rprintf("\n particle %d unit %d Ai %f Bi %f CC %f DD %f HGFarg %f - numeric approx (D <= 0)", iN, i, Ai[0], Bi[0], CC[0], DD[0], HGFarg[0]);
     params[0] = Ai[0];
     params[1] = Bi[0];
     params[2] = CC[0];
     kpi[0] = fcv_integrate(params);
    }
   }
   */
   Alphav[0] = nu[iN] + *p;
   Betav[0] = 0.5 * (Bi[0] + sqrt(Bi[0] * Bi[0] + 8 * Ai[0] * (nu[iN] + *p + 1)));
//  if((halfz2[0] > 90) && (halfz2[0] < 110)) Rprintf("\n check1: betav = %f, expArg = %f", Betav[0], 2.0 * (log(Betav[0] - Bi[0]) - M_LN2 - log(Ai[0])));
   viStar[0] = exp(2.0 * (log(Betav[0] - Bi[0]) - M_LN2 - log(Ai[0])));
   logM[0] = M_LN2 + gsl_sf_lngamma(Alphav[0]) - Alphav[0] * log(Betav[0]) - log(kpi[0]) - Ai[0] * viStar[0] - (Bi[0] - Betav[0]) * pow(viStar[0], 0.5);
//                 printf ("logm =%f, logM =%f, prop =%f\n", logm[0], logM[0], prop);
//                 printf ("logM = %f, viStar = %f, DD = %f, kpi = %f\n", logM[0], viStar[0], DD[0], kpi[0]);
   accept[0] = 0.0;
   GetRNGstate();
   while(accept[0] == 0.0){
    sqrtv[0] = rgamma(Alphav[0], 1.0 / Betav[0]);
    prop[0] = sqrtv[0] * sqrtv[0];
    u[0] = runif(0.0, 1.0);
//     logm[0] = Alphav[0] * log(Betav[0]) - M_LN2 - gsl_sf_lngamma(Alphav[0]) - log(kpi[0]) - Ai[0] * prop[0] - (Bi[0] - Betav[0]) * pow(prop[0], 0.5);
    logm[0] = M_LN2 + gsl_sf_lngamma(Alphav[0]) - Alphav[0] * log(Betav[0]) - log(kpi[0]) - Ai[0] * prop[0] - (Bi[0] - Betav[0]) * pow(prop[0], 0.5);
//  if((HGFarg[0] > 90) && (HGFarg[0] < 110)) Rprintf("\n check2: logm = %f, logM = %f", logm[0], logM[0]);
    if(u[0] < exp(logm[0] - logM[0])) accept[0] = 1.0;
//  if((HGFarg[0] > 90) && (HGFarg[0] < 110)) Rprintf("\n check3. expRes = %f", exp(logm[0] - logM[0]));
   }
   values[counter] = prop[0];
   PutRNGstate();
   propdens[iN] = propdens[iN] + (CC[0]-1) * log(prop[0]) - Ai[0] * prop[0] - Bi[0] * pow(prop[0], 0.5) - log(kpi[0]);
   counter++;
  }
 }
//                printf ("logM = %f, viStar = %f, DD = %f, kpi = %f, alphav = %f, betav = %f\n", logM[0], viStar[0], DD[0], kpi[0], Alphav[0], Betav[0]);

 gsl_matrix_free(y);
 gsl_matrix_free(psi);
 gsl_matrix_free(G);
 gsl_matrix_free(invG);
 gsl_matrix_free(ymBX);
 gsl_permutation_free(perm);
 gsl_matrix_free(Product1);
 gsl_matrix_free(Product2);
 gsl_matrix_free(Product3);
}
