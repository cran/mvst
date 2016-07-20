#include <R.h>
#include <Rmath.h>

void leftTruncNorm(double *mu, double *sigma2, double *x){
 int check1, check2;
 double alphaStar, u, muMinus, z;
 muMinus = -*mu/sqrt(*sigma2);
 if (muMinus <= 0.0){
  check1 = FALSE;
  while(check1 == FALSE){
   GetRNGstate();
   z = rnorm(0.0,1.0);
   PutRNGstate();
   check1 = (z > muMinus);
  }
 } else {
  alphaStar = 0.5 * (muMinus + sqrt(muMinus * muMinus + 4.0));
  check2 = FALSE;
  while(check2 == FALSE){
   GetRNGstate();
   z = muMinus + rexp(1/alphaStar);
   PutRNGstate();
   GetRNGstate();
   u = runif(0.0,1.0);
   PutRNGstate();
   check2 = (u <= exp(-0.5*(z-alphaStar) * (z-alphaStar)));
  }
 }
 *x = *mu + z * sqrt(*sigma2);
}

