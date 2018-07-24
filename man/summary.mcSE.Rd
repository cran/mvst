\name{summary.mcSE}

\alias{summary.mcSE}

\title{Summary function for mcSE objects.}

\description{summary method for class "mcSE".}

\usage{\method{summary}{mcSE}(object, ...)}

\arguments{
 \item{object}{an object of class mcSE.}
 \item{...}{further arguments passed to or from other methods.
}
}

\value{A numeric vector with n values of the density function, one for each row in y.}

\references{Parisi A, Liseo B (2017). Objective Bayesian Analysis for the Multivariate Skew-t Model. Statistical Methods & Applications.}

\seealso{\code{\link{mcSE}}, \code{\link{coef.mcSEsummary}}.}

\examples{
# Generate Normally distributed data
pars = list(xi=c(2,2), G=diag(2), psi=c(0.3,0.5), nu=5)
values = rmvSE(n=100, p=2, modelType='N', theta=pars)
# Estimate a Skew-t model (not run)
# fit = mcSE(y=values$y, X=NULL, N=20000, Ti=3, modelType='ST', warmUp=FALSE)
# summary(fit)
}
