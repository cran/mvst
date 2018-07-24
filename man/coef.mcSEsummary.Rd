\name{coef.mcSEsummary}

\alias{coef.mcSEsummary}

\title{Extract mcSE Model Coefficients.}

\description{The point estimates for the model parameters are obtained from mcSE summary objects.}

\usage{\method{coef}{mcSEsummary}(object, ...)}

\arguments{
 \item{object}{an object of class mcSEsummary.}
 \item{...}{further arguments passed to or from other methods.}
}

\value{A list containing the point estimates for the estimated model.}

\seealso{\code{\link{mcSE}}, \code{\link{summary.mcSE}}.}

\examples{
# Generate ST-distributed data (including the value of the latent variables)
pars = list(xi=c(2,2), G=diag(2), psi=c(0.3,0.5), nu=5)
values = rmvSE(n=100, p=2, modelType='ST', theta=pars)
# Estimate a Skew-t model (not run)
# fit = mcSE(y=values$y, X=NULL, N=20000, Ti=3, modelType='ST', warmUp=TRUE)
# stats = summary(fit)
# coef(stats)
}
