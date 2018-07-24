\name{cmlSE}

\alias{cmlSE}

\title{CML for the parameters of a p-variate Skew-Elliptical model.}

\description{Complete Maximum Likelihood for the parameters of a p-variate Skew-Elliptical model.}

\usage{cmlSE(modelType, y, z=NULL, v=NULL, X=NULL)}

\arguments{
 \item{modelType}{character indicating the model. Possible values are "N", "SN", "T", "ST"}
 \item{y}{data matrix.}
 \item{z}{(optional) numeric vector.}
 \item{v}{(optional) numeric vector.}
 \item{X}{(optional) design matrix.}
}

\value{Given the value of the latent variables z and v, the function returns a list containing the estimates for the required model. Where available, a design matrix with the value of the covariates can be provided. In this case, the parameters of a regression model with skewed errors are estimated.}

\references{Parisi, A. and Liseo, B. (2017) "Objective Bayesian Analysis for the Multivariate Skew-t Model" \emph{Statistical Methods & Applications}}

\seealso{\code{\link{mcSE}}, \code{\link{rmvSE}}.}

\examples{
## Generate artificial data
pars = list(xi=c(3,5), psi=c(2,4), G=diag(2), nu=6)
values = rmvSE(n=20, p=2, modelType='ST', theta=pars)
## CML estimates for pars
thetaHat = cmlSE(modelType='ST', y=values$y, z=values$z, v=values$v)
}
