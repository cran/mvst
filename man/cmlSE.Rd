\name{cmlSE}

\alias{cmlSE}

\title{CML for the parameters of a p-variate Skew-Elliptical model.}

\description{Complete Maximum Likelihood for the parameters of a p-variate Skew-Elliptical model.}

\usage{cmlSE(modelType, y, z=NULL, v=NULL)}

\arguments{
 \item{modelType}{character indicating the model. Possible values are "N", "SN", "T", "ST"}
 \item{y}{data matrix.}
 \item{z}{numeric vector.}
 \item{v}{numeric vector.}
}

\value{Given the value of the latent variables z and v, the function returns a vector containing the estimates for the required model.}

\seealso{\code{\link{mcSE}}, \code{\link{rmvSE}}.}

\examples{
## Generate artificial data
theta = list(xi=c(3,5), psi=c(2,4), G=diag(2), nu=6)
n = 1000
p = 2
X = rmvSE(n, p, 'ST', theta, seed=NULL)
y = X$y
latentVars = list(z=X$z, v=X$v)
## CML estimates for theta.
thetaHat = cmlSE('ST', y, z=X$z, v=X$v)
}
