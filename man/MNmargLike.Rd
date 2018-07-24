\name{MNmargLike}

\alias{MNmargLike}

\title{Marginal Likelihood for the Multivariate Normal Model.}

\description{This function computes the exact marginal likelihood for Normally distributed data, under the default priors.}

\usage{MNmargLike(y, X=NULL, LOG=FALSE)}

\arguments{
 \item{y}{data matrix.}
 \item{X}{(optional) a design matrix.}
 \item{LOG}{logical; if TRUE, the log-marginal likelihood is returned.}
}

\value{A scalar representing the marginal likelihood of a (multivariate) Normal model under the default priors for data y. If the design matrix X is provided, the function returns the marginal likelihood of a (multivariate) regression model with Normally distributed errors.}

\references{Liseo B, Parisi A (2013). Bayesian Inference for the Multivariate Skew-Normal Model: A Population Monte Carlo approach. Comput. Statist. Data Anal., 63, 125-138. ISSN 0167-9473. doi:10.1016/j.csda.2013.02.007.}

\seealso{\code{\link{rmvSE}}, \code{\link{dmvSE}}.}

\examples{
# Generate Normally distributed data
require(mvtnorm)
y = rmvnorm(100, rep(2,2), diag(2))
# Marginal likelihood (exact value)
MNmargLike(y, X=NULL, LOG=TRUE)
}
