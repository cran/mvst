\name{mcSE}
\alias{mcSE}
\title{MC sampler for a p-variate Skew-Elliptical model.}
\description{MonteCarlo sampler for a p-variate Skew-Elliptical model.}

\usage{mcSE(y, X=NULL, N, Ti, modelType='ST', warmUp=FALSE, control=list())}

\arguments{
 \item{y}{data matrix.}
 \item{X}{design matrix. (default: NULL)}
 \item{N}{number of particles per iteration.}
 \item{Ti}{number of iterations.}
 \item{modelType}{model to be estimated (default: ST). See 'Details'.}
 \item{warmUp}{logical flag (default: FALSE). See 'Deatils'}
 \item{control}{(optional) a list of control parameters. See 'Details'.}
}

\value{The function returns
 \item{particles}{a list with the resampled particles of the last iteration,}
 \item{log.py}{the vector of the estimates of the marginal likelihood (one value per iteration),}
 \item{nResampled}{the vector of the number of resampled particles (one value per iteration).}
If saveParticles is TRUE, the lists of the sampled particles, the importance weights, and the indices of the resampled particles are saved in the folder specified in outFolder, or in the default folder '/Output'. If outFolder doesn't already exists, it will be created.}

\details{Already implemented modelTypes are 'N' (Normal), 'SN' (skew-normal), 'T' (Student T), and 'ST' (skew-t, the default). To estimate a regression model, an 'X' should be added: for example, 'STX' stands for a regression model with ST errors. For these models, the argument parTypes in the control list is overridden. It is however possible to implement other models; in this case, parType is required and should contain the names of the parameters of the model. The argument warmUp allows to run preliminary iterations with a smaller number of particles, in order to speed up the algorithm. The number of these iterations, and the number of particles can be controlled using the Nwu argument in the control list.

To estimate regression models with skewed errors, it is sufficient to specify the argument X, which should contain the design matrix.

The (optional) argument control can provide a list with the following elements
\describe{
 \item{\code{seed}}{ if different from NULL, sets the random seed for replicability purposes.}
 \item{\code{parInfo}}{ data.frame containing the informations about the model parameters. Each row of the data.frame should contain the names of the parameters, the type (u: 'univariate', m: 'multivariate', M: 'matrix-variate', SM: 'symmetric matrix-variate'), the number nCols of elements, or columns, of the parameter and the number nRows of rows (eventually 1). Required if the modelType is not already implemented.}
 \item{\code{propFuncs}}{ named character vector with the names of the functions for the proposal distributions. The names of the elements of propFuncs should be the relevant elements of the set ('z', 'v', 'xi', 'psi', 'G', 'nu'). Custom proposal functions should require at least four arguments: y, the data matrix, X, the covariates (NULL if not relevant), particles, that is the list of current values of the parameters, for each particle, and priorList, a list containing the hyperparameters of the prior distributions. It should return two objects: values (the proposed values for the parameter, for each particle) and log.dq (vector with the N (log-)values of the proposal density).}
 \item{\code{logPriorFunc}}{ name of the file containing the function to compute the value (in logarithms) of the posterior density for all the particles. It should depend on the objects y, particles and priorList, while it should return N values of the posterior density, in logarithms.}
 \item{\code{Nwu}}{ numeric vector with the number of particles for each warm-up iteration. Default value is rep(2000, 3). It is however overridden if warmUp is FALSE.}
 \item{\code{priorList}}{ list of hyperparameters.}
 \item{\code{saveParticles}}{ logical flag (default: FALSE) indicating whether the value of the particles proposed in each iteration should be saved. If TRUE, the folder 'Iterations' is created in the current directory.}
 \item{\code{outFolder}}{ the folder in which the outputs are saved (if saveParticles is TRUE). The default folder is '/Output/Iterations'.}
 \item{\code{verbose}}{ logical flag (default: TRUE). If TRUE, details about the progress of the algorithm are printed.}
 }
}

\references{Parisi, A. and Liseo, B. (2017) "Objective Bayesian Analysis for the Multivariate Skew-t Model" \emph{Statistical Methods & Applications}

Azzalini, A. and Arellano-Valle, R.B. (2013) "Maximum Penalized Likelihood Estimation for Skew-normal and Skew-t Distributions" \emph{J. Statist. Plann. Inference}, \bold{143 (2)}, 419--433.}

\seealso{\code{\link{cmlSE}}, \code{\link{rmvSE}}.}

\examples{
## Generate artificial data
pars = list(xi=c(3,5), psi=c(2,4), G=diag(2), nu=6)
values = rmvSE(n=60, p=2, modelType='ST', theta=pars)
## Estimate a Skew-t model (not run)
# fit = mcSE(y=values$y, N=20000, Ti=3, modelType='ST')
# stats = summary(fit)
# coef(stats)
}
