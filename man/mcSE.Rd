\name{mcSE}
\alias{mcSE}
\title{MC sampler for a p-variate Skew-Elliptical model.}
\description{MonteCarlo sampler for a p-variate Skew-Elliptical model.}

\usage{mcSE(y, N, Ti, modelType='ST', warmUp=FALSE, control=list())}

\arguments{
 \item{y}{data matrix.}
 \item{N}{number of particles per iteration.}
 \item{Ti}{number of iterations.}
 \item{modelType}{model to be estimated (default: skew-t). See 'Details'.}
 \item{warmUp}{logical flag (default: FALSE). See 'Deatils'}
 \item{control}{(optional) a list of control parameters. See 'Details'.}
}

\value{The function returns
 \item{particles}{a list with the resampled particles of the last iteration,}
 \item{log.py}{the vector of the estimates of the marginal likelihood (one value per iteration),}
 \item{nResampled}{the vector of the number of resampled particles (one value per iteration).}
If saveParticles is TRUE, the lists of the sampled particles, the importance weights, and the indices of the resampled particles are saved in the folder specified in outFolder, or in the default folder '/Output'. If outFolder doesn't already exists, it will be created.}

\details{Already implemented modelTypes are 'N' (Normal), 'SN' (skew-normal), 'T' (Student T), and 'ST' (skew-t, the default). For these models, the argument parTypes in the control list is overridden. It is however possible to implement other models; in this case, parType is required and should contain the names of the parameters of the model. The argument warmUp allows to run preliminary iterations with a smaller number of particles, in order to speed up the algorithm. The number of these iterations, and the number of particles can be controlled using the Nwu argument in the control list.

The (optional) argument control can provide a list with the following elements
\itemize{
 \item{\code{seed}}{ if different from NULL, sets the random seed for replicability purposes.}
 \item{\code{simFlag}}{ logical flag (default: FALSE) indicating whether data are simulated or not. In the former case, real values of parameters and of latent variables should be in the files 'pars.txt', 'otherPars.txt', 'z.txt' and 'v.txt', in the directory './Samples'.}
 \item{\code{parTypes}}{ character vector with the names of the parameters. Required if the modelType is not already implemented.}
 \item{\code{propFuncs}}{ named character vector with the names of the functions for the proposal distributions. The names of the elements of propFuncs should be the relevant elements of the set ('z', 'v', 'xi', 'psi', 'G', 'nu'). Custom proposal functions should require at least two arguments: y, the data matrix and particles, the list of current values of the parameters, for each particle. It should return two objects: values (the proposed values for the parameter, for each particle) and log.dq (vector with the N (log-)values of the proposal density).}
 \item{\code{logPriorFunc}}{ name of the file containing the function to compute the log-value of the posterior density for all the particles. It should depend on the objects y, particles and priorList, while it should return N values of the posterior density, in logarithms.}
 \item{\code{Nwu}}{ (optional) numeric vector with the number of particles for each warm-up iteration. Default value is rep(2000, 3). It is however overridden if warmUp is FALSE.}
 \item{\code{priorList}}{ (optional) list of hyperparameters.}
 \item{\code{saveParticles}}{ logical flag (default: FALSE) indicating whether the value of the particles proposed in each iteration should be saved. If TRUE, the folder 'Iterations' is created in the current directory.}
 \item{\code{outFolder}}{ the folder in which the outputs are saved (if saveParticles is TRUE). The default folder is '/Output/Iterations'.}
 \item{\code{verbose}}{ logical flag (default: TRUE). If TRUE, details about the progress of the algorithm are printed.}
 }
}

\seealso{
  \code{\link{cmlSE}}, \code{\link{rmvSE}}.
}

\examples{
## Generate artificial data
theta = list(xi=c(3,5), psi=c(2,4), G=diag(2), nu=6)
n = 60
p = 2
X = rmvSE(n, p, 'ST', theta, seed=NULL)
y = X$y
## Estimate a skew-t model
fit = mcSE(y, 100, 3, modelType='ST')
summary(fit)
}
