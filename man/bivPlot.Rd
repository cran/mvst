\name{bivPlot}

\alias{bivPlot}

\title{Marginal and joint plots for bivariate data.}

\description{Scatterplot and marginal histograms for bivariate data. If theta is provided, the joint and marginal densities will be superimposed.}

\usage{bivPlot(y, modelType=NULL, theta=NULL)}

\arguments{
 \item{y}{Data matrix or vector.}
 \item{modelType}{(optional) model specification. Already implemented modelTypes are 'N' (Normal), 'SN' (skew-normal), 'T' (Student T), and 'ST' (skew-t).}
 \item{theta}{(optional) list of parameters. The list should contain elements named 'xi' (a numeric vector) or 'B' (a kxp matrix), 'G' (a pxp matrix), 'psi' (a numeric vetor, optional) and 'nu' (a scalar, optional). See Details.}
}

\value{The function draws a plot for bivariate data.}

\examples{
# Define the parameters' list
pars = list(xi=c(5,2), G=diag(2), psi=rep(1,2), nu=4)
# Generate data
values = rmvSE(200, 2, NULL, 'ST', theta=pars)
y = values$y
# Draw the data points.
bivPlot(y)
# Draw the data points and the density function.
bivPlot(y, modelType='ST', theta=pars)
}
