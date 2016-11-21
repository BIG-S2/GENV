\name{u_genv}
\alias{u_genv}

\title{
Dimension selection of the groupwise envelope subspace using 
Bayesian information criterion
}

\description{
It selects the dimension of the groupwise envelope subspace 
using Bayesian information criterion (BIC).
}

\usage{
u_genv(X, Y, Z)
}

\arguments{
  \item{X}{
  A group indicator vector of length \eqn{n}, where \eqn{n} 
  denotes the number of observations.
  }
  \item{Y}{
  An \eqn{n} by \eqn{r} matrix of multivariate responses, where 
  \eqn{r} is the number of responses.
  }
  \item{Z}{
  An \eqn{n} by \eqn{p} matrix of predictors, where \eqn{p} is the 
  number of predictors with \eqn{p < n}.
  }
}

\value{
  \item{u.bic}{
  Dimension of the groupwise envelope.  An interger between 
  \eqn{0} and \eqn{r}.
  }
  \item{bic.seq}{
  A sequence of BIC values for all possible dimensions.
  }
}

\references{
  Park, Y., Su, Z. and Zhu, H. (2016)  Groupwise envelope models 
  for Imaging Genetic Analysis.
}

\examples{
# fiberpaper data
data(fiberpaper)
Y <- as.matrix(data[1:4])
XX <- as.numeric(data[,7] > mean(data[,7]))
X <- as.vector(XX)
Z <- as.matrix(data[,5:6])
udim <- u_genv(X, Y, Z)$u.bic

# simulation data
data(simul)
udim <- u_genv(X, Y, Z)$u.bic
}

