\name{genvasy}
\alias{genvasy}

\title{
Asymptotic for groupwise envelope model
}

\description{
It calculates the asymptotic standard error matrix for elements 
in \code{mu} and \code{beta}. 
}

\usage{
genvasy(X, Y, Z, u)
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
  \item{u}{
  A given dimension of the groupwise envelope space.  It should 
  be an interger between \eqn{0} and \eqn{r}.
  }
}

\value{
  \item{asySEmu}{
  The asymptotic standard error matrix for elements in \code{mu}. 
  The multiplication by the reciprocal of square root of \eqn{n} 
  returns actual standard errors.
  }
  \item{asySEbeta}{
  The asymptotic standard error matrix for elements in \code{beta}. 
  The multiplication by the reciprocal of square root of \eqn{n} 
  returns actual standard errors.
  }
}

\references{
  Park, Y., Su, Z. and Zhu, H. (2016)  Groupwise envelope models 
  for Imaging Genetic Analysis. 
}

\examples{
data(fiberpaper)
Y <- as.matrix(data[1:4])
XX <- as.numeric(data[,7] > mean(data[,7]))
X <- as.vector(XX)
Z <- as.matrix(data[,5:6])
u <- 3
asyse <- genvasy(X, Y, Z, u)
asyse$asySEmu
asyse$asySEbeta[[1]]
asyse$asySEbeta[[2]]
}

