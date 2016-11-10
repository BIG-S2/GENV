\name{genvasy}
\alias{genvasy}

\title{
Asymptotic for groupwise envelope model
}

\description{
It calculates the asymptotic standard error matrix for elements 
in \eqn{\bm{\mu}} and \eqn{\bm{\beta}}. 
}

\usage{
genvasy(X, Y, Z, u)
}

\arguments{
  \item{X}{
  Group indicator vector of length \eqn{n}, where \eqn{n} 
  denotes the number of observations.
  }
  \item{Y}{
  An \eqn{n} by \eqn{r} matrix of multivariate responses, where 
  \eqn{r} is the number of responses.
  }
  \item{Z}{
  An \eqn{n} by \eqn{p} matrix of predictors, where \eqn{p} is the 
  number of predictors.
  }
  \item{u}{
  A given dimension of the groupwise envelope space.
  }
}

\value{
  \item{asySEmu}{
  The asymptotic standard error matrix for elements in \eqn{\bm{\mu}}. 
  The multiplication by \eqn{1/\sqrt{n}} returns actual standard errors.
  }
  \item{asySEbeta}{
  The asymptotic standard error matrix for elements in \eqn{\bm{\beta}}. 
  The multiplication by \eqn{1/\sqrt{n}} returns actual standard errors.
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

