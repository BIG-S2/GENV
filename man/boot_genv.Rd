\name{boot_genv}
\alias{boot_genv}

\title{
Bootstap standard error for groupwise envelope model
}

\description{
It computes bootstrap standard error for the regression 
coefficients in the groupwise envelope model by bootstapping 
the residuals.
}

\usage{
boot_genv(X, Y, Z, u, B)
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
  An \eqn{n} by \eqn{p} matrix of predictors, where \eqn{p} is  
  the number of predictors.
  }
  \item{u}{
  A given dimension of the groupwise envelope space.
  }
  \item{B}{
  A positive integer of bootstrap samples.
  }
}

\value{
  \item{bootsebeta}{
  List of \eqn{r} by \eqn{p} matrices.  \code{bootsebeta}[[i]]  
  denotes the standrd error for elements in \code{beta} for  
  the ith group computed by bootstrap.
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
u <- 3
B <- 200
res <- boot_genv(X, Y, Z, u, B)
res$bootsebeta[[1]]
res$bootsebeta[[2]]

# simulation data
data(simul)
u <- 1
B <- 200
res <- boot_genv(X, Y, Z, u, B)
res$bootsebeta[[1]]
res$bootsebeta[[2]]
}
