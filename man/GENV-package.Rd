\name{GENV-package}
\alias{GENV-package}
\alias{GENV}
\docType{package}

\title{
Groupwise Envelope model
}

\description{
This is a package for fitting the groupwise envelope model.
}

\details{
\tabular{ll}{
Package: \tab GENV\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-10-25\cr
License: \tab GPL-2\cr
}

Main function \code{genv} fits the groupwise envelope 
model to the responses and predictors, using the 
maximum likelihood estimation.  Function \code{genvasy} provides
the asymptotic standard error matrix for estimators. 
We can find the dimension of the envelope subspase from
\code{u_genv}.  The prediction error for a given 
dimension is obtained from \code{cv_genv}.  
}

\author{
  Yeonhee Park, Zhihua Su and Hongtu Zhu.
  
  Maintainer: Yeonhee Park <YPark3@mdanderson.org>
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
res <- genv(X, Y, Z, u=udim)

asyse <- genvasy(X, Y, Z, udim)
asyse$asySEmu
asyse$asySEbeta[[1]]
asyse$asySEbeta[[2]]

# five-fold cross validation to get prediction error for
# the groupwise envelope model with the selection dimension 
# by BIC:
PE <- cv_genv(X, Y, Z, udim, 5, 50)
}
