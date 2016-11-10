\name{cv_genv}
\alias{cv_genv}

\title{
Calculation of prediction error for a given dimension of 
groupwise envelope subspace
}

\description{
It implements \code{m}-fold cross validation to calculate the prediction 
error for a given dimension.  This function randomly splits the 
data.  
}

\usage{
cv_genv(X, Y, Z, u, m, nperm)
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
  \item{m}{
  A positive integer that is used to indicate \code{m}-fold cross 
  validation.
  }
  \item{nperm}{
  A positive integer indicating how many randomly splits to 
  perform \code{m}-fold cross validation.
  }
}

\value{
  Prediction errors for a given dimension of the groupwise 
  envelope subspace.
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
PE <- cv_genv(X, Y, Z, 3, 5, 50)

# simulation data
data(simul)
PE <- cv_genv(X, Y, Z, 3, 5, 50)
}
