\name{genv}
\alias{genv}

\title{
Fitting the groupwise envelope model
}

\description{
It fits the groupwise envelope model for multivariate linear
regression.  The groupwise envelope model allows for both 
distinct regression coefficients and distinct error structures 
for different groups.  When the dimension of the envelope is 
between \eqn{1} and \eqn{r-1}, we implement the method in 
Park et al. (2016).  When the dimension is \eqn{r}, the envelope 
model degenerates to the standard multivariate linear regression.  
When the dimension is \eqn{0}, it means that \code{X} and \code{Y} 
are uncorrelated and the fitting is different.  In addition, 
the groupwise envelope model degenerates to the envelope model 
in Cook et al. (2010) if there is a group.
}

\usage{
genv(X, Y, Z, u)
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
  \item{mu}{
  The estimator of group mean.  A \eqn{r} by \eqn{L} matrix  
  whose ith column contains the mean for the group.
  }
  \item{beta}{
  A list of \eqn{r} by \eqn{p} matrices for the estimator of
  regression coefficients.  \code{beta}[[i]] indicates the 
  estimate of regression coefficient for the ith group.
  }
  \item{Gamma}{
  An \eqn{r} by \eqn{u} matrix for the orthogonal basis of the 
  groupwise envelope subspace.
  }
  \item{Gamma0}{
  An \eqn{r} by \eqn{(r-u)} matrix for the orthogonal basis 
  of the complement of the groupwise envelope subspace.
  }
  \item{eta}{
  A \eqn{u} by \eqn{p} matrix for the coordinates of \code{beta}  
  with respect to \code{Gamma}.
  }
  \item{Omega}{
  A \eqn{u} by \eqn{u} matrix for the coordinates of \code{Sigma}  
  with respect to \code{Gamma}.
  }
  \item{Omega0}{
  A \eqn{(r-u)} by \eqn{(r-u)} matrix for the coordinates of \code{Sigma}  
  with respect to \code{Gamma0}.
  }
  \item{Sigma}{
  A list of the estimator of error covariance matrix. \code{Sigma}[[i]] 
  contains the estimated covariance matrix for the ith group.
  }
  \item{loglik}{
  The maximized log likelihood function.
  }
  \item{Yfit}{
  A \eqn{n} by \eqn{r} matrix for fitted response.
  }
  \item{ng}{
  A \eqn{L} by \eqn{1} vector of the number of observations in 
  each group.
  }
  \item{paranum}{
  The number of parameter.
  }
}

\references{
  Cook, R. D., Li, B. and Chiaromonte, F. (2010) Envelope models 
  for parsimonious and efficient multivariate linear regression 
  (with discussion). \emph{Statistica Sinica} \eqn{\mathbf{20}}{%
  20}, 927-1010.

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

# simulation data
data(simul)
udim <- u_genv(X, Y, Z)$u.bic
res <- genv(X, Y, Z, u=udim)
}
