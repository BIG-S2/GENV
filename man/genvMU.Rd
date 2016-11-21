\name{genvMU}
\alias{genvMU}

\title{
Estimate of orthogonal basis of the groupwise envelope subspace
}

\description{
It gives estimator of \code{Gamma} for groupwise envelope model.
}

\usage{
genvMU(M, U, MU, u, n, ng, L)
}

\arguments{
  \item{M}{
  A matrix M for the non-Grassmann manifold optimization problem in 
  Cook et al. (2016)
  }
  \item{U}{
  A matrix U for the non-Grassmann manifold optimization problem in 
  Cook et al. (2016)
  }
  \item{MU}{
  Sum of matrix M and U.
  }
  \item{u}{
  A given dimension of the groupwise envelope space.  It should 
  be an interger between \eqn{0} and \eqn{r}.
  }
  \item{n}{
  The number of observations.
  }
  \item{ng}{
  A \eqn{L} by \eqn{1} vector of the number of observations in 
  each group.
  }
  \item{L}{
  The number of groups.  
  }
}

\value{
  \item{Gamma}{
  An \eqn{r} by \eqn{u} matrix for the orthogonal basis of the  
  groupwise envelope subspace.
  }
  \item{Gamma0}{
  An \eqn{r} by \eqn{(r-u)} matrix for the orthogonal basis of 
  the complement of the groupwise envelope subspace.
  }
}

\references{
  Cook, R. D., Forzani. L. and Su, Z. (2016) A note on fast 
  envelope estimation. \emph{Journal of Multivariate Analysis}. 
  \eqn{\mathbf{150}}{%
  150}, 42-54. 
  
  Park, Y., Su, Z. and Zhu, H. (2016)  Groupwise envelope 
  models for Imaging Genetic Analysis.
}


