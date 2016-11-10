\name{genvMU}
\alias{genvMU}

\title{
Estimate of orthogonal basis of the groupwise envelope subspace
}

\description{
It gives estimator of \eqn{\bm{\Gamma}} for groupwise envelope model.
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
  A given dimension of the envelope subspace.
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

\details{
For each group \eqn{l}, let 
\eqn{\hat{\bm{\Sigma}}_{\code{res}, (l)} = (1/n_{(l)}) \code{Y}_{(l)c}^{\top} \mathbf{Q}_{\code{X}_{(l)}} \code{Y}_{(l)c}} 
and 
\eqn{\hat{\bm{\Sigma}}_{\mathbf{Y}} = (1/n) \sum_{l=1}^L \code{Y}_{(l)c}^{\top}\code{Y}_{(l)c}}, 
where \eqn{\code{X}_{(l)} \in \R^{n_{(l)} \times p}} is the centered 
data matrix for \eqn{\code{X}}, \eqn{\code{Y}_{(l)c} \in \R^{n_{(l)} \times r}} 
is the centerd data matrix for \eqn{\code{Y}} and 
\eqn{\mathbf{Q}_{\mathbf{Y}} = \mathbf{I}-\mathbf{P}_{\mathbf{Y}}} 
with \eqn{\mathbf{P}_{\mathbf{Y}}} denoting the projection 
matrix onto the spanned space by \eqn{\mathbf{Y}}.
The envelope subspace is estimated by solving the spanned space 
by \eqn{\bm{\Gamma}} minimizing 
\deqn{
\sum_{l=1}^L \frac{n_{(l)}}{n} \log|\bm{\Gamma}^{\top} \hat{\bm{\Sigma}}_{\code{res}, (l)} \bm{\Gamma}|
+\log|\bm{\Gamma}^{\top} \hat{\bm{\Sigma}}_{\mathbf{Y}}^{-1} \bm{\Gamma}|.
} 
This is the optimization problem over the Grassmann manifold 
of \eqn{r \times u} matrices.  It can be converted to a 
non-Grassmann manifold optimization problems, which is faster 
and more reliable (Cook et al. 2016). \cr

The objective function in Grassmann manifold optimization is 
written as
\deqn{
\sum_{l=1}^L \bigg\{ \frac{n_{(l)}}{n} \log|\bm{\Gamma}^{\top} M_{(l)} \bm{\Gamma}|
+ (1/L) \log|\bm{\Gamma}^{\top} (M_{(l)} + U_{(l)})^{-1} \bm{\Gamma}| \bigg\},
}
where \eqn{M_{(l)} = \hat{\bm{\Sigma}}_{\code{res}, (l)}} and 
\eqn{U_{(l)} = \hat{\bm{\Sigma}}_{\mathbf{Y}} - \hat{\bm{\Sigma}}_{\code{res}, (l)}}.  
Then, we can follow the way to convert the Grassmann manifold 
optimization to non-Grassmann manifold optimization problem and 
apply the algorithm to minimize the objective function like in 
Cook et al. (2016).
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
  \eqn{\mathbf{150}}, 42-54. 
  
  Park, Y., Su, Z. and Zhu, H. (2016)  Groupwise envelope 
  models for Imaging Genetic Analysis.
}


