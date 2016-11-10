\name{genv}
\alias{genv}

\title{
Fitting the groupwise envelope model
}

\description{
It fits the groupwise envelope model for multivariate linear
regression.  The groupwise envelope model allows for both 
distinct regression coefficients and distinct error structures 
for different groups.
}

\usage{
genv(X, Y, Z, u)
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

\details{
It fits the groupwise envelope model to the responses and
predictors, using the maximum likelihood estimation.  When the
dimension of the envelope is between \eqn{1} and \eqn{r-1}, we implement 
the method in Park et al. (2016).  When the dimension is \eqn{r}, 
the envelope model degenerates to the standard multivariate 
linear regression.  When the dimension is \eqn{0}, it means that 
\code{X} and \code{Y} are uncorrelated and the fitting is 
different. \cr

A standard model in multivariate linear regression model is 
given by
\deqn{\mathbf{Y} = \bm{\mu} + \bm{\beta} \mathbf{X} + \bm{\epsilon},}
where \eqn{\mathbf{Y} \in \R^r} is response vector, 
\eqn{\mathbf{X}} denotes a \eqn{p \times 1} vector of covariates 
and \eqn{\bm{\epsilon}} follows a distribution with mean \eqn{0}
and positive definite covariance matrix 
\eqn{\bm{\Sigma} \in \R^{r \times r}}.  However, this
model is not sufficient to explain the association between 
responses and covatiates differ across different groups.  
Suppose that we observe data in \eqn{L} different groups.  For
each \eqn{l=1,\ldots,L}, the \eqn{l-}th group has \eqn{n_{(l)}} 
observations and the total sample size is 
\eqn{n = \sum_{l=1}^L n_{(l)}}.  By incorporating such group 
information, we reformulate the model as
\deqn{\mathbf{Y}_{(l)j} = \bm{\mu}_{(l)}+\bm{\beta}_{(l)}\mathbf{X}_{(l)j}
+\bm{\epsilon}_{(l)j}
}
for \eqn{l=1, \ldots, L} and \eqn{j=q, \ldots, n_{l}}, where
\eqn{\mathbf{Y}_{(l)j} \in \R^r} is the \eqn{j}th 
observed response vector in the \eqn{l-}th population, 
\eqn{\bm{\mu}_{(l)}\in \R^r} is the mean of the 
\eqn{l-}th population, \eqn{\mathbf{X}_{(l)j} \in \R^p} 
is the \eqn{j}th observed covariate vector in the \eqn{l-}th 
population, \eqn{\bm{\beta}_{(l)}\in \R^{r \times p}} 
contains the regression coefficients for the \eqn{l-}th 
population and \eqn{\bm{\epsilon}_{(l)j}} follows some 
distribution with mean \eqn{\mathbf{0}} and covariance matrix 
\eqn{\bm{\Sigma}_{(l)}}. \cr

Let \eqn{u} denote the dimension of the envelope subspace, 
\eqn{\bm{\Gamma} \in \R^{r \times u}} be an orthogonal 
basis of the subspace and \eqn{\bm{\Gamma}_0 \in \R^{r \times (r-u)}} 
be an orthogonal complemment of \eqn{\bm{\Gamma}}.  The coordinate 
form of the groupwise envelope model is given by
\deqn{
\mathbf{Y}_{(l)j} = \bm{\mu}_{(l)}+\bm{\Gamma}\bm{\eta}_{(l)}\mathbf{X}_{(l)j}
+\bm{\epsilon}_{(l)j}
}
and
\deqn{
\bm{\Sigma}_{(l)} = \bm{\Gamma}\bm{\Omega}_{(l)}\bm{\Gamma}^{\top} +
\bm{\Gamma}_0 \bm{\Omega}_0 \bm{\Gamma}_0^{\top}
}
for each \eqn{l = 1, \ldots, L}, where 
\eqn{\bm{\beta}_{(l)} = \bm{\Gamma} \bm{\eta}_{(l)}}, 
\eqn{\bm{\eta}_{(l)} \in \R^{u \times p}} carries the 
coordinate of \eqn{\bm{\beta}_{(l)}} with respect to 
\eqn{\bm{\Gamma}} and \eqn{\bm{\Omega}_{(l)} \in \R^{u \times u}} 
and \eqn{\bm{\Omega}_0 \in \R^{(r-u) \times (r-u)}} 
are symmetric matrices that carry the corrdinates of 
\eqn{\bm{\Sigma}_{(l)}} with respect to \eqn{\bm{\Gamma}} and 
\eqn{\bm{\Gamma}_0}, respectively.  The groupwsie envelope 
model degenerates to the envelope model is Cook et al. (2010) 
if \eqn{L=1}. 
}

\value{
  \item{mu}{
  The estimator of group mean.  A \eqn{r} by \eqn{L} matrix  
  whose ith column contains the mean for the group.
  }
  \item{beta}{
  A list of \eqn{r} by \eqn{p} matrices for the estimator of regression 
  coefficients.  \code{beta}[[i]] indicates the estimate of regression
  coefficient for the ith group.
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
  (with discussion). \emph{Statistica Sinica} \eqn{\mathbf{20}}, 
  927-1010.

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
