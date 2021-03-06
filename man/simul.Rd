\name{simul}
\alias{simul}
\docType{data}

\title{
Simulation data
}

\description{
This is a simulation data.
}

\usage{data(simul)}

\format{
  A data frame with 100 observations on the following 3 variables.
  \describe{
    \item{\code{X}}{
    A numeric vector for group indicator.  It takes 
    either 0 or 1.  The first group for 0 has 40 observations and 
    the second group for 1 has 60 observations.
    }
    \item{\code{Y}}{
    A numeric response marix with 10 response
    variables.
    }
    \item{\code{Z}}{
    A numeric matrix with 10 covariates.
    }
  }
}

\details{
We consider two groups which have 40\% and 60\% of the 
observations.  We set true dimension of the envelope 
subspace 1.  The orthogonal basis \eqn{\bm{\Gamma}}{%
\Gamma} of the envelope subspace and its complement 
\eqn{\bm{\Gamma}_0}{%
\Gamma_0} are obtained by normalizing an \eqn{r} by 
\eqn{r} matrix of independent normal variates.  The first
group has mean \eqn{3} and the second group has mean 
\eqn{10}.  Generate \eqn{\eta_1}{%
\eta_1} from the chi squared variates with degrees of 
freedom \eqn{10} and \eqn{\eta_{2}}{%
\eta_2} from the chi squared variates with degrees of 
freedom \eqn{5}.  We also generate \eqn{\Omega_1=a^2, 
\Omega_2=b^2} and \eqn{\Omega_0 = A A^{\top}}{%
\Omega_0 = A A^T}, where \eqn{a} and \eqn{b} are both 
standard normal variates and \eqn{A} is a \eqn{(r-u)} 
by \eqn{(r-u)} matrix of independent normal variates 
with mean \eqn{5} and variance \eqn{1}.  The predictors 
are independent normal variates with mean zero and 
variance \eqn{25} for the first group and independent 
normal variates with mean zero and variance \eqn{100} 
for the second group.
}

\references{
  Park, Y., Su, Z. and Zhu, H. (2016)  Groupwise envelope models 
  for Imaging Genetic Analysis. 
}

\examples{
data(simul)
}

\keyword{datasets}
