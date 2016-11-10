\name{fiberpaper}
\alias{fiberpaper}
\docType{data}

\title{
Data for pulp and paper property
}

\description{
This is a data for pulp and paper property.  The original data
contains measurements of properties of pulp fibers and the paper.
It has 62 observations on 8 variables. 
}

\usage{data(fiberpaper)}

\format{
  A data frame with 62 by 8 matrix.
  \describe{
    \item{\code{data}}{
    A 62 by 8 matrix.  The first column denotes breaking length, 
    the second column denotes elastic modulus, the third column
    denotes stress at failure, the fourth column denotes burst 
    strength, the fifth column denotes arithmetic fiber length, 
    the sixth column denotes long fiber fraction, the seventh
    column denotes fine fiber fraction and the eighth column
    denotes zero span tensile.
    }
  }
}

\references{
  Johnson, R. A. and Wichern, D. W. (2007)  Applied Multivariate
  Statistical Analysis, 6th edition.
}

\examples{
data(fiberpaper)
Y <- as.matrix(data[1:4])
XX <- as.numeric(data[,7] > mean(data[,7]))
X <- as.vector(XX)
Z <- as.matrix(data[,5:6])
}

\keyword{datasets}
