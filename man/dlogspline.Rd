\name{dlogspline}
\alias{dlogspline}
\alias{plogspline}
\alias{qlogspline}
\alias{rlogspline}
\title{
Logspline Density Estimation
}
\description{
  Density, cumulative probability, quantiles and random samples from
  a logspline density.
}
\usage{
dlogspline(x, fit)
plogspline(q, fit)
qlogspline(p, fit)
rlogspline(n, fit)
}
\arguments{
  \item{x, q}{
    vector of quantiles. Missing values (NAs) are allowed.
  }
  \item{p}{
    vector of probabilities. Missing values (NAs) are allowed.
  }
  \item{n}{
    sample size. If \code{length(n)} is larger than 1, then
    \code{length(n)} random values are returned.
  }
  \item{fit}{
    a list like the output from \code{logspline.fit}.
  }
}
\value{
  densities (\code{dlogspline}), probabilities (\code{plogspline}),
  quantiles (\code{qlogspline}),
  or a random sample (\code{rlogspline}) from a logspline density.
}
\details{
  Elements of \code{q} or \code{p} that are missing will cause the
  corresponding elements of the result to be missing.
}
\seealso{
  \code{\link{logspline.fit}}, \code{\link{logspline.plot}},
  \code{\link{logspline.summary}}
}
\keyword{smooth}
