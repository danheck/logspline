\name{logspline.plot}
\alias{logspline.plot}
\title{
Plots a Logspline Density, Distribution Function, Hazard Function or Survival Function.
}
\description{
  Plots a logspline density, distribution function, hazard function or
  survival function.
}
\synopsis{
logspline.plot(fit, n=100, what="d", xlim, xlab, ylab, type, \dots)
}
\usage{
logspline.plot(fit, n=100, what="d", \dots)
}
\arguments{
\item{fit}{
a list like the output from \code{logspline.fit}.
}
\item{n}{
the number of equally spaced points at which to plot the density.
}
\item{what}{
what should be plotted: \code{d} (density), \code{p} (distribution function), \code{s} (survival
function) or \code{h} (hazard function).
}
\item{...}{
all regular plotting options as desired.
}}
\value{
This function produces a plot of a logspline fit at n equally
spaced points roughly covering the support of the density. (Use
\code{xlim=c(from,to)} to change the range of these points.)
}
\seealso{
  \code{\link{logspline.fit}}, \code{\link{logspline.summary}},
  \code{\link{dlogspline}}, \code{\link{plogspline}},
  \code{\link{qlogspline}}, \code{\link{rlogspline}}.
}
\keyword{smooth}
\keyword{hplot}
