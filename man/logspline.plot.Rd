\name{logspline.plot}
\alias{logspline.plot}
\title{
Plots a Logspline Density, Distribution Function, Hazard Function or Survival Function.
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
\examples{
fit <- logspline.fit(y)      
logspline.plot(fit)
}
\seealso{
\code{\link{logspline.fit}}, \code{\link{logspline.summary}}, \code{\link{dlogspline}}, \code{\link{plogspline}}, \code{\link{qlogspline}},
\code{\link{rlogspline}}.
}
% Converted by Sd2Rd version 0.2-a5.
