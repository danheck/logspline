\name{logspline.fit}
\alias{logspline.fit}
\title{
Logspline Density Estimation
}
\section{DESCRIPTION}{
logspline.fit: density estimation using splines to
approximate the log-density.
}
\usage{
logspline.fit(uncensored, right, left, interval, 
              lbound= -Inf, ubound=Inf,
              nknots, knots, penalty, delete=T)
}
\arguments{
\item{uncensored}{
vector of uncensored observations from the distribution whose density is
to be estimated. If there are no uncensored observations, this argument can
be omitted. However, either \code{uncensored} or \code{interval} must be specified.
}
\item{right}{
vector of right censored observations from the distribution
whose density is to be estimated. If there are no right censored
observations, this argument can be omitted.
}
\item{left}{
vector of left censored observations from the distribution
whose density is to be estimated. If there are no left censored
observations, this argument can be omitted.
}
\item{interval}{
two column matrix of lower and upper bounds of observations
that are interval censored from the distribution whose density is
to be estimated. If there are no interval censored observations, this
argument can be omitted.
}
\item{lbound,ubound}{
lower/upper bound for the support of the density. For example, if there
is a priori knowledge that the density equals zero to the left of 0,
and has a discontinuity at 0,
the user could specify \code{lbound=0}. However, if the density is 
essentially zero near 0, one does not need to specify \code{lbound}.
}
\item{nknots}{
forces the method to start with nknots knots (\code{delete=T}) or to fit a
density with \code{nknots} knots (\code{delete=F}). The method has an automatic rule
for selecting \code{nknots} if this parameter is not specified.
}
\item{knots}{
ordered vector of values (that should cover the complete range of the
observations), which forces the method to start with these knots (\code{delete=T})
or to fit a density with these knots (\code{delete=F}). Overrules \code{nknots}.
If \code{knots} is not specified, a default knot-placement rule is employed.
}
\item{penalty}{
the parameter to be used in the AIC criterion. The method chooses
the number of knots that minimizes `
-2*loglikelihood+penalty*(number of knots-1)'.
The default is to use \code{penalty=log(samplesize)} as in BIC. The effect of
this parameter is summarized in \code{logspline.summary}.
}
\item{delete}{
should stepwise knot deletion be employed?
}}
\value{
The output is organized to serve as input for \code{logspline.plot},
\code{logspline.summary}, \code{dlogspline}, \code{plogspline}, \code{qlogspline} and
\code{rlogspline}.
The function returns a list with the following members:

\item{call}{
the command that was executed.
}
\item{knots}{
vector of the locations of the knots in the logspline model.
}
\item{coef}{
coefficients of the spline. The first coefficient is the constant term,
the second is the linear term and the k-th (k>2) is the coefficient
of (x-t(k-2))^3+ (where ^3+ means the positive part of the third power,
and t(k-2) means knot k-2). If a coefficient is zero the corresponding
knot was deleted from the model.
}
\item{bound}{
first element: 0 - \code{lbound} was -infinity, 1 it was something else; second
element: \code{lbound}, if specified; third element: 0 - \code{ubound} was infinity,
1 it was something else; fourth element: \code{ubound}, if specified.
}
\item{logl}{
the k-th element is the log-likelihood of the fit with k+2 knots.
}
\item{penalty}{
the penalty that was used.
}
\item{sample}{
the sample size that was used.
}
\item{delete}{
was stepwise knot deletion employed?
}}
\references{
Kooperberg, C. and Stone, C. J. (1991).  A Study of Logspline Density
Estimation.
\emph{Computational Statistics and Data Analysis}
\bold{12}, 327-347.


Kooperberg, C. and Stone, C. J. (1992). Logspline Density Estimation for
Censored Data. 
\emph{Journal of Computational and Graphical Statistics, }
\bold{1}, 301-328.
}
\seealso{
\code{\link{logspline.plot}}, \code{\link{logspline.summary}},
\code{\link{dlogspline}}, \code{\link{plogspline}}, \code{\link{qlogspline}}, \code{\link{rlogspline}}. 
}
\examples{
estimate the density of a positive random variable when some of the
observations are right censored, and generate 50 random numbers from
the estimated density:


fit <- logspline.fit(y, y2, lbound=0)
rlogspline(50, fit)
}
% Converted by Sd2Rd version 0.2-a5.
