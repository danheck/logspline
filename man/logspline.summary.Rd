\name{logspline.summary}
\alias{logspline.summary}
\title{
Summarizes Various Logspline Density Estimates
}
\usage{
logspline.summary(fit)
}
\arguments{
\item{fit}{
a list like the output from logspline.fit
}}
\value{
This function produces only printed output. The main body
is a table with five columns: the first column is a possible number
of knots for the fitted model;
the second column is the log-likelihood for the fit;
the third column is -2*loglikelihood + penalty*(number of knots-1),
which is the AIC criterion - logspline.fit selected the model with
the smallest value of AIC;
the fourth and fifth columns give the
endpoints of the interval of values of penalty that would yield the
model with the indicated number of knots. (NAs imply that the model is
not optimal for any choice of penalty.) At the bottom of the table the
number of knots corresponding to the selected model is reported, as is
the value of penalty that was used.
}
\seealso{
\code{\link{logspline.fit}}, \code{\link{logspline.plot}}, \code{\link{dlogspline}}, \code{\link{plogspline}}, \code{\link{qlogspline}}, \code{\link{rlogspline}}.
}
\examples{
fit <- logspline.fit(y) 
logspline.summary(fit)
}
% Converted by Sd2Rd version 0.2-a5.
