%\encoding{latin1}
\name{lynden}
\alias{lynden}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ NPMLE computation with the second Efron-Petrosian algorithm }
\description{
  This function computes the NPMLE for the cumulative distribution function of  \code{X} observed under one-sided (right or left) and two-sided (double) truncation.
It provides simple bootstrap pointwise confidence limits too. This function allows for ties in the samples of \code{X}, \code{U} and \code{V}.

}
\usage{
lynden(X, U = NA, V = NA, error = NA, nmaxit = NA,
		 boot = TRUE, B = NA, alpha = NA,
			display.F = FALSE, display.S = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ Numeric vector with the values of the target variable.}
  \item{U}{Numeric vector with the values of the left truncation variable. If there are no truncation values from the left, put \code{U}=\code{NA}.}
  \item{V}{ Numeric vector with the values of the right truncation variable. If there are no truncation values from the right, put \code{V}=\code{NA}.}
   \item{error}{ Numeric value. Maximum pointwise error when estimating the density associated to X (f) in two consecutive steps. If this is missing, it is $1e{-06}$.}
  \item{nmaxit}{ Numeric value. Maximum number of iterations. If this is missing, it is set to \code{nmaxit} =100 . }
  \item{boot}{ Logical. If TRUE (default), the simple bootstrap method is applied to lifetime distribution estimation.
 Pointwise confidence bands are provided. }
  \item{B}{ Numeric value. Number of bootstrap resamples . The default \code{NA} is equivalent to \code{B} =500 . }
  \item{alpha}{ Numeric value. (1-\code{alpha}) is the nominal coverage for the pointwise confidence intervals. }
  \item{display.F}{ Logical. Default is FALSE. If TRUE, the estimated cumulative distribution function associated to \code{X}, (F) is plotted. }
  \item{display.S}{ Logical. Default is FALSE. If TRUE, the estimated survival function associated to \code{X}, (S) is plotted. }


}
\details{
  The NPMLE for the cumulative distribution function is computed by the second algorithm proposed in \cite{Efron and Petrosian (1999)}. This is an iterative algorithm which converges to the
NMPLE after a number of iterations.
If the second (respectively third) argument is missing, the Lynden-Bell estimator for right-truncated (respectively left-truncated)
 data is obtained.
Note that individuals with NAs in the three first arguments will be automatically excluded.

}
\value{
A list containing the following values:


  \item{time }{The timepoint on the curve. }
  \item{n.event }{The number of events that ocurred at time \code{t}. }
  \item{events}{The total number of events. }
  \item{NJ }{The number of individuals in risk considering the left truncation times. }
  \item{density }{The estimated density values. }
  \item{cumulative.df }{The estimated cumulative distribution values. }
  \item{truncation.probs }{The probability of \code{X} falling into each truncation interval. }
  \item{hazard }{The estimated hazard values. }
  \item{S0 }{\code{error} reached in the algorithm. }
  \item{Survival }{The estimated survival values. }
  \item{n.iterations }{The number of iterations used by this algorithm. }
  \item{B }{Number of bootstrap resamples computed. }
  \item{alpha }{The nominal level used to construct the confidence intervals. }
  \item{upper.df }{The estimated upper limits of the confidence intervals for F. }
  \item{lower.df }{The estimated lower limits of the confidence intervals for F. }
  \item{upper.Sob }{The estimated upper limits of the confidence intervals for S. }
  \item{lower.Sob }{The estimated lower limits of the confidence intervals for S. }
   \item{sd.boot }{The bootstrap standard deviation of F estimator.}
\item{boot.repeat }{The number of resamples done in each bootstrap call to ensure the existence and uniqueness of the bootstrap NPMLE. }
}
\references{ Efron B and Petrosian V (1999) Nonparametric methods for doubly truncated data. \emph{Journal of the American Statistical Association} \bold{94}, 824-834.

Lynden-Bell D (1971) A method of allowing for known observational selection in small samples applied to 3CR quasars. \emph{Monograph National Royal Astronomical Society} \bold{155}, 95-118.}

\author{ Carla Moreira, Jacobo de Uña-Álvarez and Rosa Crujeiras }

\seealso{\code{\link{efron.petrosian}} }
\examples{

#  Generating data which are doubly truncated

set.seed(4321)

n<-25
X<-runif(n,0,1)
U<-runif(n,0,0.25)
V<-runif(n,0.75,1)

for (i in 1:n){
	while (X[i]<U[i]|X[i]>V[i]){
	U[i]<-runif(1,0,0.25)
	X[i]<-runif(1,0,1)
	V[i]<-runif(1,0.75,1)
					   }
		 }
res<-lynden(X=X, U=U, V=V, boot=FALSE, display.F=TRUE, display.S=TRUE)

# Generating data which are right truncated

set.seed(4321)

n<-25
X<-runif(n,0,1)
V<-runif(n,0.75,1)

for (i in 1:n){
	while (X[i]>V[i]){
	X[i]<-runif(1,0,1)
	V[i]<-runif(1,0.75,1)
			     }
		 }
res<-lynden(X=X,U=NA, V=V, boot=FALSE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ iteration }
\keyword{ nonparametric }% __ONLY ONE__ keyword per line
