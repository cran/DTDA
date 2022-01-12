%\encoding{latin1}
\name{hazardDT}
\alias{hazardDT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of the kernel hazard function under random double truncation
}
\description{
This function provides the nonparametric kernel hazard estimation for a variable which is observed under random double truncation, which is defined as a convolution of a kernel function with the NPMLE of the cumulative hazard. Least square cross validation bandwidth selection procedure is implemented too.
}
\usage{
hazardDT(X, U, V, bw = "LSCV", from, to, n, wg = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Numeric vector with the values of the target variable.
}
\item{U}{
Numeric vector with the values of the left truncation variable.

}
\item{V}{
Numeric vector with the values of the rigth truncation variable.

}
  \item{bw}{
The smoothing bandwidth to be used, but can also be a character string giving a rule to choose the bandwidth. This must be  \code{"LSCV"}.
}

\item{from}{
the left  point of the grid at which the density is to be estimated. The default is min(X)+1e-04.
}

\item{to}{
the rigth  point of the grid at which the density is to be estimated. The default is max(X)-1e-04.
}
\item{n}{
number of evaluation points on a equally spaced grid.
}
\item{wg}{
Numeric vector of random weights to correct for double truncation. Default weights correspond to the Efron-Petrosian NPMLE.
}

}

\details{
The nonparametric kernel density estimation for a variable which is observed under random double truncation is computed  as proposed in \cite{Moreira et al.(2021)}. As usual in kernel smoothing, the estimator is  obtained as a convolution between a kernel function and an appropriate  estimator of the cumulative hazard. Gaussian kernel is used.
The automatic bandwidth selection procedures for the kernel hazard estimator is the least square cross validation, presented in \cite{Moreira et al. (2021)}.
}

\value{
A list containing the following values:
\item{x }{the n coordinates of the points where the hazard is estimated. }
 \item{y}{the estimated hazard values.}
 \item{bw}{the bandwidth used.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}

\references{
Moreira C, de Uña-Álvarez J, Santos AC and Barros H (2021) Smoothing Methods to estimate the hazard rate under double truncation. \url{https://arxiv.org/abs/2103.14153}.
}
\author{ Carla Moreira, Jacobo de Uña-Álvarez and Rosa Crujeiras }
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{densityDT}} }
\examples{

set.seed(4321)

n<-100
X <- runif(n, 0, 1)
U <- runif(n,-1/3, 1)
V <- U + 1/3
for (i in 1:n){

	while (U[i] > X[i] | V[i] < X[i]){
	X[i] <- runif(1, 0, 1)
  U[i] <- runif(1, -1/3, 1)
	V[i] <- U[i] + 1/3
	}

}


vxhazard1<-hazardDT(X,U,V,bw=0.3,0,1,500)
vxhazard2<-hazardDT(X,U,V,bw="LSCV",0,1,500)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Kernel hazard estimation }
 \keyword{Bandwidth selection }
  \keyword{Double truncation }
