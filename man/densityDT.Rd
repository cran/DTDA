%\encoding{latin1}
\name{densityDT}
\alias{densityDT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Estimation of a kernel density function under random double truncation
}
\description{
This function provides the nonparametric kernel density estimation of a doubly truncated random variable. A bandwidth value is required.

}
\usage{
densityDT(X, U, V, bw = "DPI2", from, to, n, wg = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
numeric vector with the values of the target variable.
}
  \item{U}{
numeric vector with the values of the left truncation variable.

}
  \item{V}{
numeric vector with the values of the rigth truncation variable.

}
  \item{bw}{
The smoothing bandwidth to be used, but can also be a character string giving a rule to choose the bandwidth. This must be one of \code{"NR"}, \code{"DPI1"}, \code{"DPI2"}, \code{"LSCV"} or \code{"SBoot"} with default \code{"DPI2"}.
}

  \item{from}{
the left  point of the grid at which the density is to be estimated.
}

 \item{to}{
the rigth  point of the grid at which the density is to be estimated.
}
  \item{n}{
number of evaluation points on a equally spaced grid.
}
  \item{wg}{
Numeric vector of non-negative initial solution, with the same length as \code{X}. Default value is the solution obtained with Efron and Petrosian algotithm.
}


}

\details{
The nonparametric kernel density estimation for a variable which is observed under random double truncation is computed  as proposed in \cite{Moreira and de Uña-Álvarez (2012)}. As usual in kernel smoothing, the estimator is  obtained as a convolution between a kernel function and an appropriate  estimator of the cumulative df. Gaussian kernel is used.
The automatic bandwidth selection procedures for the kernel density
estimator are those proposed in \cite{Moreira and Van Keilegom (2013)}. The automatic bandwidth selection alternatives are appropriate modifications, i.e, taking into account the double truncation issue,  of the
normal reference rule, two types of plug-in procedures, the least squares cross-validation and a bootstrap based method proposed in \cite{Cao et al. (1994)} and \cite{Sheater and Jones (1991)} for the complete data.
}

\value{
A list containing the following values:
\item{x }{the n coordinates of the points where the density is estimated. }
 \item{y}{the estimated density values.}
 \item{bw}{the bandwidth used.}
}
\references{



Cao R, Cuevas A and González-Manteiga W (1994). A comparative study of several smoothing methods in density estimation. \emph{Computational Statistics and Data
Analysis} \bold{17}, 153-176.


Moreira C and de Uña-Álvarez J (2012) Kernel density estimation with doubly
truncated data. \emph{Electronic Journal of Statistics} \bold{6}, 501-521.

Moreira C and Van Keilegom I (2013) Bandwidth selection for kernel density estimation with
doubly truncated data. \emph{Computational Statistics and Data Analysis} \bold{61}, 107-123.


Sheather S and Jones M (1991) A reliable data-based bandwidth selection method for kernel density estimation. \emph{Journal of the Royal Statistical Society: Series
B} \bold{53}, 683-690.

Silverman BW (1986) Density Estimation. London: Chapman and Hall.

}
\author{ Carla Moreira, de de Uña-Álvarez and Rosa Crujeiras }

%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{hazardDT}} }
\examples{


n<-50
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


vxDens1<-densityDT(X,U,V,bw="DPI1",0,1,500)
vxDens2<-densityDT(X,U,V,bw="DPI2",0,1,500)
vxDens3<-densityDT(X,U,V, bw=0.5,0,1,500)
vxDens4<-densityDT(X,U,V,bw="LSCV",0,1,500)


data(Quasars)
densityDT(Quasars[,1],Quasars[,2],Quasars[,3],bw="DPI1",-2.5,2.2,500)
densityDT(Quasars[,1],Quasars[,2],Quasars[,3], bw=0.5,-2.5,2.2,500)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
 \keyword{ ~Kernel density estimation }
 \keyword{ ~Bandwidth selection }
  \keyword{ ~Double truncation }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.

