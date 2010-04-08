\encoding{latin1}
\name{DTDA-package}
\alias{DTDA-package}
\alias{DTDA}
\docType{package}
\title{
Doubly truncated data analysis
}
\description{
This package implements different algorithms for analyzing randomly truncated data, one-sided and two-sided (i.e. doubly) truncated data. 
It incorporates the iterative methods introduced by \cite{Efron and Petrosian (1999)} and \cite{Shen (2008)}.
Estimation of the lifetime distribution function and truncation times distributions is possible, together with the corresponding pointwise confidence limits based on bootstrap methods.
Plots of cumulative distributions and survival functions are provided.
}
\details{
\tabular{ll}{
Package: \tab DTDA\cr
Type: \tab Package\cr
Version: \tab 1.2-2\cr
Date: \tab 2010-04-08\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
Missing data are allowed. Registries with missing data are simply removed. This package incorporates the functions \code{efron.petrosian}, \code{lynden} and \code{ shen}, which call the iterative methods introduced by \cite{Efron and Petrosian (1999)} and \cite{Shen (2008)}.

For a complete list of functions, use \code{library(help="DTDA").} }
\author{
Carla Moreira, Jacobo de Uña-Álvarez and Rosa Crujeiras

Maintainer: Carla Moreira \email{carlamgmm@gmail.com}
}
\references{
Efron, B. and Petrosian, V.(1999) Nonparametric methods for doubly truncated data. \emph{Journal of the American Statistical Association} \bold{94}, 824-834.

Lynden-Bell, D. (1971) A method of allowing for known observational selection in small samples applied to 3CR quasars.  \emph{Monograph National Royal Astronomical Society} \bold{155}, 95-118.

Shen, P-S. (2008) Nonparametric analysis of doubly truncated data. \emph{Annals of the Institute of Statistical Mathematics} DOI 10.1007/s10463-008-0192-2.
}
\keyword{ iteration} 
\keyword{ nonparametric}

\section{Acknowledgements}{ Work supported by the research grants MTM 2008-02139 and MTM 2008-0310 of the Spanish Ministry of Science and Innovation, and Xunta de
 Galicia grants PGIDIT07PXIB300191PR and PGIDIT06PXIB207009RR.}