%\encoding{latin1}
\name{DTDA-package}
\alias{DTDA-package}
\alias{DTDA}
\docType{package}
\title{
Doubly Truncated Data Analysis
}
\description{
Implementation of different algorithms for analyzing randomly truncated data, one-sided and two-sided (i.e. doubly) truncated data.
%It incorporates the iterative methods introduced by \cite{Efron and Petrosian (1999)} and \cite{Shen (2010)}.
The package allows for the estimation of the distribution function of the doubly truncated (target) variable. The package provides estimators for the distribution of the truncation variables too. Pointwise confidence limits based on bootstrap methods are included.
Automatic plots of cumulative distributions and survival functions are provided.
The package also implements the kernel density estimator for doubly truncated data with different bandwidth selectors. The hazard rate function with least square cross-validation selector type is also included. Real datasets are provided within the package. Besides the right-truncated AIDS data, eight doubly truncated datasets are available: Childhood Cancer Data, AIDS Blood Transfusion Data, Equipment-S Rounded Failure Time Data, Quasars Data, Parkinson's Disease Data (early and late onset groups), and Acute Coronary Syndrome Data (complete and reduced).
}
\details{
\tabular{ll}{
Package: \tab DTDA\cr
Type: \tab Package\cr
Version: \tab 3.0\cr
Date: \tab 2021-04-08\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
Missing data are allowed. Registries with missing data are simply removed. This package incorporates the functions \code{efron.petrosian}, \code{lynden}, \code{shen} to implement the iterative methods to compute the NPMLE for doubly truncated data; \code{densityDT}, \code{hazardDT} to calculate the kernel estimators of the density and hazard functions, respectively. The function \code{rsim.DT}, allows to simulate doubly truncated data in two different settings.
For a complete list of functions, use \code{library(help="DTDA").} }
\author{
Carla Moreira, Jacobo de Uña-Álvarez and Rosa Crujeiras

Maintainer: Carla Moreira \email{carlamgmm@gmail.com}
}
%\references{


%}

\keyword{iteration}
\keyword{nonparametric}
\keyword{Kernel estimation}
\keyword{double truncation}

\section{Acknowledgements}{
\itemize{
\item{Carla Moreira was partially financed by Portuguese Funds through FCT (Fundação para a Ciência e a Tecnologia) within the Projects UIDB/00013/2020 and UIDP/00013/2020.
}
\item{Jacobo de Uña-Álvarez was supported by the Grant MTM2017-89422-P (MINECO/AEI/FEDER, UE).}
\item{Rosa Crujeiras was supported by Grant MTM2016-76969-P.
}
\item{The authors also thank Hugo S. Oliveira, from University of Porto, for his helpful in improving the programming to reduce the computational burden.}
}}
