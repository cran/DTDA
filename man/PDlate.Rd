%\encoding{latin1}
\name{PDlate}
\alias{PDlate}
\docType{data}
\title{
Parkinson's Disease Data: late onset
}
\description{
The sample consists of DNA from 100 Caucasian Parkinson's Disease (PD) patients with late onset PD (age 63-87 years).  To remove the selection bias related to survival, the study was limited to patients diagnosed from PD who had their DNA sample taken within eight years after onset. Consequently, the age of onset is doubly truncated by the age at blood sampling and this time minus 8 years. This is a situation of interval sampling, the sampling interval being subject-specific.
}
\usage{data("PDlate")}
\format{
  A data frame with 99 observations on the following 5 variables.
  \describe{
    \item{\code{X}}{a numeric vector, age at onset of PD (in years).}
    \item{\code{U}}{a numeric vector, age at blood sampling minus 8 years. }
    \item{\code{V}}{a numeric vector, age at blood sampling.}
    \item{\code{SNP_A10398G}}{a factor with allels levels \code{A} and \code{G}.}
    \item{\code{SNP_PGC1a}}{a factor with allels levels \code{A}, \code{AG} and \code{G}.}
  }
}
\details{
\cite{Clark et al., 2011} hypothesized that the rs8192678 PGC-1a single nucleotide polymorphism (SNP) and the A10398G mitochondrial SNP may influence risk or age of onset of PD. To test these hypotheses, genomic DNA samples from human blood samples were obtained from the National Institute of Neurological Disorders and Stroke(NINDS) Human Genetics DNA and Cell Line Repositoryat the Coriell Institute for Medical Research (Camden, NewJersey).
}

\source{Mandel M, de Uña-Álvarez J,  Simon DK and Betensky R (2018). Inverse Probability Weighted Cox Regression for Doubly Truncated Data. \emph{Biometrics} \bold{74}, 481-487.
}

\references{
Clark J, Reddy S,  Zheng K,  Betensky RA and Simon DK (2011) Association of PGC-1alphapolymorphisms with age of onset and risk of Parkinson's disease. \emph{BMC Medical Genetics} \bold{12}, 69.

}
\seealso{\code{\link{PDearly}} }

\examples{
data(PDlate)
str(PDlate)
}
\keyword{datasets}
