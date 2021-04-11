%\encoding{latin1}
\name{ACSred}
\alias{ACSred}
\docType{data}
\title{
Acute Coronary Syndrome reduced data
}
\description{
The data include information of 917 patients with confirmed diagnosis of type 1 (primary spontaneous) acute coronary syndrome (ACS). Patients were consecutively admitted to the Cardiology Department of two tertiary hospitals in Portugal between August 2013 and December 2014. The age at diagnosis is doubly truncated because of the interval sampling.
}
\usage{data(ACSred)}
\format{
 A data frame with 917 observations on the following 5 variables.
  \describe{
    \item{\code{X}}{a numeric vector, age at diagnosis (in years).}
    \item{\code{U}}{a numeric vector, the elapsed time
(in years) between birth and the beggining of the study (August
2013).}
    \item{\code{V}}{a numeric vector, the elapsed time
(in years) between birth and end of the study (December
2014).}
    \item{\code{Sex}}{a numeric vector, sex of the participants (0 = female, 1 = male).}
    \item{\code{diagnosis}}{a numeric vector, type of diagnosis at discharge	1 - STEMI (ST elevation myocardial infarction ) and 2 - NSTEACS (all others diagnoses).}
  }
}
\details{The age at diagnosis \code{X} is doubly truncated   due to the interval sampling. The length of the sampling interval (\code{V}-\code{U}) is 1.42 years. The NPMLE of the cumulative distribution function of \code{X} for the complete data does not exist or is not unique for this dataset. This dataset is a reduced sample of the original ACS data, guaranteeing the existence and uniqueness of the NPMLE, according to \cite{Xiao and Hudgens (2020)}.
}
%\source{
%%%  ~~ reference to a publication or URL from which the data were obtained %~~
%}
\references{
Araújo C, Laszczynska O, Viana M, Melão F, Henriques A,  Borges A,  Severo M, Maciel MJ, Moreira I, Azevedo A (2018) Sex differences in presenting symptoms of acute coronary syndrome: the EPIHeart cohort study. \emph{BMJ Open } \bold{8}.

Xiao J and  Hudgens MG (2020) On nonparametric maximum likelihood estimation with double truncation. \emph{Biometrika} \bold{106}, 989-996.
}

\seealso{\code{\link{ACS}} }
\examples{
data(ACSred)
str(ACSred)
}
\keyword{datasets}
