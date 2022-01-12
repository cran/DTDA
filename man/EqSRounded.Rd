\name{EqSRounded}
\alias{EqSRounded}
\docType{data}
\title{
Equipment-S Rounded Failure Time Data
}
\description{
Digitized data from Figure 2 in \cite{Ye and Tang 2016}. The dataset contains (rounded) observations of 174 failure times of certain devices, observed under interval sampling. Right-runcation is years between installation and 2011 and left truncation corresponds to right-truncation time minus 15 years. The failure time is doubly truncated due to the interval sampling.
}
\usage{data("EqSRounded")}
\format{
  A data frame with 174 observations on the following 3 variables.
  \describe{
\item{\code{X}}{a numeric vector, time to failure in years.}
\item{\code{U}}{a numeric vector, the number of years between installation and 2011 minus 15 years.}
\item{\code{V}}{a numeric vector, the number of years between installation and 2011.}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
Digitalization of the data plot in the original paper of \cite{Ye and Tang 2016}.
}
\references{
Ye ZS and Tang LC (2016) Augmenting the unreturned for field data with information on returned failures only. \emph{Technometrics} \bold{58}, 513-523.
}
\examples{
data(EqSRounded)
str(EqSRounded)
}
\keyword{datasets}
