\name{AIDS}
\alias{AIDS}
\docType{data}
\title{ AIDS Blood Transfusion Data, right truncated
}
\description{
The data
include information on the infection and induction times for 258 adults who
were infected with HIV virus and developed AIDS by June 30, 1996. The data consist on the time in years, measured from April 1, 1978, when adults were infected by the virus from a contaminated blood transfusion, and the waiting time to development of AIDS, measured
from the date of infection. The induction times are right-truncated.
}
\usage{data(AIDS)}
\format{
  A data frame with 258 observations on the following 3 variables.
  \describe{
    \item{\code{INFTime}}{a numeric vector, the infection time (years).}
    \item{\code{INDTime}}{a numeric vector, the induction time (years).}
    \item{\code{V}}{a numeric vector, the time from HIV infection to the end of the study (years).}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
Klein and Moeschberger (1997) Survival Analysis Techniques for Censored and truncated data.
Springer.
}
\references{
Lagakos SW and Barraj LM and de Gruttola V (1988) Nonparametric Analysis of Truncated Survival
Data, with Applications to AIDS. \emph{Biometrika} \bold{75}, 515--523.
}
\examples{
data(AIDS)
str(AIDS)
}
\keyword{datasets}
