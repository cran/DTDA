\name{AIDS.DT}
\alias{AIDS.DT}
\docType{data}
\title{
AIDS Blood Transfusion Data
}
\description{
The data include information of transfusions cases of transfusion-related AIDS, corresponding to individuals diagnosed prior to July 1, 1986.  Only 295 patients with consistent data, for which the infection could be attributed to
a single transfusion or short series were included. Since HIV was unknown before 1982, this implies that cases developing AIDS prior to this date were not reported, leading to a doubly truncated data. The incubation time is doubly truncated due to the interval sampling.
}
\usage{data(AIDS.DT)}

\format{
  A data frame with 295 observations on the following 4 variables.
\describe{
\item{\code{X}}{a numeric vector, the induction or incubation time: time elapsed from HIV infection to AIDS (in months).}
\item{\code{U}}{a numeric vector, time from 1982 to HIV infection (in months). }
\item{\code{V}}{a numeric vector,time from HIV infection to July 1, 1986 (in months).}
\item{\code{AGE}}{a numeric vector, age of the individual at diagnosis (in months).}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
Kalbfleisch JD and Lawless JF
}
\references{
Kalbfleisch JD and Lawless JF (1989) Inference based on retrospective ascertainment: An analysis of the data on transfusion-related AIDS. \emph{Journal of the American Statistical Association} \bold{84}, 360--372.
}
\examples{
data(AIDS.DT)
str(AIDS.DT)
}
\keyword{datasets}


