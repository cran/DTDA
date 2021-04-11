%\encoding{latin1}
\name{ChildCancer}
\alias{ChildCancer}
\docType{data}
\title{
Childhood Cancer Data
}
\description{
This dataset corresponds to all children diagnosed from cancer between January
1, 1999  and December 31, 2003  in the region of North Portugal. The database includes information of 406 children with complete records on the age at diagnosis. Because of the interval sampling, the age at diagnosis is doubly truncated by the time  from birth to the end of the study, and time from birth to the beggining of the study (time in days). The age at diagnosis is doubly truncated due to the interval sampling.
}
\usage{data("ChildCancer")}
\format{
  A data frame with 406 observations on the following 8 variables.
  \describe{
    \item{\code{X}}{a numeric vector, age at diagnosis(time in days).}
    \item{\code{U}}{a numeric vector, time from birth to the beggining of the study (time in days).}
\item{\code{V}}{a numeric vector, time from birth to the end of the study (time in days).}
\item{\code{ICCGroup}}{a numeric vector,  cancer types
 identified according to the International Classification of Childhood Cancer
(ICCC). 1=Leukaemias, 2=Limphoma, 3=Nervous System Tumour, 4=Neuroblastoma, 5=Retinoblastoma, 6=Renal, 7=Hepatic, 8=Bone, 9=Soft Tissues, 10=Germ Cell, 11=Melanoma and other epitelial tumours, 12=Other Tumours.}
\item{\code{Status}}{a numeric vector,  the status indicator at the end of the study:  0=alive, 1=dead.}
\item{\code{SurvTime}}{a numeric vector,  the survival time (time from birth to death or end of the study.}
\item{\code{Residence}}{a numeric vector,  districts of residence. 1=Braga, 2=Bragança,3=Porto,4=Viana do Castelo, 5=Vila Real.}
\item{\code{Sex}}{a numeric vector, sex of the participants (1 = female, 2 = male).}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
The childhood cancer data were gathered from the IPO (Registo Oncológico do Norte) service, kindly provided by Doctor Maria José Bento.
}
\references{
Moreira C and  de Uña-Álvarez J (2010) Bootstrapping the NPMLE for doubly truncated data. \emph{Journal of Nonparametric Statistics} \bold{22}, 567-583.
}
\examples{
data(ChildCancer)
str(ChildCancer)
}
\keyword{datasets}
