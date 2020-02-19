\name{Quasars}
\alias{Quasars}
\docType{data}
\title{ Quasars
}
\description{
The original data set studied by Efron and Petrosian (1999), comprised independently collected
quadruplets \eqn{(z_i;m_i; a_i; b_i), i = 1, \ldots n}, where \eqn{z_i} is the redshift of the ith quasar and
\eqn{m_i} is the apparent magnitude. Due to experimental constraints, the distribution of each
luminosity in the log-scale \eqn{(y_i = t(z_i, m_i))} is truncated to a known interval \eqn{[a_i; b_i]}, where \eqn{t}
represents a transformation which depends on the cosmological model assumed (see Efron
and Petrosian (1999) for details). Quasars with apparent magnitude above \eqn{b_i} were too dim to
yield dependent redshifts, and hence they were excluded from the study. The lower limit \eqn{a_i}
was used to avoid confusion with non quasar stellar objects.
}
\usage{data(Quasars)}
\format{
  A data frame with 210 observations on the following 3 variables.
  \describe{
    \item{\code{y (adj lum)}}{a numeric vector, the log lominosity values.}
    \item{\code{u (lower)}}{a numeric vector, lower truncation limits.}
    \item{\code{v (upper)}}{a numeric vector, upper truncation limits.}
  }
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
Vahé Petrosian and Bradley Efron.
}
\references{
Boyle, B.J. and Fong, R. and Shanks, T. and Peterson, B.A. (1990)  \emph{Monograph National Royal Astronomical Society} \bold{243}, 1.

Efron, B. and Petrosian, V.(1999) Nonparametric methods for doubly truncated data. \emph{Journal of the American Statistical Association} \bold{94}, 824-834.

}
\examples{
data(Quasars)
 str(Quasars)
}
\keyword{datasets}
