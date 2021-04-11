\name{Quasars}
\alias{Quasars}
\docType{data}
\title{ Quasars Data
}
\description{
The original dataset studied by Efron and Petrosian (1999) comprised independlently collected quadruplets of the redshift and the apparent magnitude of a quasar object. Due to experiemtnal constraints, the distribution of each luminosity in a log-scale is truncated to a known interval.
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
\details{
Quadruplets in the original data set studied by \cite{Efron and Petrosian (1999)} are of the form \eqn{(z_i;m_i; a_i; b_i), i = 1, \ldots n}, where \eqn{z_i} is the redshift of the ith quasar and
\eqn{m_i} is the apparent magnitude. Due to experimental constraints, the distribution of each luminosity in the log-scale \eqn{(y_i = t(z_i, m_i))} is truncated to a known interval \eqn{[a_i; b_i]}, where \eqn{t}
represents a transformation which depends on the cosmological model assumed (see \cite{Efron
and Petrosian (1999)} for details). Quasars with apparent magnitude above \eqn{b_i} were too dim to yield dependent redshifts, and hence they were excluded from the study. The lower limit \eqn{a_i}
was used to avoid confusion with non quasar stellar objects.

}
\source{
Vahé Petrosian and Bradley Efron.
}
\references{
Boyle BJ,  Fong R, Shanks, T and Peterson, BA (1990) A catalogue of faint, UV-excess objects. \emph{Monograph National Royal Astronomical Society} \bold{243}, 1-56.

Efron B and Petrosian V (1999) Nonparametric methods for doubly truncated data. \emph{Journal of the American Statistical Association} \bold{94}, 824-834.

}
\examples{
data(Quasars)
 str(Quasars)
}
\keyword{datasets}
