%\encoding{latin1}
\name{rsim.DT}
\alias{rsim.DT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Random generation functions of doubly truncated data
}
\description{
Random generation functions of doubly truncated data with two different patterns of observational bias.
}
\usage{
rsim.DT(n,tau, model=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{n}{ number of observations to generate.}
\item{tau}{ length of the observational window.}
\item{model}{model to be simulated. Number 1 or 2 corresponding to different patterns of observacional bias.}
}
\details{
If \code{model}=1, \eqn{U\sim Unif(-\code{tau},1)} and V= U+ \code{tau}. If \code{model}=2, \eqn{U \sim Unif(0,1)^2\times (\code{tau}+1)-\code{tau}} and V= U+ \code{tau}. In \code{model}=1 there is no observational bias due double truncation while in \code{model}=2 double truncation induces observational bias.
}
\value{
A matrix with \code{n} unit length rows representing the generated values from a doubly truncated data with triplets \eqn{(X, U and V)}, in which \eqn{(U \leq X \leq V)}.
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{ Carla Moreira, Jacobo de Uña-Álvarez and Rosa Crujeiras
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

set.seed(4321)
rsim.DT(500,1/2, model=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
 \keyword{double truncation }
 \keyword{observational bias }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
