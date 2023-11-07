\name{k.truncate.Mobius-methods}
\docType{methods}
\alias{k.truncate.Mobius}
\alias{k.truncate.Mobius-methods}
\alias{k.truncate.Mobius,set.func,numeric-method}
\alias{k.truncate.Mobius,Mobius.set.func,numeric-method}

\title{k-order truncation of the Möbius representation of a set function.}

\description{Truncates the \enc{Möbius}{Mobius} representation of a set function by
  considering that the values of subsets whose cardinal is superior to \code{k}
  are zero. The result is at most \code{k}-additive. }

\section{Methods}{
\describe{

\item{object = "set.func", k = "numeric" }{Returns an object of class \code{Mobius.set.func}. The \enc{Möbius}{Mobius} representation of the set function is first computed and then \code{k}-truncated.}

\item{object = "Mobius.set.func", k = "numeric" }{Returns an object of class \code{Mobius.set.func}.}
}}

\references{
  M. Grabisch (1997), \emph{k-order additive discrete fuzzy measures and their
    representation}, Fuzzy Sets and Systems 92(2), pages 167-189.

  M. Grabisch (2000), \emph{The interaction and \enc{Möbius}{Mobius}
  representations of fuzzy measures on finites spaces, k-additive
  measures: a survey}, in: Fuzzy Measures and Integrals: Theory and
  Applications, M. Grabisch, T. Murofushi, and M. Sugeno Eds, Physica
  Verlag, pages 70-93.
}

\seealso{
  	\code{\link{Mobius.set.func-class}}, 
	\cr \code{\link{set.func-class}}, 
	\cr \code{\link{Mobius-methods}}. 
}

\examples{
## a set function
mu <- set.func(c(0,1,1,1,2,2,2,3))
mu

## 2-truncate it
k.truncate.Mobius(mu,2)
## 2-truncate it
k.truncate.Mobius(Mobius(mu),2)
}

\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}