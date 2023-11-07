\name{is.kadditive-methods}
\docType{methods}
\alias{is.kadditive}
\alias{is.kadditive-methods}
\alias{is.kadditive,Mobius.set.func,numeric-method}
\alias{is.kadditive,card.set.func,numeric-method}
\alias{is.kadditive,set.func,numeric-method}

\title{Test method}

\description{Tests whether a set function is \code{k}-additive, i.e., if its
  \enc{Möbius}{Mobius} function vanishes for subsets of more than k elements. The set
  function can be given either under the form of an object of class
  \code{set.func}, \code{card.set.func} or
  \code{Mobius.set.func}. }

\section{Methods}{
\describe{

\item{object = "Mobius.set.func", k = "numeric", epsilon = "numeric",
  epsilon = "numeric" }{Returns an object of class \code{logical}.}

\item{object = "card.set.func", k = "numeric", epsilon = "numeric" }{Returns an object of class \code{logical}.}

\item{object = "set.func", k = "numeric", epsilon = "numeric" }{Returns an object of class \code{logical}.}
}}

\details{
  In order to test whether a coefficient is equal to zero, its
  absolute value is compared with \code{epsilon} whose default
  value is \code{1e-9}.
}

\references{
  M. Grabisch (1997), \emph{k-order additive discrete fuzzy measures and their
    representation}, Fuzzy Sets and Systems 92(2), pages 167--189.

  M. Grabisch (2000), \emph{The interaction and \enc{Möbius}{Mobius}
  representations of fuzzy measures on finites spaces, k-additive
  measures: a survey}, in: Fuzzy Measures and Integrals: Theory and
  Applications, M. Grabisch, T. Murofushi, and M. Sugeno Eds, Physica
  Verlag, pages 70-93.
  
  
}

\seealso{
  	\code{\link{Mobius.set.func-class}}, 
	\cr \code{\link{card.set.func-class}}, 
	\cr \code{\link{set.func-class}}, 
	\cr \code{\link{Mobius-methods}}, \cr
	\cr \code{\link{k.truncate.Mobius-methods}}. 
}

\examples{
## a set function
mu <- set.func(c(0,1,1,1,2,2,2,3))
mu
is.kadditive(mu,2)
is.kadditive(mu,1)

## the Mobius representation of a set function, 2-additive by construction 
a <- Mobius.set.func(c(0,1,2,1,3,1,2,1,2,3,1),4,2)
is.kadditive(a,2)
is.kadditive(a,1)
}

\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}