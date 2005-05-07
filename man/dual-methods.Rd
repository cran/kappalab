\name{dual-methods}
\docType{methods}
\alias{dual}
\alias{dual-methods}
\alias{dual,capacity-method}
\alias{dual,card.capacity-method}
\alias{dual,card.set.func-method}
\alias{dual,set.func-method}

\title{The dual (or conjugate) transform}

\description{Computes the dual (also called the \emph{conjugate}) of a set function. The
  dual of the dual of a game gives the original game back.}

\section{Methods}{
\describe{

\item{object = "capacity" }{Returns an object of class \code{capacity}.}

\item{object = "card.capacity" }{Returns an object of class \code{card.capacity}.}

\item{object = "card.set.func" }{Returns an object of class \code{card.game}.}

\item{object = "set.func" }{Returns an object of class \code{game}. }
}}

\references{
  T. Murofushi and M. Sugeno (2000), \emph{Fuzzy measures and fuzzy integrals},
  in: M. Grabisch, T. Murofushi, and M. Sugeno Eds, Fuzzy Measures and
  Integrals: Theory and Applications, Physica-Verlag, pages 3-41.
}

\seealso{
  	\code{\link{capacity-class}}, 
	\cr \code{\link{card.capacity-class}}, 
	\cr \code{\link{set.func-class}},
	\cr \code{\link{card.set.func-class}}. 
}

\examples{
## a game
mu <- game(c(0,-7:7))
mu

## its dual 
dual(mu)
## and mu again 
dual(dual(mu))

## a similar example with the upper capacity
mu <- capacity(c(0,rep(1,15)))
mu
dual(mu)
dual(dual(mu))

## a similar example with an object of class card.capacity
mu <- upper.capacity(6)
mu
dual(mu)
dual(dual(mu))

## the dual of a set function is a game
mu <- set.func(-7:8)
mu
dual(mu)
mu <- card.set.func(-2:5)
dual(mu)
}

\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}
