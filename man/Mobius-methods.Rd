\name{Mobius-methods}
\docType{methods}
\alias{Mobius}
\alias{Mobius-methods}
\alias{Mobius,capacity-method}
\alias{Mobius,card.set.func-method}
\alias{Mobius,game-method}
\alias{Mobius,set.func-method}

\title{The Möbius transform}

\description{Computes the \enc{Möbius}{Mobius} transform of a set function. The \enc{Möbius}{Mobius} transform is the inverse of the zeta transform.}

\section{Methods}{
\describe{

\item{object = "capacity" }{Returns an object of class \code{Mobius.capacity}.}

\item{object = "card.set.func" }{Returns an object of class \code{Mobius.card.set.func}.}

\item{object = "game" }{Returns an object of class \code{Mobius.game}.}

\item{object = "set.func" }{Returns an object of class \code{Mobius.set.func}.}
}}

\references{
  G-C. Rota (1964), \emph{On the foundations of combinatorial theory. I. Theory
    of \enc{Möbius}{Mobius} functions}, Z. Wahrscheinlichkeitstheorie und
  Verw. Gebiete 2, pages 340-368.
  
  
  A. Chateauneuf and J-Y. Jaffray (1989), \emph{Some characterizations of lower 
  probabilities and other monotone capacities through the use of \enc{Möbius}{Mobius} inversion},
  Mathematical Social Sciences 17, pages 263-283.
  
  M. Grabisch, J-L. Marichal and M. Roubens (2000), \emph{Equivalent
    representations of set functions}, Mathematics of Operations
  Research 25:2, pages 157-178.
  
}

\seealso{
	\code{\link{Mobius.capacity-class}}, 
	\cr \code{\link{Mobius.card.set.func-class}}, 
	\cr \code{\link{Mobius.game-class}},
	\cr \code{\link{Mobius.set.func-class}},
  	\cr \code{\link{capacity-class}}, 
	\cr \code{\link{card.set.func-class}}, 
	\cr \code{\link{game-class}},
	\cr \code{\link{set.func-class}}, 
	\cr \code{\link{zeta-methods}}. 
}

\examples{
## a capacity
mu <- capacity(0:15)
mu

## its Mobius transform
a <- Mobius(mu)
a

## its zeta transform
zeta(a)

## a similar example with a game object
mu <- game(c(0,-2:12))
mu
Mobius(mu)
zeta(Mobius(mu))

## a similar example with a set.func object
mu <- set.func(-7:8)
mu
Mobius(mu)
zeta(Mobius(mu))

## a similar example with a card.set.func object
mu <- card.set.func(-3:4)
mu
Mobius(mu)
zeta(Mobius(mu))
}

\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}