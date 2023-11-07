\name{zeta-methods}
\docType{methods}
\alias{zeta}
\alias{zeta-methods}
\alias{zeta,Mobius.capacity-method}
\alias{zeta,Mobius.card.set.func-method}
\alias{zeta,Mobius.game-method}
\alias{zeta,Mobius.set.func-method}

\title{The zeta transform}


\description{Computes the zeta transform of a set function given under
  the form of its \enc{Möbius}{Mobius} transform. The zeta transform is the inverse of
  the \enc{Möbius}{Mobius} transform.
}

\section{Methods}{
\describe{

\item{object = "Mobius.capacity" }{Returns an object of class \code{capacity}.}

\item{object = "Mobius.card.set.func" }{Returns an object of class \code{card.set.func}.}

\item{object = "Mobius.game" }{Returns an object of class \code{game}.}

\item{object = "Mobius.set.func" }{Returns an object of class \code{set.func}.}
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
	\code{\link{capacity-class}}, 
	\cr \code{\link{card.set.func-class}}, 
	\cr \code{\link{game-class}},
	\cr \code{\link{set.func-class}},
  	\cr \code{\link{Mobius.capacity-class}}, 
	\cr \code{\link{Mobius.card.set.func-class}}, 
	\cr \code{\link{Mobius.game-class}},
	\cr \code{\link{Mobius.set.func-class}}, 
	\cr \code{\link{Mobius-methods}}. 
}

\examples{
## the Mobius transform of a capacity
a <- Mobius.capacity(c(rep(0,15),1),4,4)
a

## its zeta transform
zeta(a)

## let us check
Mobius(zeta(a))

## a similar example with a Mobius.card.set.func object
mu <- card.set.func(-3:4)
a <- Mobius(mu)
zeta(a)
}

\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}
