\name{Sipos.integral-methods}
\docType{methods}
\alias{Sipos.integral}
\alias{Sipos.integral-methods}
\alias{Sipos.integral,Mobius.game,numeric-method}
\alias{Sipos.integral,card.game,numeric-method}
\alias{Sipos.integral,game,numeric-method}

\title{Sipos integral}

\description{Computes the Sipos integral (also called \emph{symmetric Choquet
  integral}) of a real-valued function with respect to a game. The game can be
  given either under the form of an object of class \code{game},
  \code{card.game} or \code{Mobius.game}.}

\section{Methods}{
  \describe{
    
\item{object = "game", f = "numeric" }{The Sipos or symmetric Choquet integral of \code{f}
  is computed from a game.}

\item{object = "Mobius.game", f = "numeric" }{The Sipos or symmetric Choquet integral of
\code{f} is computed from the \enc{Möbius}{Mobius} transform of a game.}

\item{object = "card.game", f = "numeric" }{The Sipos or symmetric Choquet integral of
\code{f} is computed from a cardinal game.}

}}

\references{
  M. Grabisch and Ch. Labreuche (2002), The symmetric and asymmetric Choquet
  integrals on finite spaces for decision making, Statistical Papers 43, pages
  37-52.
}

\seealso{
  \code{\link{game-class}},
  \cr \code{\link{Mobius.game-class}},
  \cr \code{\link{card.game-class}}. 
}

\examples{
## a normalized capacity
mu <- capacity(c(0:13/13,1,1))

## and its Möbius transform
a <- Mobius(mu)

## a discrete function f
f <- c(0.1,-0.9,-0.3,0.8)

## the Sugeno integral of f w.r.t mu
Sipos.integral(mu,f)
Sipos.integral(a,f)

## a similar example with a cardinal capacity
mu <- uniform.capacity(4)
Sipos.integral(mu,f)
}

\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}
