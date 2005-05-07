\name{Sugeno.integral-methods}
\docType{methods}
\alias{Sugeno.integral}
\alias{Sugeno.integral-methods}
\alias{Sugeno.integral,Mobius.game,numeric-method}
\alias{Sugeno.integral,card.game,numeric-method}
\alias{Sugeno.integral,game,numeric-method}

\title{Sugeno integral}

\description{Computes the Sugeno integral of a \bold{non negative} function with
  respect to a game. Moreover, if the game is a capacity, the range of the function must be contained into
  the range of the capacity. The game can be given either under the form of an
  object of class \code{game}, \code{card.game} or \code{Mobius.game}.}

\section{Methods}{
\describe{

  \item{object = "Mobius.game", f = "numeric" }{The Sugeno integral of
    \code{f} is computed from the \enc{Möbius}{Mobius} transform of a game.}

   \item{object = "game", f = "numeric" }{The Sugeno integral of \code{f}
    is computed from a game.}
  
  \item{object = "card.game", f = "numeric" }{The Sugeno integral of
    \code{f} is computed from a cardinal game.}

}}

\references{
  M. Sugeno (1974), \emph{Theory of fuzzy integrals and its applications}, Tokyo
  Institute of Technology, Tokyo, Japan.
  
  J-L. Marichal (2000), \emph{On Sugeno integral as an aggregation function}, Fuzzy Sets
  and Systems 114, pages 347-365.
  
  J-L. Marichal (2001), \emph{An axiomatic approach of the discrete Sugeno integral as a
  tool to aggregate interacting criteria in a qualitative framework}, IEEE
  Transactions on Fuzzy Systems 9:1, pages 164-172.
  
  T. Murofushi and M. Sugeno (2000), \emph{Fuzzy measures and fuzzy integrals},
  in: M. Grabisch, T. Murofushi, and M. Sugeno Eds, Fuzzy Measures and
  Integrals: Theory and Applications, Physica-Verlag, pages 3-41.
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
f <- c(0.1,0.9,0.3,0.8)

## the Sugeno integral of f w.r.t mu
Sugeno.integral(mu,f)
Sugeno.integral(a,f) 

## a similar example with a cardinal capacity
mu <- uniform.capacity(4)
Sugeno.integral(mu,f)
}

\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}
