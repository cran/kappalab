\name{Choquet.integral-methods}
\docType{methods}
\alias{Choquet.integral}
\alias{Choquet.integral-methods}
\alias{Choquet.integral,Mobius.game,numeric-method}
\alias{Choquet.integral,card.game,numeric-method}
\alias{Choquet.integral,game,numeric-method}

\title{Choquet integral}

\description{Computes the Choquet integral of a discrete function with
  respect to a game. The game can be given either under the form of an
  object of class \code{game}, \code{card.game} or
  \code{Mobius.game}. If the integrand is not positive, this function
  computes what is known as the \emph{asymmetric Choquet integral}.}


\section{Methods}{
\describe{

  \item{object = "Mobius.game", f = "numeric" }{The Choquet integral of
    \code{f} is computed from the \enc{Möbius}{Mobius} transform of a game.}

   \item{object = "game", f = "numeric" }{The Choquet integral of
    \code{f} is computed from a game.}

  \item{object = "card.game", f = "numeric" }{The Choquet integral of
    \code{f} is computed from  a cardinal game.}
}}

\references{
  
  G. Choquet (1953), \emph{Theory of capacities}, Annales de l'Institut
  Fourier 5, pages 131-295.

  D. Denneberg (2000), \emph{Non-additive measure and integral, basic concepts and their
  role for applications}, in: M. Grabisch, T. Murofushi, and M. Sugeno Eds, Fuzzy
  Measures and Integrals: Theory and Applications, Physica-Verlag, pages 42-69.
  
  M. Grabisch, T. Murofushi, M. Sugeno Eds (2000), \emph{Fuzzy Measures and
    Integrals: Theory and Applications}, Physica-Verlag.
  
  M. Grabisch and Ch. Labreuche (2002), \emph{The symmetric and asymmetric Choquet
  integrals on finite spaces for decision making}, Statistical Papers 43, pages
  37-52.
  
  M. Grabisch (2000), \emph{A graphical interpretation of the Choquet integral}, IEEE
  Transactions on Fuzzy Systems 8, pages 627-631.
  
  J.-L. Marichal (2000), \emph{An axiomatic approach of the discrete Choquet integral as
  a tool to aggregate interacting criteria}, IEEE Transactions on Fuzzy Systems
  8:6, pages 800-807.
  
  Murofushi and M. Sugeno (1993), \emph{Some quantities represented by the Choquet
  integral}, Fuzzy Sets and Systems 56, pages 229-235.
  
  Murofushi and M. Sugeno (2000), \emph{Fuzzy measures and fuzzy integrals}, in: M.
  Grabisch, T. Murofushi, and M. Sugeno Eds, Fuzzy Measures and Integrals: Theory
  and Applications, Physica-Verlag, pages 3-41.
  
}

\seealso{
  	\code{\link{game-class}}, \cr \code{\link{Mobius.game-class}},
  	\cr \code{\link{card.game-class}}. 
}

\examples{
## a normalized capacity
mu <- capacity(c(0:13/13,1,1))

## and its Möbius transform
a <- Mobius(mu)

## a discrete positive function f
f <- c(0.1,0.9,0.3,0.8)

## the Choquet integral of f w.r.t mu
Choquet.integral(mu,f)
Choquet.integral(a,f)

## a similar example with a cardinal capacity
mu <- uniform.capacity(4)
Choquet.integral(mu,f)
}

\keyword{methods}
