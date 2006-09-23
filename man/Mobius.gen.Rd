\name{Mobius.set.func}
\alias{Mobius.set.func}
\alias{Mobius.game}
\alias{Mobius.capacity}
\alias{additive.capacity}

\title{Create objects representing the Möbius transform of a set function.}

\description{These functions create objects of class
  \code{Mobius.set.func}, \code{Mobius.game}, or \cr
  \code{Mobius.capacity} from an object of class \code{numeric}.}

\usage{
Mobius.set.func(object, n, k)
Mobius.game(object, n, k)
Mobius.capacity(object, n, k)
additive.capacity(v)
}

\arguments{
  \item{n}{An object of class \code{numeric} of length 1 equal to the
    cardinal of the set on which the \enc{Möbius}{Mobius} representation is defined.}
  
  \item{k}{An object of class \code{numeric} of length 1 equal to the
    order of truncation of the \enc{Möbius}{Mobius} representation. Coefficients whose
    cardinal is superior to \code{k} are considered to be zero.}
  
  \item{object}{Object of class \code{numeric} of length
      \code{choose(n,0) + ... + choose(n,k)} representing the
      coefficients of the truncated \enc{Möbius}{Mobius} transform in
      "natural" order.}
    
  \item{v}{An object of class \code{numeric} containing the coefficients
    of the capacity on singletons.}  
}

\value{
 Return objects of class \code{Mobius.set.func}, \code{Mobius.game}, or \code{Mobius.capacity}.
}


\seealso{
	\code{\link{Mobius.capacity-class}},  
	\cr \code{\link{Mobius.game-class}},
	\cr \code{\link{Mobius.set.func-class}},
	\cr \code{\link{k.truncate.Mobius-methods}}.
}

\examples{
Mobius.set.func(c(0,1,1,2,1,2,1,2,2,3,2),4,2)
Mobius.game(c(0,1,1,2,1,2,1,2,2,3,2),4,2)
Mobius.capacity(c(0,1,1,2,1,2,1,2,2,3,2),4,2)
additive.capacity(c(1,2,3,4))
}
\keyword{math}
