\name{card.set.func}
\alias{card.set.func}
\alias{card.game}
\alias{card.capacity}
\alias{lower.capacity}
\alias{upper.capacity}
\alias{uniform.capacity}

\title{Create objects representing cardinal set functions.}

\description{These functions create objects of class
  \code{card.set.func}, \code{card.game}, or \code{card.capacity}
  from an object of class \code{numeric}.}

\usage{
card.set.func(object)
card.game(object)
card.capacity(object)
lower.capacity(n)
upper.capacity(n)
uniform.capacity(n)
}

\arguments{
  \item{object}{An object of class \code{numeric} containing the
    coefficients of the cardinal set function. The coefficient at
    position \code{k} corresponds to the value of the cardinal set
    function for subsets of size \code{k}.}
  
  \item{n}{An object of class \code{numeric} of length 1 containing
    the cardinal of the set on which the cardinal set function is
    defined.}
}

\value{
Return objects of class \code{card.set.func}, \code{card.game}, or \code{card.capacity}.
}

\seealso{
	\code{\link{card.capacity-class}}, 
	\cr \code{\link{card.game-class}},
	\cr \code{\link{card.set.func-class}}.
}

\examples{
card.set.func(4:-2)
card.game(c(0,-2:2))
card.capacity(0:5)
lower.capacity(3)
upper.capacity(4)
uniform.capacity(5)
}
\keyword{math}