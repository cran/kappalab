\name{Mobius.card.set.func}
\alias{Mobius.card.set.func}


\title{Creates an object representing the Möbius transform of a cardinal set function.}

\description{Creates objects of class \code{Mobius.card.set.func}: from
  an object of class \code{numeric}.}

\usage{
Mobius.card.set.func(object)
}

\arguments{
   \item{object}{An object of class \code{numeric} containing the
    coefficients of the cardinal \enc{Möbius}{Mobius} representation. The coefficient
    at position \code{k} corresponds to the value of the cardinal
    \enc{Möbius}{Mobius} representation for subsets of size \code{k}.}
}

\value{
 Returns an object of class \code{Mobius.card.set.func}.
}


\seealso{
	\code{\link{Mobius.card.set.func-class}}.
}

\examples{
Mobius.card.set.func(4:-2)
}
\keyword{math}
