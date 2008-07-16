\name{is.normalized-methods}
\docType{methods}
\alias{is.normalized}
\alias{is.normalized-methods}
\alias{is.normalized,Mobius.capacity-method}
\alias{is.normalized,capacity-method}
\alias{is.normalized,card.capacity-method}

\title{Test method}

\description{Tests whether a capacity is normalized, i.e., if its value on the
  universal set is 1. The capacity can be
  given either under the form of an object of class \code{capacity},
  \code{card.capacity} or \code{Mobius.capacity}.}

\section{Methods}{
\describe{

\item{object = "Mobius.capacity" }{Returns a \code{logical}.}

\item{object = "capacity" }{Returns a \code{logical}.}

\item{object = "card.capacity" }{Returns a \code{logical}.}
}}

\seealso{
  \code{\link{capacity-class}},
  \cr \code{\link{Mobius.capacity-class}},
  \cr \code{\link{card.capacity-class}}.
}

\examples{
## a capacity
mu <- capacity(0:15)
is.normalized(mu)
normalize(mu)

## its Mobius transform
a <- Mobius(mu)
is.normalized(a)
normalize(a)

## a cardinal capacity
mu <- uniform.capacity(7)
is.normalized(mu)
}

\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}
