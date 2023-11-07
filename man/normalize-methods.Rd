\name{normalize-methods}
\docType{methods}
\alias{normalize}
\alias{normalize-methods}
\alias{normalize,Mobius.capacity-method}
\alias{normalize,capacity-method}
\alias{normalize,card.capacity-method}

\title{Normalizes a capacity.}
\description{Returns the normalized version of a capacity. The capacity
  can be given either under the form of an object of class
  \code{capacity}, \code{card.capacity} or \code{Mobius.capacity}. In
  the case of objects of class \code{capacity} and \code{card.capacity},
  the normalization is performed by dividing each coefficient of the
  capacity by the value of the capacity on the universe. In the case of
  objects of class \code{Mobius.capacity}, the normalization is
  performed by dividing each coefficient of the \enc{Möbius}{Mobius} transform by the
  sum of all the coefficients.  
}
\section{Methods}{
\describe{

\item{object = "Mobius.capacity" }{Returns an object of class \code{Mobius.capacity}.}

\item{object = "capacity" }{Returns an object of class \code{capacity}.}

\item{object = "card.capacity" }{Returns an object of class \code{card.capacity}.}
}}


\seealso{
  \code{\link{capacity-class}},
  \cr \code{\link{Mobius.capacity-class}},
  \cr \code{\link{card.capacity-class}},
  \cr \code{\link{Mobius-methods}},
  \cr \code{\link{is.normalized-methods}}. 
}

\examples{
## a capacity
mu <- capacity(0:15)

## normalize it
is.normalized(mu)
normalize(mu)

## a similar example with a Mobius.capacity object
a <- Mobius(mu)
is.normalized(a)
normalize(a)
zeta(normalize(a))

## a similar example with a card.capacity object
mu <- card.capacity(0:6)
is.normalized(mu)
normalize(mu)
Mobius(normalize(mu))

}

\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}