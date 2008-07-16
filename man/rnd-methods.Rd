\name{rnd-methods}
\docType{methods}
\alias{rnd}
\alias{rnd-methods}
\alias{rnd,superclass.set.func-method}

\title{Rounding of set function coefficients}

\description{Rounds the coefficients of a set function to the specified
  number of decimal places (default 0).
}

\section{Methods}{
  \describe{
    
\item{x = "superclass.set.func", digits = "numeric" }{Returns an object of the same class as
  \code{x} whose \code{data} attribute is rounded.}
}}

\examples{
## a capacity
mu <- capacity(0:15/15)
mu
rnd(mu,2)

## a similar example with a Mobius.capacity object
a <- Mobius(mu)
a
rnd(a,1)

## a similar example with a card.capacity object
mu <- uniform.capacity(6)
mu
rnd(mu)
}
\keyword{methods}

