\name{variance-methods}
\docType{methods}
\alias{variance}
\alias{variance-methods}
\alias{variance,Mobius.capacity-method}
\alias{variance,capacity-method}
\alias{variance,card.capacity-method}

\title{Normalized variance of a capacity}

\description{Computes the normalized variance of a capacity. The
  capacity can be given either under the form of an object of class \code{capacity}, a \code{card.capacity} or \code{Mobius.capacity}.}

\section{Methods}{
\describe{

\item{object = "Mobius.capacity" }{The normalized variance is computed from the \enc{Möbius}{Mobius} transform of a capacity.}

\item{object = "capacity" }{The normalized variance is computed directly from a capacity.}

\item{object = "card.capacity" }{The normalized variance is computed from a cardinal capacity.}

}}


\references{
I. Kojadinovic (2005), \emph{Minimum variance capacity identification}, submitted.
}

\seealso{
  	\code{\link{capacity-class}}, \cr \code{\link{Mobius.capacity-class}},
  	\cr \code{\link{card.capacity-class}}, \cr \code{\link{Mobius-methods}}. 
}

\examples{
## a capacity
mu <- capacity(c(0,0,0,0:12)/12)

## its Möbius transform
a <- Mobius(mu)

## their normalized variance
variance(mu)
variance(a)

## similar examples with card.capacity objects
mu <- lower.capacity(4)
variance(mu)
mu <- uniform.capacity(4)
variance(mu)
}

\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}
