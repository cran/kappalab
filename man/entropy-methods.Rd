\name{entropy-methods}
\docType{methods}
\alias{entropy}
\alias{entropy-methods}
\alias{entropy,Mobius.capacity-method}
\alias{entropy,capacity-method}
\alias{entropy,card.capacity-method}

\title{Normalized entropy of a capacity}

\description{Computes the normalized entropy of a capacity. The capacity
  can be given either under the form of an object of class
  \code{capacity}, \code{card.capacity} or \code{Mobius.capacity}.}

\section{Methods}{
\describe{

  \item{object = "Mobius.capacity" }{The normalized entropy is computed
    from the
    \enc{Möbius}{Mobius} transform of a capacity.}

\item{object = "capacity" }{The normalized entropy is computed directly from a capacity.}

\item{object = "card.capacity" }{The normalized entropy is computed from a cardinal capacity.}

}}

\references{
  J-L. Marichal (2002), \emph{Entropy of discrete Choquet capacities}, European
  Journal of Operational Research, 3:137, 2002, pages 612-624.

  I. Kojadinovic, J-L. Marichal and M. Roubens (2005), \emph{An axiomatic
    approach to the definition of the entropy of a discrete Choquet
    capacity}, Information Sciences 172, pages 131-153.
}

\seealso{
  	\code{\link{capacity-class}}, \cr \code{\link{Mobius.capacity-class}},
  	\cr \code{\link{card.capacity-class}}, \cr \code{\link{Mobius-methods}}. 
}

\examples{
## a capacity
mu <- capacity(c(0,0,0:13))

## its Mobius transform
a <- Mobius(mu)

## their normalized entropy
entropy(mu)
entropy(a)

## similar examples with card.capacity objects
mu <- lower.capacity(4)
entropy(mu)
mu <- uniform.capacity(4)
entropy(mu)
}

\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}
