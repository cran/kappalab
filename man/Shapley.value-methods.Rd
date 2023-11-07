\name{Shapley.value-methods}
\docType{methods}
\alias{Shapley.value}
\alias{Shapley.value-methods}
\alias{Shapley.value,Mobius.set.func-method}
\alias{Shapley.value,card.set.func-method}
\alias{Shapley.value,set.func-method}

\title{The Shapley value}

\description{Computes the Shapley value (\code{n} indices) of a set
  function. The set function can be given either under the form of an
  object of class \code{set.func}, \code{card.set.func} or \code{Mobius.set.func}.}

\section{Methods}{
\describe{

\item{object = "Mobius.set.func" }{The Shapley value is computed from the \enc{Möbius}{Mobius} transform of a set function.}

\item{object = "card.set.func" }{The Shapley value is computed from a cardinal set function.}

\item{object = "set.func" }{The Shapley value is computed from a general set function.}
}}

\references{
  L.S. Shapley (1953), \emph{A value for \code{n}-person games},
  Contributions to the theory of games, vol. 2, Annals of Mathematics
  Studies, no. 28, pages 307-317, Princeton University Press, Princeton,
  N. J. .
}

\seealso{
  \code{\link{Mobius.set.func-class}},
  \cr \code{\link{card.set.func-class}},
  \cr \code{\link{set.func-class}},
  \cr \code{\link{Mobius-methods}}. 
}

\examples{
## a set function
mu <- set.func(c(0:13/13,1,1))

## the Shapley value
Shapley.value(mu)

## the efficiency property should be satisfied
sum(Shapley.value(mu))

## a similar example using a Mobius.set.func object
a <- Mobius(mu)
Shapley.value(a)

## a similar example using a card.set.func object
mu <- upper.capacity(6)
Shapley.value(mu)
## the efficiency property should be satisfied
Shapley.value(mu)*6
}

\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}