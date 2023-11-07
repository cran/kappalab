\name{interaction.indices-methods}
\docType{methods}
\alias{interaction.indices}
\alias{interaction.indices-methods}
\alias{interaction.indices,Mobius.set.func-method}
\alias{interaction.indices,card.set.func-method}
\alias{interaction.indices,set.func-method}

\title{The Shapley interaction indices}

\description{
Computes the Shapley interaction index for pairs of elements with
respect to a set function. The set function can be given either under
the form of an object of class \code{set.func}, \code{card.set.func}
or \code{Mobius.set.func}
.}

\section{Methods}{
\describe{

\item{object = "Mobius.set.func" }{The Shapley interaction indices are computed from the \enc{Möbius}{Mobius} transform of a set function.}

\item{object = "card.set.func" }{The Shapley interaction indices are computed from a cardinal set function.}

\item{object = "set.func" }{The Shapley interaction indices are computed from a general set function.}
}}

\references{
  G. Owen (1971/72), \emph{Multilinear extensions of games}, Management Sci. 18, pages 64--79.
     
  T. Murofushi and S. Soneda (1993), \emph{Techniques for reading fuzzy measures
    (III): interaction index}, 9th Fuzzy System Symposium, pages 693--696, Saporo, Japan.

  M. Grabisch, J-L. Marichal and M. Roubens (2000), \emph{Equivalent
    representations of set functions}, Mathematics of Operations
  Research 25(2), pages 157--178.

}

\seealso{
  \code{\link{Mobius.set.func-class}},
  \cr \code{\link{card.set.func-class}},
  \cr \code{\link{set.func-class}},
  \cr \code{\link{Mobius-methods}}. 
}

\examples{
## a set function
mu <- set.func(c(-7:6,6,6))

## the associated interaction indices
interaction.indices(mu)


## a similar example using a Mobius.set.func object
a <- Mobius(mu)
interaction.indices(a)

## a similar example using a card.set.func object
mu <- upper.capacity(6)
interaction.indices(mu)
}

\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}
