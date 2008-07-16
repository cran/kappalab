\name{veto-methods}
\docType{methods}
\alias{veto}
\alias{veto-methods}
\alias{veto,Mobius.capacity-method}
\alias{veto,capacity-method}
\alias{veto,card.capacity-method}

\title{Veto indices}

\description{Computes the veto indices of a Choquet integral from the
  underlying \bold{normalized} capacity. The capacity can be given either under the form
  of an object of class \code{capacity}, \code{card.capacity} or \code{Mobius.capacity}.}

\section{Methods}{
\describe{

\item{object = "Mobius.capacity" }{The veto indices are computed from the \enc{Möbius}{Mobius} transform of a capacity.}

\item{object = "capacity" }{The veto indices are computed directly from a capacity.}

\item{object = "card.capacity" }{The veto indices are computed from a cardinal capacity.}
}}

\references{
 J.-L. Marichal (2000), \emph{Behavioral analysis of aggregation in
    multicriteria decision aid}, in: Preferences and Decisions under
  Incomplete Knowledge, J. Fodor and B. De Baets and P. Perny Eds,
  Physica-Verlag, pages 153-178.

  J.-L. Marichal (2004), \emph{Tolerant or intolerant character of
    interacting criteria in aggregation by the Choquet integral},
  European Journal of Operational Research 155:3, pages 771-791.
}

\seealso{
  	\code{\link{capacity-class}}, \cr \code{\link{Mobius.capacity-class}},
  	\cr \code{\link{card.capacity-class}}, \cr \code{\link{Mobius-methods}}. 
}

\examples{
## a capacity
mu <- capacity(c(0:13,13,13)/13)

## its Mobius transform
a <- Mobius(mu)

## their veto indices
veto(mu)
veto(a)

## the same with a card.capacity object
mu <- lower.capacity(4)
veto(mu)
}

\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}
