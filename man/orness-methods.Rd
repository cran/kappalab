\name{orness-methods}
\docType{methods}
\alias{orness}
\alias{orness-methods}
\alias{orness,Mobius.capacity-method}
\alias{orness,capacity-method}
\alias{orness,card.capacity-method}
\title{Orness degree}

\description{Computes the orness degree of a Choquet integral from the
  underlying \bold{normalized} capacity. The capacity can be given either under the form
  of an object of class \code{capacity}, \code{card.capacity} or
  \code{Mobius.capacity}.}

\section{Methods}{
\describe{

\item{object = "Mobius.capacity" }{The orness degree is computed from the \enc{Möbius}{Mobius} transform of a capacity.}

\item{object = "capacity" }{The orness degree is computed directly from a capacity.}

\item{object = "card.capacity" }{The orness degree is computed from a cardinal capacity.}
}}

\references{
  J.-L. Marichal (2000), \emph{Behavioral analysis of aggregation in
    multicriteria decision aid}, in Preferences and Decisions under
  Incomplete Knowledge, J. Fodor and B. De Baets and P. Perny Eds,
  Physica-Verlag, pages 153-178, 2000.
  
}

\seealso{
  \code{\link{capacity-class}},
  \cr \code{\link{Mobius.capacity-class}},
  \cr \code{\link{card.capacity-class}},
  \cr \code{\link{Mobius-methods}}. 
}

\examples{
## the upper capacity
mu <- capacity(c(0,rep(1,15)))

## the Choquet integral w.r.t mu behaves like the maximum operator
f <- c(0.1,0.1,0,0.9)
Choquet.integral(mu,f)

## its orness is 1
orness(mu)

## the same example with a Mobius.capacity object
a <- Mobius(mu)
Choquet.integral(a,f)
orness(a)

## the same example with a card.capacity object
mu <- upper.capacity(4)
Choquet.integral(mu,f)
orness(mu)
}

\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}
