\name{capacity-class}
\docType{class}
\alias{capacity-class}
%\alias{as.card.capacity,capacity-method}
%\alias{conjugate,capacity-method}
%\alias{entropy,capacity-method}
%\alias{favor,capacity-method}
%\alias{is.normalized,capacity-method}
%\alias{Mobius,capacity-method}
%\alias{orness,capacity-method}
%\alias{variance,capacity-method}
%\alias{veto,capacity-method}

\title{Class "capacity"}

\description{Class representing a capacity, i.e. a monotone set function
vanishing at the empty set (also called \emph{fuzzy measure}, \emph{non-additive measure},
\emph{monotone measure}).}

\section{Objects from the Class}{Objects can be mainly created by calls
  to the functions \code{capacity} and \code{entropy.capa.ident}.}

\section{Slots}{
  \describe{
    
     \item{\code{n}:}{Object of class \code{numeric} of length 1 equal to the
      number of elements of the set on which the capacity is defined.}
    
    \item{\code{subsets}:}{Object of class \code{numeric} of length
      \code{2^n} containing the power set of the underlying set in
      "natural" order. The subsets are coded as integers.}
    
    \item{\code{data}:}{Object of class \code{numeric} of length
      \code{2^n} containing the coefficients of the capacity in binary
      order. We necessarily have \code{data[1] = 0} plus some additional monotonicity
      constraints.}
    
  }
}

\section{Extends}{
Class \code{game}, directly.
Class \code{superclass.capacity}, directly.
Class \code{set.func}, by class \code{game}.
Class \code{superclass.set.func}, by class \code{game}.
}

\section{Methods}{
  \describe{
    \item{as.card.capacity}{\code{signature(object = "capacity")}}
    \item{conjugate}{\code{signature(object = "capacity")}}
    \item{entropy}{\code{signature(object = "capacity")}}
    \item{favor}{\code{signature(object = "capacity")}}
    \item{is.normalized}{\code{signature(object = "capacity")}}
    \item{normalize}{\code{signature(object = "capacity")}}
    \item{Mobius}{\code{signature(object = "capacity")}}
    \item{orness}{\code{signature(object = "capacity")}}
    \item{variance}{\code{signature(object = "capacity")}}
    \item{veto}{\code{signature(object = "capacity")}}
  }
}

\seealso{
  \code{\link{capacity}},
  \cr \code{\link{entropy.capa.ident}},
  \cr \code{\link{as.card.capacity-methods}},
  \cr \code{\link{conjugate-methods}}, 
  \cr \code{\link{entropy-methods}},
  \cr \code{\link{favor-methods}}, 
  \cr \code{\link{is.normalized-methods}},
  \cr \code{\link{normalize-methods}},
  \cr \code{\link{Mobius-methods}}, 
  \cr \code{\link{orness-methods}},
  \cr \code{\link{variance-methods}}, 
  \cr \code{\link{veto-methods}}. 
}

\examples{
## a capacity
mu <- capacity(c(0:13,13,13)/13)

## the attributes of the object
mu@n
mu@data
mu@subsets

## a test
is.normalized(mu)
normalize(mu)

## a conversion that should not work
## as.card.capacity(mu)

## some transformations
conjugate(mu)
Mobius(mu)
## let us check ...
zeta(Mobius(mu))

## some summary indices
orness(mu)
veto(mu)
favor(mu)
variance(mu)
entropy(mu)
## the same
summary(mu)
}

\keyword{classes}
