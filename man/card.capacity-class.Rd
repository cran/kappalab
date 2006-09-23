\name{card.capacity-class}
\docType{class}
\alias{card.capacity-class}
%\alias{as.capacity,card.capacity-method}
%\alias{conjugate,card.capacity-method}
%\alias{entropy,card.capacity-method}
%\alias{favor,card.capacity-method}
%\alias{is.normalized,card.capacity-method}
%\alias{orness,card.capacity-method}
%\alias{variance,card.capacity-method}
%\alias{veto,card.capacity-method}

\title{Class "card.capacity"}

\description{Class representing a cardinal capacity, i.e. a capacity whose
values 
depend only on the cardinality of subsets (also called \emph{symmetric capacity}).}

\section{Objects from the Class}{Objects can be created by calls to the functions \code{card.capacity}, \code{lower.capacity}, \code{upper.capacity}, \code{uniform.capacity}.}

\section{Slots}{
  \describe{

    \item{\code{n}:}{Object of class \code{numeric} of length 1 equal
      to the number of elements of the set on which the cardinal game is defined.}
    
    \item{\code{data}:}{Object of class \code{numeric} of length
      \code{n+1} containing the coefficients of the cardinal capacity. We necessarily have
      \code{data[1]=0} and \code{data[i+1]-data[i]>0}.} 
  }
}

\section{Extends}{
Class \code{card.game}, directly.
Class \code{superclass.capacity}, directly.
Class \code{card.set.func}, by class \code{card.game}.
Class \code{superclass.set.func}, by class \code{card.game}.
}

\section{Methods}{
  \describe{
    \item{as.capacity}{\code{signature(object = "card.capacity")} }
    \item{conjugate}{\code{signature(object = "card.capacity")} }
    \item{entropy}{\code{signature(object = "card.capacity")} }
    \item{favor}{\code{signature(object = "card.capacity")} }
    \item{is.normalized}{\code{signature(object = "card.capacity")} } 
    \item{normalize}{\code{signature(object = "card.capacity")} }
    \item{orness}{\code{signature(object = "card.capacity")} }
    \item{variance}{\code{signature(object = "card.capacity")} }
    \item{veto}{\code{signature(object = "card.capacity")} }
  }
}

\seealso{
	\code{\link{capacity-class}}, \cr \code{\link{card.capacity}},
  	\cr \code{\link{as.capacity-methods}}, \cr \code{\link{conjugate-methods}}, 
	\cr \code{\link{entropy-methods}}, \cr \code{\link{favor-methods}}, 
	\cr \code{\link{is.normalized-methods}}, \cr \code{\link{orness-methods}}, 
	\cr \code{\link{variance-methods}}, \cr \code{\link{veto-methods}}.
}

\examples{
## a capacity
mu <- card.capacity(0:6/6)
## the same
mu <- uniform.capacity(6)

# the attributes of the object
mu@n
mu@data

## a test
is.normalized(mu)
normalize(mu)

## a transformation
conjugate(mu)

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
