\name{Mobius.card.set.func-class}
\docType{class}
\alias{Mobius.card.set.func-class}
%\alias{as.Mobius.set.func,Mobius.card.set.func-method}
%\alias{to.data.frame,Mobius.card.set.func-method}
%\alias{zeta,Mobius.card.set.func-method}

\title{Class "Mobius.card.set.func"}

\description{Class representing the \enc{Möbius}{Mobius} transform of a cardinal
  set function.}


\section{Objects from the Class}{
  Objects can be created by calls to the function
  \code{Mobius.card.set.func}.}

\section{Slots}{
  \describe{

    \item{\code{n}:}{Object of class \code{numeric} of length 1 equal
      to the number of elements of the set on which the \enc{Möbius}{Mobius}
      representation of a cardinal set function is defined.}
    
    \item{\code{data}:}{Object of class \code{numeric} of length
      \code{n+1} containing the coefficients of the \enc{Möbius}{Mobius}
      representation of a cardinal set function.}
  }
}

\section{Extends}{
  Class \code{superclass.set.func}, directly.}

\section{Methods}{
  \describe{

    \item{as.Mobius.set.func}{\code{signature(object =
	"Mobius.card.set.func")}}

    \item{as.set.func}{\code{signature(object =
	"Mobius.card.set.func")}}

    \item{as.card.set.func}{\code{signature(object = "Mobius.card.set.func")}}
    
    \item{to.data.frame}{\code{signature(object = "Mobius.card.set.func")}}
    
    \item{zeta}{\code{signature(object = "Mobius.card.set.func")}}}
  
  
}

\seealso{
  \code{\link{Mobius.set.func-class}},
  \cr \code{\link{Mobius.card.set.func}},
  \cr \code{\link{as.Mobius.set.func-methods}}, 
  \cr \code{\link{as.card.set.func-methods}},
  \cr \code{\link{as.set.func-methods}},
  \cr \code{\link{zeta-methods}}, 
  \cr \code{\link{to.data.frame-methods}}.
}
\examples{
## the Mobius representation of a cardinal set function
a <- Mobius.card.set.func(-3:2)

# the attributes of the object
a@n
a@data

## some transformations
as.set.func(a)
zeta(a)
## let us check ...
Mobius(zeta(a))
}
\keyword{classes}
