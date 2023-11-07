\name{Mobius.capacity-class}
\docType{class}
\alias{Mobius.capacity-class}
%\alias{entropy,Mobius.capacity-method}
%\alias{favor,Mobius.capacity-method}
%\alias{is.normalized,Mobius.capacity-method}
%\alias{orness,Mobius.capacity-method}
%\alias{variance,Mobius.capacity-method}
%\alias{veto,Mobius.capacity-method}
%\alias{zeta,Mobius.capacity-method}

\title{Class "Mobius.capacity"}

\description{Class representing the \enc{Möbius}{Mobius} transform of a capacity.}

\section{Objects from the Class}{Objects can be mainly created by calls
  to the functions \code{Mobius.capacity},
  \cr \code{mini.var.capa.ident},
  \code{ls.sorting.capa.ident},
  and \code{least.squares.capa.ident}.}

\section{Slots}{
  \describe{

     \item{\code{n}:}{Object of class \code{numeric} of length 1 containing the
      number of elements of the set on which the \enc{Möbius}{Mobius} transform is
      defined. }
    
    \item{\code{k}:}{Object of class \code{numeric} of length 1 containg the order of
      truncation of the \enc{Möbius}{Mobius} transform: the value of subsets whose cardinal is
      superior to k is put to zero.}
    
    \item{\code{subsets}:}{Object of class \code{numeric} containing
      the "\code{k} power set" of the underlying
      set in "natural" order . The subsets are encoded as integers.}
    
    \item{\code{data}:}{Object of class \code{numeric} of length
      \code{choose(n,0) + ... + choose(n,k)} representing the
      coefficients of a truncated \enc{Möbius}{Mobius} transform of a capacity in "natural" order.}
    
  }
}

\section{Extends}{
Class \code{Mobius.game}, directly.
Class \code{superclass.capacity}, directly. \cr
Class \code{Mobius.set.func}, by class \code{Mobius.game}.
Class \code{superclass.set.func}, by class \code{Mobius.game}.
}

\section{Methods}{
  \describe{
    \item{entropy}{\code{signature(object = "Mobius.capacity")} }
    \item{favor}{\code{signature(object = "Mobius.capacity")} }
    \item{is.normalized}{\code{signature(object = "Mobius.capacity")} }
    \item{normalize}{\code{signature(object = "Mobius.capacity")} }
    \item{orness}{\code{signature(object = "Mobius.capacity")} }
    \item{variance}{\code{signature(object = "Mobius.capacity")} }
    \item{veto}{\code{signature(object = "Mobius.capacity")} }
    \item{zeta}{\code{signature(object = "Mobius.capacity")} }
  }
}

\seealso{
  \code{\link{capacity-class}},
  \cr \code{\link{entropy-methods}},
  \cr \code{\link{favor-methods}}, 
  \cr \code{\link{is.normalized-methods}},
  \cr \code{\link{orness-methods}}, 
  \cr \code{\link{variance-methods}},
  \cr \code{\link{veto-methods}},
  \cr \code{\link{zeta-methods}},
  \cr \code{\link{mini.var.capa.ident}},
  \cr \code{\link{least.squares.capa.ident}},
  \cr \code{\link{ls.sorting.capa.ident}}.

}

\examples{
## a capacity
mu <- capacity(c(0,0,0:13))
## and its Mobius representation
a <- Mobius(mu)
a

# the attributes of object a
a@n
a@k
a@data
a@subsets

## a test
is.normalized(a)
## normalize it
normalize(a)

## a transformation
zeta(a)
## Let us check ...
Mobius(zeta(a))

## some summary indices
orness(a)
veto(a)
favor(a)
variance(a)
entropy(a)
## the same
summary(a)
}

\keyword{classes}
