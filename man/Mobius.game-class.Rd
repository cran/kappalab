\name{Mobius.game-class}
\docType{class}
\alias{Mobius.game-class}
%\alias{Choquet.integral,Mobius.game,numeric-method}
%\alias{Sipos.integral,Mobius.game,numeric-method}
%\alias{Sugeno.integral,Mobius.game,numeric-method}
%\alias{zeta,Mobius.game-method}

\title{Class "Mobius.game"}

\description{Class representing the \enc{Möbius}{Mobius} transform of a game.}

\section{Objects from the Class}{Objects can be created by calls to the function \code{Mobius.game}.}

\section{Slots}{
  \describe{

    \item{\code{n}:}{Object of class \code{numeric} of length 1 containing the
      number of elements of the set on which the \enc{Möbius}{Mobius} transform is
      defined. }
    
    \item{\code{k}:}{Object of class \code{numeric} of length 1 containg the order
      of truncation of the \enc{Möbius}{Mobius} transform: subsets whose cardinal is
      superior to k are considered to be zero.}
    
   \item{\code{subsets}:}{Object of class \code{numeric} containing
      the "\code{k} power set" of the underlying
      set in "natural" order . The subsets are encoded as integers.}
    
    \item{\code{data}:}{Object of class \code{numeric} of length
      \code{choose(n,0) + ... + choose(n,k)} representing the
      coefficients of a truncated \enc{Möbius}{Mobius} transform of a game in "natural" order.}
    
  }
}

\section{Extends}{
Class \code{Mobius.set.func}, directly.
Class \code{superclass.set.func}, by class \code{Mobius.set.func}.
}

\section{Methods}{
  \describe{

    \item{Choquet.integral}{\code{signature(object = "Mobius.game", f =
	"numeric")} }
    
    \item{Sipos.integral}{\code{signature(object = "Mobius.game", f =
	"numeric")} }
    
    \item{Sugeno.integral}{\code{signature(object = "Mobius.game", f =
	"numeric")} }
    
    \item{zeta}{\code{signature(object = "Mobius.game")} }
    
  }
}

\seealso{
  \code{\link{game-class}},
  \cr \code{\link{Mobius.game}},
  \cr \code{\link{Choquet.integral-methods}},
  \cr \code{\link{Sipos.integral-methods}}, 
  \cr \code{\link{Sugeno.integral-methods}},
  \cr \code{\link{zeta-methods}},
  \cr \code{\link{expect.Choquet.norm-methods}}.
}

\examples{
## a game (which is a capacity)
mu <- game(c(0,rep(1,15)))
## and its Mobius representation
a <- Mobius(mu)

# the attributes of object a
a@n
a@k
a@data
a@subsets

## a transformation
zeta(a)
## let us check ...
Mobius(zeta(a))

## integral calculations 
f <- c(0.2,0.3,0.1,0.7)
Choquet.integral(a,f)
Sugeno.integral(a,f)
f <- c(0.2,-0.3,0.1,-0.7)
Sipos.integral(a,f)
}

\keyword{classes}
