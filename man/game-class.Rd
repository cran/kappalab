\name{game-class}
\docType{class}
\alias{game-class}
%\alias{as.card.game,game-method}
%\alias{Choquet.integral,game,numeric-method}
%\alias{Mobius,game-method}
%\alias{Sipos.integral,game,numeric-method}
%\alias{Sugeno.integral,game,numeric-method}

\title{Class "game"}

\description{Class representing a game, i.e. a set function vanishing at
the empty set (also called \emph{non monotonic fuzzy measure}).}

\section{Objects from the Class}{Objects can be created by calls to the function \code{game}.}

\section{Slots}{
  \describe{
    \item{\code{n}:}{Object of class \code{numeric} of length 1 equal to the
      number of elements of the set on which the game is defined.}
    
    \item{\code{subsets}:}{Object of class \code{numeric} of length
      \code{2^n} containing the power set of the underlying set in
      "natural" order. The subsets are coded as integers.}
    
    \item{\code{data}:}{Object of class \code{numeric} of length
      \code{2^n} containing the coefficients of the game in binary
      order.  We necessarily have \code{data[1] = 0}.}
    
  }
}

\section{Extends}{
Class \code{set.func}, directly.
Class \code{superclass.set.func}, by class \code{set.func}.
}

\section{Methods}{
  \describe{
    \item{as.card.game}{\code{signature(object = "game")}}
    \item{Choquet.integral}{\code{signature(object = "game", f = "numeric")}}
    \item{Mobius}{\code{signature(object = "game")}}
    \item{Sipos.integral}{\code{signature(object = "game", f = "numeric")}}
    \item{Sugeno.integral}{\code{signature(object = "game", f =
	"numeric")}}
    \item{pdf.Choquet.unif}{\code{signature(object = "game", f = "numeric")}}
    \item{cdf.Choquet.unif}{\code{signature(object = "game", f = "numeric")}}
    \item{expect.Choquet.unif}{\code{signature(object = "game")}}
    \item{sd.Choquet.unif}{\code{signature(object = "game")}}
    \item{expect.Choquet.norm}{\code{signature(object = "game")}}
    \item{sd.Choquet.norm}{\code{signature(object = "game")}}
  }
}

\seealso{	
  \code{\link{game}},
  \cr \code{\link{as.card.game-methods}},
  \cr \code{\link{Choquet.integral-methods}}, 
  \cr \code{\link{Mobius-methods}},
  \cr \code{\link{Sipos.integral-methods}}, 
  \cr \code{\link{Sugeno.integral-methods}},
  \cr \code{\link{pdf.Choquet.unif-methods}}, 
  \cr \code{\link{cdf.Choquet.unif-methods}},
  \cr \code{\link{expect.Choquet.unif-methods}}, 
  \cr \code{\link{sd.Choquet.unif-methods}},
  \cr \code{\link{expect.Choquet.norm-methods}}, 
  \cr \code{\link{sd.Choquet.norm-methods}}. 
}

\examples{
## a game (which is a capacity)
mu <- game(c(0,rep(1,15)))

## the attributes of the object
mu@n
mu@data
mu@subsets

## a conversion
as.card.game(mu)

## a transformation
Mobius(mu)
## let us check ...
zeta(Mobius(mu))

## integral calculations 
f <- c(0.2,0.3,0.1,0.7)
Choquet.integral(mu,f)
Sugeno.integral(mu,f)
f <- c(0.2,-0.3,0.1,-0.7)
Sipos.integral(mu,f)
}

\keyword{classes}
