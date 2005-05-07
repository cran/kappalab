\name{card.game-class}
\docType{class}
\alias{card.game-class}
%\alias{as.game,card.game-method}
%\alias{Choquet.integral,card.game,numeric-method}
%\alias{Sipos.integral,card.game,numeric-method}
%\alias{Sugeno.integral,card.game,numeric-method}

\title{Class "card.game"}

\description{Class representing a cardinal game, i.e. a game whose values depend
only on the cardinality of subsets.}

\section{Objects from the Class}{Objects can be created by calls to the function \code{card.game}.}

\section{Slots}{
  \describe{

    \item{\code{n}:}{Object of class \code{numeric} of length 1 equal
      to the number of elements of the set on which the cardinal game is defined.}
    
    \item{\code{data}:}{Object of class \code{numeric} of length
      \code{n+1} containing the coefficients of the cardinal game. We necessarily have
      \code{data[1]=0}.}
    
  }
}

\section{Extends}{
Class \code{card.set.func}, directly.
Class \code{superclass.set.func}, by class \code{card.set.func}.
}

\section{Methods}{
  \describe{
    \item{as.game}{\code{signature(object = "card.game")} }
    \item{Choquet.integral}{\code{signature(object = "card.game", f = "numeric")} }
    \item{Sipos.integral}{\code{signature(object = "card.game", f = "numeric")} }
    \item{Sugeno.integral}{\code{signature(object = "card.game", f = "numeric")} }
  }
}

%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ }

\seealso{
  	\code{\link{game-class}}, \cr \code{\link{card.game}},
  	\cr \code{\link{as.game-methods}}, \cr \code{\link{Choquet.integral-methods}}, 
	\cr \code{\link{Sipos.integral-methods}}, \cr \code{\link{Sugeno.integral-methods}},
}

\examples{
## a cardinal game (which is a capacity)
mu <- card.game(c(0,rep(1,4)))

# the attributes of the object
mu@n
mu@data

## a conversion
as.game(mu)

## integral calculations 
f <- c(0.2,0.3,0.1,0.7)
Choquet.integral(mu,f)
Sugeno.integral(mu,f)
f <- c(0.2,-0.3,0.1,-0.7)
Sipos.integral(mu,f)
}

\keyword{classes}
