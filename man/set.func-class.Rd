\name{set.func-class}
\docType{class}
\alias{set.func-class}
%\alias{as.capacity,set.func-method}
%\alias{as.card.set.func,set.func-method}
%\alias{as.game,set.func-method}
%\alias{dual,set.func-method}
%\alias{interaction.indices,set.func-method}
%\alias{is.cardinal,set.func-method}
%\alias{is.kadditive,set.func,numeric-method}
%\alias{is.monotone,set.func-method}
%\alias{k.truncate.Mobius,set.func,numeric-method}
%\alias{Mobius,set.func-method}
%\alias{Shapley.value,set.func-method}
%\alias{to.data.frame,set.func-method}

\title{Class "set.func"}

\description{Class representing a set function.}

\section{Objects from the Class}{Objects can be created by calls to the function \code{set.func}.}

\section{Slots}{
  \describe{
    
    \item{\code{n}:}{Object of class \code{numeric} of length 1 equal to
      the number of elements of the set on which the set function is defined.}
    
    
    \item{\code{subsets}:}{Object of class \code{numeric} of length
      \code{2^n} containing the power set of the underlying set in
      "natural" order. The subsets are coded as integers.}
    
    \item{\code{data}:}{Object of class \code{numeric} of length
      \code{2^n} containing the coefficients of the set function in binary order.}
    
  }
}

\section{Extends}{
Class \code{superclass.set.func}, directly.
}

\section{Methods}{
  \describe{

    \item{show}{\code{signature(object = "set.func")} }
    
    \item{as.game}{\code{signature(object = "set.func")}}

    \item{as.capacity}{\code{signature(object = "set.func")}}

    \item{as.card.set.func}{\code{signature(object = "set.func")}}

    \item{as.Mobius.set.func}{\code{signature(object = "set.func")}}

    \item{as.Mobius.card.set.func}{\code{signature(object = "set.func")}}

    \item{dual}{\code{signature(object = "set.func")}}

    \item{interaction.indices}{\code{signature(object = "set.func")}}

    \item{is.cardinal}{\code{signature(object = "set.func")}}

    \item{is.kadditive}{\code{signature(object = "set.func", k = "numeric")}}

    \item{is.monotone}{\code{signature(object = "set.func")}}

    \item{k.truncate.Mobius}{\code{signature(object = "set.func", k =
	"numeric")}}
    
    \item{Mobius}{\code{signature(object = "set.func")}}
    
    \item{Shapley.value}{\code{signature(object = "set.func")}}
    
    \item{to.data.frame}{\code{signature(object = "set.func")}}
    
  }
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ }

\seealso{
  \code{\link{set.func}},
  \cr \code{\link{as.game-methods}},	
  \cr \code{\link{as.capacity-methods}},
  \cr \code{\link{as.card.set.func-methods}}, 
  \cr \code{\link{as.Mobius.set.func-methods}},
  \cr \code{\link{as.Mobius.card.set.func-methods}},
  \cr \code{\link{dual-methods}}, 
  \cr \code{\link{interaction.indices-methods}},
  \cr \code{\link{is.cardinal-methods}}, 
  \cr \code{\link{is.kadditive-methods}},
  \cr \code{\link{is.monotone-methods}},
  \cr \code{\link{k.truncate.Mobius-methods}},
  \cr \code{\link{Mobius-methods}}, 
  \cr \code{\link{Shapley.value-methods}},
  \cr \code{\link{to.data.frame-methods}}. 
}

\examples{
## a set function
mu <- set.func(c(1:8,8:1)/8)

## the attributes of the object
mu@n
mu@data
mu@subsets

## some conversions that cannot work
\dontrun{as.game(mu)}
\dontrun{as.capacity(mu)}
\dontrun{as.card.set.func(mu)}

## some tests
is.cardinal(mu)
is.kadditive(mu,2)
is.monotone(mu)

## some transformations
dual(mu)
Mobius(mu)
k.truncate.Mobius(mu,2)

## summary 
Shapley.value(mu)
interaction.indices(mu)
# the same
summary(mu)

## save the set function to a file
d <- to.data.frame(mu)
write.table(d,"my.set.func.csv",sep="\t")

# finally, some conversions that should work
mu <- set.func(c(0,1,1,1,2,2,2,3))
as.game(mu)
as.capacity(mu)
as.card.set.func(mu)
}

\keyword{classes}
