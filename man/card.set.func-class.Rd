\name{card.set.func-class}
\docType{class}
\alias{card.set.func-class}
%\alias{as.card.capacity,card.set.func-method}
%\alias{as.card.game,card.set.func-method}
%\alias{as.set.func,card.set.func-method}
%\alias{conjugate,card.set.func-method}
%\alias{interaction.indices,card.set.func-method}
%\alias{is.cardinal,card.set.func-method}
%\alias{is.kadditive,card.set.func,numeric-method}
%\alias{is.monotone,card.set.func-method}
%\alias{Mobius,card.set.func-method}
%\alias{Shapley.value,card.set.func-method}
%\alias{to.data.frame,card.set.func-method}

\title{Class "card.set.func"}

\description{Class representing a cardinal set function, i.e. whose values
depend only on the cardinality of subsets.}

\section{Objects from the Class}{Objects can be created by calls to the function \code{card.set.func}.}

\section{Slots}{
  \describe{

    \item{\code{n}:}{Object of class \code{numeric} of length 1 equal
      to the number of elements of the set on which the cardinal set
      function is defined.}
    
    \item{\code{data}:}{Object of class \code{numeric} of length
      \code{n+1} containing the coefficients of the cardinal set function.}
    
  }
}

\section{Extends}{
Class \code{superclass.set.func}, directly.
}

\section{Methods}{
  \describe{
    \item{as.card.game}{\code{signature(object = "card.set.func")}}
    
    \item{as.card.capacity}{\code{signature(object = "card.set.func")}}
    
    \item{as.set.func}{\code{signature(object = "card.set.func")} }
    
    \item{as.Mobius.card.set.func}{\code{signature(object =
	"card.set.func")} }
    
    \item{as.Mobius.set.func}{\code{signature(object = "card.set.func")}}
    
    \item{conjugate}{\code{signature(object = "card.set.func")}}
    
    \item{interaction.indices}{\code{signature(object =
	"card.set.func")}}
    
    \item{is.cardinal}{\code{signature(object = "card.set.func")}}
    
    \item{is.kadditive}{\code{signature(object = "card.set.func", k =
	"numeric")}}
    
    \item{is.monotone}{\code{signature(object = "card.set.func")}}
    
    \item{Mobius}{\code{signature(object = "card.set.func")}}
    
    \item{Shapley.value}{\code{signature(object = "card.set.func")}}
    
    \item{to.data.frame}{\code{signature(object = "card.set.func")}}
    
  }
}

\seealso{
  \code{\link{set.func-class}},
  \cr \code{\link{card.set.func}},
  \cr \code{\link{as.card.game-methods}}, 
  \cr \code{\link{as.card.capacity-methods}},
  \cr \code{\link{as.set.func-methods}},
  \cr \code{\link{as.Mobius.set.func-methods}},
  \cr \code{\link{as.Mobius.card.set.func-methods}},
  \cr \code{\link{conjugate-methods}}, 
  \cr \code{\link{interaction.indices-methods}},
  \cr \code{\link{is.cardinal-methods}}, 
  \cr \code{\link{is.kadditive-methods}},
  \cr \code{\link{is.monotone-methods}}, 
  \cr \code{\link{Mobius-methods}},
  \cr \code{\link{Shapley.value-methods}}, 
  \cr \code{\link{to.data.frame-methods}}. 
  
}

\examples{
## a cardinal set function
mu <- card.set.func(-3:2)

# the attributes of the object
mu@n
mu@data

## some conversions that cannot work
\dontrun{as.card.game(mu)}
\dontrun{as.card.capacityfunc(mu)}

## the following should work
as.set.func(mu)

## some tests
is.cardinal(mu)
is.kadditive(mu,2)
is.monotone(mu)

## some transformations
conjugate(mu)
Mobius(mu)
## let us check ...
zeta(Mobius(mu))

## summary 
Shapley.value(mu)
interaction.indices(mu)
# the same
summary(mu)

## save the set function to a file
d <- to.data.frame(mu)
\dontrun{write.table(d,"my.card.set.func.csv",sep="\t")}

# finally, some other conversions that should work
mu <- card.set.func(0:5)
as.card.game(mu)
as.card.capacity(mu)
}

\keyword{classes}
