\name{Mobius.set.func-class}
\docType{class}
\alias{Mobius.set.func-class}
%\alias{as.Mobius.capacity,Mobius.set.func-method}
%\alias{as.Mobius.card.set.func,Mobius.set.func-method}
%\alias{as.Mobius.game,Mobius.set.func-method}
%\alias{interaction.indices,Mobius.set.func-method}
%\alias{is.cardinal,Mobius.set.func-method}
%\alias{is.kadditive,Mobius.set.func,numeric-method}
%\alias{is.monotone,Mobius.set.func-method}
%\alias{k.truncate.Mobius,Mobius.set.func,numeric-method}
%\alias{Shapley.value,Mobius.set.func-method}
%\alias{to.data.frame,Mobius.set.func-method}
%\alias{zeta,Mobius.set.func-method}

\title{Class "Mobius.set.func"}

\description{Class representing the \enc{Möbius}{Mobius} transform of a set function.}

\section{Objects from the Class}{Objects can be created by calls to the
  function \code{Mobius.set.func}.
}

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
      coefficients of a truncated \enc{Möbius}{Mobius} transform of a set function in "natural" order.}
    
  }
}

\section{Extends}{
  Class \code{superclass.set.func}, directly.
}

\section{Methods}{
  \describe{

    \item{show}{\code{signature(object = "Mobius.set.func")} }

    \item{as.Mobius.card.set.func}{\code{signature(object =
	"Mobius.set.func")} }

    \item{as.card.set.func}{\code{signature(object =
	"Mobius.set.func")} }

    \item{as.set.func}{\code{signature(object = "Mobius.set.func")}}

    \item{as.Mobius.game}{\code{signature(object = "Mobius.set.func")}}

    \item{as.Mobius.capacity}{\code{signature(object = "Mobius.set.func")}}

    \item{interaction.indices}{\code{signature(object = "Mobius.set.func")}}

    \item{is.cardinal}{\code{signature(object = "Mobius.set.func")}}

    \item{is.kadditive}{\code{signature(object = "Mobius.set.func", k =
	"numeric")} }
    
    \item{is.monotone}{\code{signature(object = "Mobius.set.func")}}

    \item{k.truncate.Mobius}{\code{signature(object = "Mobius.set.func",
	k = "numeric")}}
    
    \item{Shapley.value}{\code{signature(object = "Mobius.set.func")}}
    
    \item{to.data.frame}{\code{signature(object = "Mobius.set.func")}}
    
    \item{zeta}{\code{signature(object = "Mobius.set.func")}}
    
  }
}

%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ }

\seealso{
	\code{\link{set.func-class}}, 
	\cr \code{\link{Mobius.set.func}},
	\cr \code{\link{as.Mobius.card.set.func-methods}},
	\cr \code{\link{as.card.set.func-methods}},
	\cr \code{\link{as.set.func-methods}},  
	\cr \code{\link{as.Mobius.game-methods}},
	\cr \code{\link{as.Mobius.capacity-methods}}, 
	\cr \code{\link{interaction.indices-methods}}, 
	\cr \code{\link{is.cardinal-methods}}, 
	\cr \code{\link{is.kadditive-methods}},
	\cr \code{\link{is.monotone-methods}}, 
	\cr \code{\link{k.truncate.Mobius-methods}}, 
	\cr \code{\link{Shapley.value-methods}}, 
	\cr \code{\link{to.data.frame-methods}},
	\cr \code{\link{zeta-methods}}.
}

\examples{
## the Mobius transform of a set function directly
a <- Mobius.set.func(1:16,4,4)

## the attributes of the object
a@n
a@k
a@data
a@subsets

## a set function
mu <- set.func(7:-8)
## and its Mobius transform
a <- Mobius(mu)

## some conversions that cannot work
## as.game(a)
## as.capacity(a)
## as.card.set.func(a)

## some tests
is.cardinal(a)
is.kadditive(a,2)
is.monotone(a)

## some transformations
zeta(a)
k.truncate.Mobius(a,2)

## summary 
Shapley.value(a)
interaction.indices(a)
# the same
summary(a)

## save the Mobius transform to a file
d <- to.data.frame(a)
\dontrun{write.table(d,"my.Mobius.set.func.csv",sep="\t")}

# finally, some conversions that should work
mu <- set.func(c(0,1,1,1,2,2,2,3))
a <- Mobius(mu)
as.Mobius.game(a)
as.Mobius.capacity(a)
as.Mobius.card.set.func(a)
}

\keyword{classes}
