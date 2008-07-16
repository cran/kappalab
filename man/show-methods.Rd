\name{show-methods}
\docType{methods}
\alias{show-methods}
%\alias{show,ANY-method}
%\alias{show,traceable-method}
%\alias{show,ObjectsWithPackage-method}
%\alias{show,MethodDefinition-method}
%\alias{show,MethodWithNext-method}
%\alias{show,genericFunction-method}
%\alias{show,classRepresentation-method}
\alias{show,superclass.set.func-method}
\alias{show,set.func-method}
\alias{show,Mobius.set.func-method}
\alias{show,summary.superclass.set.func-method}
\alias{show,summary.superclass.capacity-method}

\title{Methods for Function show in Package `kappalab'}

\description{Displays the contents of an object.}

\section{Methods}{
\describe{

  \item{object = "superclass.set.func" }{Generic \code{show} method for \code{superclass.set.func}.}

  \item{object = "set.func" }{Displays a \code{set.func} object.}

  \item{object = "Mobius.set.func" }{Displays a \code{Mobius.set.func} object.}

  \item{object = "summary.superclass.set.func" }{Displays a \code{summary.superclass.set.func} object.}

  \item{object = "summary.superclass.capacity" }{Displays a
    \code{summary.superclass.capacity} object.}
  }

}

\seealso{
  \code{\link{superclass.set.func-class}},
  \cr \code{\link{set.func-class}},
  \cr \code{\link{Mobius.set.func-class}},
  \cr \code{\link{summary.superclass.set.func-class}},
  \cr \code{\link{summary.superclass.capacity-class}}. 
}

\examples{
## a set function
mu <- set.func(0:15/15)
show(mu)
## the same
mu

## a Mobius transform
a <- Mobius.set.func(0:10,4,2)
show(a)
a

## a cardinal capacity
mu <- lower.capacity(5)
show(mu)
mu
}

\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}
