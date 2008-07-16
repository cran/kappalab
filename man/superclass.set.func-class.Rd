\name{superclass.set.func-class}
\docType{class}
\alias{superclass.set.func-class}
%\alias{show,superclass.set.func-method}
%\alias{summary,superclass.set.func-method}

\title{Class "superclass.set.func"}

\description{Virtual class, superclass of all \code{*.set.func} classes. Used to define common \code{show} and \code{summary} methods for \code{*.set.func} classes.}

\section{Objects from the Class}{A virtual Class: No objects may be created from it.}

\section{Slots}{
   \describe{
     \item{\code{n}:}{Object of class \code{numeric} of length 1 equal
       to the number of elements of the set on which the
       set function is defined.}
     
     \item{\code{data}:}{Object of class \code{numeric} containing the
       coefficients of the set function.}
    
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "superclass.set.func")} }
    \item{rnd}{\code{signature(x = "superclass.set.func", digits = "numeric")} }
    \item{summary}{\code{signature(object = "superclass.set.func")} }
  }
}

\seealso{
  \code{\link{show-methods}},
  \cr \code{\link{rnd-methods}},
  \cr \code{\link{summary-methods}}, 
  \cr \code{\link{set.func-class}},
  \cr \code{\link{Mobius.set.func-class}}, 
  \cr \code{\link{card.set.func-class}},
  \cr \code{\link{Mobius.card.set.func-class}}.
}
\examples{
## three set functions
mu1 <- card.set.func(-2:4)
mu2 <- set.func(c(-2,-2,-2:11/11))
mu3 <- Mobius.set.func(c(-7:6,6,6),4,4)

## print mu1
show(mu1)
## the same
mu1
## the others
mu2
mu3

## round mu2
rnd(mu2,2)

## compute indices to summarize them
summary(mu1)
summary(mu2)
summary(mu3)
}
\keyword{classes}
