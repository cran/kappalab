\name{summary.superclass.set.func-class}
\docType{class}
\alias{summary.superclass.set.func-class}
%\alias{show,summary.superclass.set.func-method}

\title{Class "summary.superclass.set.func"}

\description{Class containing most of the indices that can be used to summarize an object of class \code{superclass.set.func} and not of class \code{superclass.capacity}.}

\section{Objects from the Class}{Objects are created by calls of the form \code{summary(...)} on objects of class \code{superclass.set.func} and not of class \code{superclass.capacity}.}

\section{Slots}{
  \describe{
    \item{\code{Shapley.value}:}{Object of class \code{numeric} containing the Shapley value of a set function.}
    \item{\code{interaction.indices}:}{Object of class \code{matrix} containing the Shapley interaction indices of a set function.}
  }
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "summary.superclass.set.func")}}
  }
}

\seealso{
  	\code{\link{superclass.set.func-class}}, \cr \code{\link{superclass.capacity-class}},
	\cr \code{\link{summary-methods}},
	\cr \code{\link{show-methods}}, \cr \code{\link{Shapley.value-methods}}, 
	\cr \code{\link{interaction.indices-methods}}.
}

\examples{
## a capacity
mu <- set.func(c(0:13,13,13)/13)
a <- Mobius(mu)

## its summary
summary(mu)
summary(a)
}

\keyword{classes}
