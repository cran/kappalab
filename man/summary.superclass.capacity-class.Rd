\name{summary.superclass.capacity-class}
\docType{class}
\alias{summary.superclass.capacity-class}
%\alias{show,summary.superclass.capacity-method}

\title{Class "summary.superclass.capacity"}

\description{Class containing most of the indices that can be used to summarize an object of class \code{superclass.capacity}.}

\section{Objects from the Class}{Objects are created by calls of the form \code{summary(...)} on objects of class \code{superclass.capacity}.}

\section{Slots}{
  \describe{
    \item{\code{Shapley.value}:}{Object of class \code{numeric} containing the Shapley value of a set function.}
    \item{\code{interaction.indices}:}{Object of class \code{matrix} containing the Shapley interaction indices of a set function.}
    \item{\code{orness}:}{Object of class \code{numeric} containg the orness degree of the Choquet integral corresponding to the underlying capacity.}
    \item{\code{veto}:}{Object of class \code{numeric} containing the veto indices of the Choquet integral corresponding to the underlying capacity.}
    \item{\code{favor}:}{Object of class \code{numeric} containing the favor indices of the Choquet integral corresponding to the underlying capacity.}
    \item{\code{variance}:}{Object of class \code{numeric} containing the variance of the underlying capacity.} 
    \item{\code{entropy}:}{Object of class \code{numeric} containing the variance of the underlying capacity.}
  }
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "summary.superclass.capacity")}: ... }
  }
}

\seealso{
	\code{\link{superclass.capacity-class}}, \cr \code{\link{summary-methods}},
	\cr \code{\link{Shapley.value-methods}}, \cr \code{\link{interaction.indices-methods}}, 
	\cr \code{\link{orness-methods}}, \cr \code{\link{veto-methods}},
	\cr \code{\link{favor-methods}}, \cr \code{\link{variance-methods}},
	\cr \code{\link{entropy-methods}}.
}

\examples{
## a capacity
mu <- capacity(c(0:13,13,13)/13)
a <- Mobius(mu)

## its summary
summary(mu)
summary(a)
}

\keyword{classes}
