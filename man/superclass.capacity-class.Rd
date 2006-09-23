\name{superclass.capacity-class}
\docType{class}
\alias{superclass.capacity-class}
%\alias{summary,superclass.capacity-method}

\title{Class "superclass.capacity"}

\description{Virtual class, superclass of all \code{*.capacity} classes. Used to define a common \code{summary} method for \code{*.capacity} classes. Contains no slots.}

\section{Objects from the Class}{A virtual Class: No objects may be created from it.}

\section{Methods}{
  \describe{
    \item{summary}{\code{signature(object = "superclass.capacity")}}
  }
}

%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ }

\seealso{
	\code{\link{summary-methods}},
  	\cr \code{\link{capacity-class}}, \cr \code{\link{Mobius.capacity-class}}, 
	\cr \code{\link{card.capacity-class}}.
}

\examples{
## three capacities
mu1 <- uniform.capacity(5)
mu2 <- capacity(c(0,0,0:13))
a <- Mobius(mu2)

## compute indices to summarize them
summary(mu1)
summary(mu2)
summary(a)
}

\keyword{classes}
