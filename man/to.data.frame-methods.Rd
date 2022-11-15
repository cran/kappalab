\name{to.data.frame-methods}
\docType{methods}
\alias{to.data.frame}
\alias{to.data.frame-methods}
\alias{to.data.frame,Mobius.card.set.func-method}
\alias{to.data.frame,Mobius.set.func-method}
\alias{to.data.frame,card.set.func-method}
\alias{to.data.frame,set.func-method}

\title{Puts a set function under the form of a data.frame}

\description{Puts a set function under the form of a
  \code{data.frame}. This function can be used to write a set function to a
  file.}

\section{Methods}{
\describe{

\item{object = "Mobius.card.set.func" }{Returns a \code{data.frame}.}

\item{object = "Mobius.set.func" }{Returns a \code{data.frame}.}

\item{object = "card.set.func" }{Returns a \code{data.frame}.}

\item{object = "set.func" }{Returns a \code{data.frame}.}
}}

\examples{
## the Mobius representation of set function
a <- Mobius.set.func(-1:-16,4,4)

## to data.frame
d <- to.data.frame(a)
\dontrun{write.table(d,"my.set.func.csv",sep="\t")}
}

\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}
