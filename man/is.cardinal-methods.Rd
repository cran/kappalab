\name{is.cardinal-methods}
\docType{methods}
\alias{is.cardinal}
\alias{is.cardinal-methods}
\alias{is.cardinal,Mobius.set.func-method}
\alias{is.cardinal,card.set.func-method}
\alias{is.cardinal,set.func-method}

\title{Test method}

\description{Tests whether a set function is cardinal, i.e., if its values
  depend only on the cardinality of subsets. The set function
  can be given either under the form of an object of class
  \code{set.func}, \code{card.set.func} or \code{Mobius.set.func}.}

\section{Methods}{
\describe{

\item{object = "Mobius.set.func" }{Returns an object of class \code{logical}.}

\item{object = "card.set.func" }{Returns an object of class \code{logical}.}

\item{object = "set.func" }{Returns an object of class \code{logical}.}
}}

\examples{
is.cardinal(set.func(-7:8))
is.cardinal(uniform.capacity(8))
is.cardinal(Mobius.game(0:10,4,2))
}

\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}
