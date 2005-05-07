\name{set.func}
\alias{set.func}
\alias{game}
\alias{capacity}

\title{Create objects of class "set.func", "game", or "capacity".}

\description{These functions create objects of class \code{set.func},
  \code{game}, or \code{capacity} from objects of class \code{numeric}.}

\usage{
set.func(object)
game(object)
capacity(object)
}

\arguments{
  \item{object}{An object of class \code{numeric} containing the
    coefficients of the set function in "natural" order. Its length has
    to be a power of 2.}
}

\value{
Return objects of class \code{set.func}, \code{game}, or \code{capacity}.
}

\seealso{
	\code{\link{capacity-class}}, 
	\cr \code{\link{game-class}},
	\cr \code{\link{set.func-class}}.
}

\examples{
set.func(-12:3)
game(c(0,rep(-1,15)))
capacity(0:15)
}
\keyword{math}