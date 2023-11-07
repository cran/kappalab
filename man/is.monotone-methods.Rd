\name{is.monotone-methods}
\docType{methods}
\alias{is.monotone}
\alias{is.monotone-methods}
\alias{is.monotone,Mobius.set.func-method}
\alias{is.monotone,card.set.func-method}
\alias{is.monotone,set.func-method}

\title{Test method}

\description{Tests whether a set function is monotone with respect to set
  inclusion. The set function
  can be given either under the form of an object of class
  \code{set.func}, \code{card.set.func} or \code{Mobius.set.func}.
}

\section{Methods}{
\describe{

  \item{object = "Mobius.set.func", verbose = "logical", epsilon =
    "numeric" }{
    Returns an object of class \code{logical}. If \code{verbose=TRUE},
    displays the violated monotonicity constraints, if any.}

  \item{object = "card.set.func", verbose = "logical", epsilon =
    "numeric" }{
    Returns an object of class  \code{logical}. If
    \code{verbose=TRUE}, displays the violated monotonicity
    constraints, if any.}

\item{object = "set.func", verbose = "logical", epsilon = "numeric" }{
  Returns an object of class \code{logical}. If \code{verbose=TRUE},
  displays the violated monotonicity constraints, if any.}
}}

\details{
  For objects of class \code{set.func} or \code{card.set.func}, the
  monotonicity constraints are considered to be satisfied
  (cf. references hereafter) if the following inequalities are satisfied
  \deqn{\mu(S \cup i) - \mu(S) \ge -epsilon}{mu(S U i) - mu(S) >=
    -epsilon}
  for all \eqn{S} and all \eqn{i}.
  For objects of class \code{Mobius.set.func}, it is
  required that a similar condition with respect to the \enc{Möbius}{Mobius}
  representation be satisfied (cf. references hereafter).
}

\references{
  A. Chateauneuf and J-Y. Jaffray (1989), \emph{Some characterizations of
    lower probabilities and other monotone capacities through the use of
    \enc{Möbius}{Mobius} inversion}, Mathematical Social Sciences 17:3, pages
  263--283.

  M. Grabisch (2000), \emph{The interaction and \enc{Möbius}{Mobius} representations of fuzzy
    measures on finites spaces, k-additive measures: a survey}, in:
  Fuzzy Measures and Integrals: Theory and Applications, M. Grabisch,
  T. Murofushi, and M. Sugeno Eds, Physica Verlag, pages 70-93.
}

\seealso{
  	\code{\link{Mobius.set.func-class}}, 
	\cr \code{\link{card.set.func-class}}, 
	\cr \code{\link{set.func-class}}. 
}
      
\examples{
## a monotone set function
mu <- set.func(c(0,1,1,1,2,2,2,3))
mu
is.monotone(mu)

## the Mobius representation of a monotone set function
a <- Mobius.set.func(c(0,1,2,1,3,1,2,1,2,3,1),4,2)
is.monotone(a)

## non-monotone examples
mu <- set.func(c(0,-7:7))
is.monotone(mu,verbose=TRUE)
a <- Mobius(mu)
is.monotone(a,verbose=TRUE)
}

\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}
