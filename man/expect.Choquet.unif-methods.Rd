\name{expect.Choquet.unif-methods}
\docType{methods}
\alias{expect.Choquet.unif}
\alias{sd.Choquet.unif}
\alias{expect.Choquet.unif-methods}
\alias{sd.Choquet.unif-methods}
\alias{expect.Choquet.unif,game-method}
\alias{sd.Choquet.unif,game-method}

\alias{expect.Choquet.norm}
\alias{sd.Choquet.norm}
\alias{expect.Choquet.norm-methods}
\alias{sd.Choquet.norm-methods}
\alias{expect.Choquet.norm,game-method}
\alias{sd.Choquet.norm,game-method}

\alias{expect.Choquet.norm,Mobius.game-method}

\title{Expectation and standard deviation of the Choquet integral in the
uniform and normal cases}

\description{Methods for computing the expectation and standard
  deviation of the Choquet integral in the standard uniform and standard
  normal cases.}

\section{Methods}{
\describe{

   \item{object = "game"}{Returns the expectation or the standard
     deviation of the Choquet integral.}
}}

\references{
  J-L. Marichal and I. Kojadinovic (2007), \emph{The distribution of linear
  combinations of lattice polynomials from the uniform distribution},
  submitted. 
}

\seealso{
  	\code{\link{game-class}}, \code{\link{Mobius.game-class}}. 
}


\examples{

## a capacity
mu <- capacity(c(0,0.1,0.6,rep(0.9,4),1))

## the expectation and the standard deviation
## of the Choquet integral in the uniform case
expect.Choquet.unif(mu)
sd.Choquet.unif(mu)

## the same but empirically
m <- 10000
ch <- numeric(m)
for (i in 1:m) {
     f <- runif(3) 
     ch[i] <- Choquet.integral(mu,f)
}
mean(ch)
sd(ch)

## the expectation and the standard deviation
## of the Choquet integral in the normal case
expect.Choquet.norm(mu)
sd.Choquet.norm(mu)
expect.Choquet.norm(Mobius(mu))

## the same but empirically
for (i in 1:m) {
     f <- rnorm(3) 
     ch[i] <- Choquet.integral(mu,f)
}
mean(ch)
sd(ch)

}

\keyword{methods}
