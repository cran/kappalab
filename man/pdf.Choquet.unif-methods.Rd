\name{pdf.Choquet.unif-methods}
\docType{methods}
\alias{pdf.Choquet.unif}
\alias{cdf.Choquet.unif}
\alias{pdf.Choquet.unif-methods}
\alias{cdf.Choquet.unif-methods}
\alias{pdf.Choquet.unif,game,numeric-method}
\alias{cdf.Choquet.unif,game,numeric-method}

\title{Distribution of the Choquet integral for evaluations uniformly
  distributed on the unit hypercube}

\description{Methods for computing the probability density and cumulative
  distribution functions of the Choquet integral with respect to a game
  for evaluations uniformly distributed on the unit hypercube.}


\section{Methods}{
\describe{

   \item{object = "game", y = "numeric" }{Returns the value of the
     p.d.f. or the c.d.f. at \code{y}.}
}}

\references{
  J-L. Marichal and I. Kojadinovic (2007), \emph{The distribution of linear
  combinations of lattice polynomials from the uniform distribution},
  submitted. 
}

\seealso{
  	\code{\link{game-class}}. 
}


\examples{

## a capacity
mu <- capacity(c(0,0.1,0.6,rep(0.9,4),1))
## the cdf of the Choquet integral at 0.7
cdf.Choquet.unif(mu,0.7)

## the same but empirically
m <- 10000
ch <- numeric(m)
for (i in 1:m) {
     f <- runif(3) 
     ch[i] <- Choquet.integral(mu,f)
}
sum(ifelse(ch<=0.7,1,0))/m
}

\keyword{methods}
