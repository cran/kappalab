\name{entropy.capa.ident}
\alias{entropy.capa.ident}

\title{Unsupervised identification of a capacity from profiles}

\description{This function estimates a capacity using as argument a set
  of data under the form: datum=(score on attribute 1, ..., score on attribute
  n). The approach roughly consists in replacing the subjective notion
  of \emph{importance} of a subset of attributes by that of
  \emph{information content} of a subset of attributes, which is estimated
  from the data by means of a parametric entropy measure. For
  more details, see the references hereafter.
}


\usage{
entropy.capa.ident(d, entropy = "renyi", parameter = 1)
}

\arguments{
  \item{d}{An object of class \code{data.frame} containing the
    discretized data. Each column of the \code{data.frame} must
    be a \code{factor}. Each line corresponds to a datum.}
  
  \item{entropy}{An object of class \code{character} containg the name
    of the parametric entropy measure to be used for the estimation. The
  allowed values are "renyi" and "havrda.charvat".}
  
  \item{parameter}{An object of class \code{numeric} containing he value
    of the parameter of the choosen entropy. The parameter value must be
  a positive real number. If equal to 1, the Shannon entropy is used.}
  
}

\value{Returns an object of class \code{capacity}.}

\references{
  I. Kojadinovic (2004), \emph{Estimation of the weights of interacting
    criteria from the set of profiles by means of information-theoretic
    functionals}, European Journal of Operational Research 155:3, pages 741-751.
  
  I. Kojadinovic (2005), \emph{Unusupervised aggregation of commensurate
    correlated attributes by means of the Choquet integral and entropy
    functionals}, International Journal of Intelligent Systems, in press.
}
  
\seealso{
  \code{\link{capacity-class}},
  \cr \code{\link{least.squares.capa.ident}},
  \cr \code{\link{mini.var.capa.ident}}.  
}

\examples{
## a set of randomly generated data
## for instance, marks on a [0,20] scale
p <- data.frame(matrix(runif(500,0,20),100,5))
names(p) <- c("Stat","Prob","Alg","Cal","Eng")

## discretization
p[p <= 5] <- 1  
p[p > 5 & p <= 10] <- 2 
p[p > 10 & p <= 15] <- 3 
p[p > 15] <- 4

d <- data.frame(factor(p[[1]]),
                factor(p[[2]]),
                factor(p[[3]]),
                factor(p[[4]]),
                factor(p[[5]]))

## associated unsupervised capacity
mu <- entropy.capa.ident(d)
mu
}
\keyword{math}