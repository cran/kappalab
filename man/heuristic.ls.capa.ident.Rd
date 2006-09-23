\name{heuristic.ls.capa.ident}
\alias{heuristic.ls.capa.ident}

\title{Heuristic least squares capacity identification}

\description{Creates an object of class \code{capacity} or \code{game}
  by means of the heuristic algorithm proposed by Grabisch (1995). More
  precisely, given a set of data under the form:
  datum=(score on criterion 1, ..., score on criterion
  n, overall score), this function heuristically determines a capacity
  that should be close to minimizing the sum of squared errors between overall scores
  as given by the data and the output of the Choquet integral for those
  data. See reference hereafter for more details.}

\usage{
heuristic.ls.capa.ident(n, mu, C, g, Integral="Choquet", maxiter = 500,
                        alpha = 0.01, epsilon = 1e-6) 
}

\arguments{
  \item{n}{Object of class \code{numeric} containing the
    number of elements of the set on which the object of class
    \code{Mobius.capacity} is to be defined.}

  \item{mu}{Object of class \code{capacity} containing the initial capacity. It
  should be initialized to the equidistributed capacity.}
 
  \item{C}{Object of class \code{matrix} containing the
    \code{n}-column criteria matrix. Each line of this matrix
    corresponds to a vector of the form: (score on criterion 1, ..., score on
    criterion n).}

  \item{g}{Object of class \code{numeric} containg the global
    scores associated with the vectors given in the criteria matrix.}

  \item{Integral}{Object of class \code{character} indicating whether the
    model is based on the asymmetric Choquet integral (\code{Integral =
      "Choquet"}) or the symmetric Choquet integral (\code{Integral =
      "Sipos"}).}

  \item{maxiter}{Maximum number of iterations. }
  
  \item{alpha}{Object of class \code{numeric} containing the coefficient of the
  gradient (multiplicative coefficient of the partial derivative). Its value
  lies in [0,1]. The higher the value, the larger the move of the capacity at
  each iteration.}

  \item{epsilon}{Object of class \code{numeric} used for defining a stopping
  criterion. Iterations are stopped if one of the following cases happen: 1)
  itmax iterations have been performed 2) the error criterion is growing up 3)
  normax < epsilon, where normax is the maximum absolute normalized model error,
  i.e., the maximum over all data of the quantity |e|/outmax, where e=desired
  output - model output, and outmax is the highest value of the absolute value
  of the model output}
}

\details{
  The algorithm is explained in details in the reference hereafter.
 }
 
\value{
  The function returns a list structured as follows:

  \item{solution}{Object of class \code{capacity} if the solution is
    montone or of class \code{game} otherwise.}
  
  \item{n.iter}{Number of iterations taken by the algorithm.}

  \item{residuals}{Differences between the provided global evaluations
    and those returned by the obtained model.}

  \item{mse}{Mean square error between the provided global evaluations
    and those returned by the obtained model.}

}

\references{
  M. Grabisch (1995), \emph{A new algorithm for identifying fuzzy measures and
    its application to pattern recognition}, Int. Joint Conf. of the 4th
  IEEE Int. Conf. on Fuzzy Systems and the 2nd Int. Fuzzy Engineering
  Symposium, Yokohama, Japan, 145-150.
}

\seealso{
  \code{\link{capacity-class}},
  \cr \code{\link{least.squares.capa.ident}},
  \cr \code{\link{lin.prog.capa.ident}},
  \cr \code{\link{mini.var.capa.ident}},
  \cr \code{\link{mini.dist.capa.ident}},
  \cr \code{\link{ls.sorting.capa.ident}},
  \cr \code{\link{ls.ranking.capa.ident}},
  \cr \code{\link{entropy.capa.ident}}.  
}

\examples{

## number of criteria
n <- 4

## the number of alternatives
n.a <- 1000

## a randomly generated 5-criteria matrix
C <- matrix(rnorm(n*n.a,10,2),n.a,n)

## the corresponding global scores
g <- numeric(n.a)

## generate a random capacity
x <- runif(2^n-1)
for (i in 2:(2^n-1))
    x[i] <- x[i] + x[i-1]
mu <- normalize(capacity(c(0,x)))
for (i in 1:n.a)
  g[i] <- Choquet.integral(mu,C[i,])

## the initial capacity
## here the uniform capacity
mu.in <- as.capacity(uniform.capacity(n))

## the solution 
hlsc <- heuristic.ls.capa.ident(n,mu.in,C,g)
mu.sol <- hlsc$solution

## the difference between mu and mu.sol
mu@data - mu.sol@data

hlsc
}

\keyword{math}