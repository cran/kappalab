\name{least.squares.capa.ident}
\alias{least.squares.capa.ident}

\title{Least squares capacity identification}

\description{Creates an object of class \code{Mobius.capacity} by means of an
  approach grounded on least squares optimization. More precisely, given a set
  of data under the form: datum=(score on criterion 1, ..., score on criterion
  n, overall score), and possibly additional linear constraints expressing
  preferences, importance of criteria, etc., this function determines, if it
  exists, a capacity minimizing the sum of squared errors between overall scores
  as given by the data and the output of the Choquet integral for those data,
  and compatible with the additional linear constraints. The existence is
  ensured if no additional constraint is given. The problem is solved using
  quadratic programming.}

\usage{
least.squares.capa.ident(n, k, C, g, Integral="Choquet",
A.Shapley.preorder = NULL, A.Shapley.interval = NULL,
A.interaction.preorder = NULL, A.interaction.interval = NULL,
A.inter.additive.partition = NULL,  sigf = 7, maxiter = 40,
epsilon = 1e-6)
}

\arguments{
  \item{n}{Object of class \code{numeric} containing the
    number of elements of the set on which the object of class
    \code{Mobius.capacity} is to be defined.}
 
  \item{k}{Object of class \code{numeric} imposing that the solution is at
    most a k-additive capacity (the \enc{Möbius}{Mobius} transform of subsets whose cardinal is
    superior to \code{k} vanishes).}

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
    
   \item{A.Shapley.preorder}{Object of class \code{matrix} containing the
    constraints relative to the preorder of the criteria. Each line
    of this 3-column matrix corresponds to one constraint of the type
    "the Shapley importance index of criterion \code{i} is greater than
    the Shapley importance index of criterion \code{j} with preference threshold
    \code{delta.S}". A line is structured as follows: the first element
    encodes \code{i}, the second \code{j}, and the third element contains
    the preference threshold \code{delta.S}.}
  
  \item{A.Shapley.interval}{Object of class \code{matrix} containing the
    constraints relative to the quantitative importance of the
    criteria. Each line of this 3-column matrix corresponds to one
    constraint of the type "the Shapley importance index of criterion
    \code{i} lies in the interval \code{[a,b]}". The interval
    \code{[a,b]} has to be included in \code{[0,1]}. A line of the
    matrix is structured as follows: the first element encodes \code{i},
    the second \code{a}, and the third \code{b}.}
  
  \item{A.interaction.preorder}{Object of class \code{matrix}
    containing the constraints relative to the preorder of the pairs of
    criteria in terms of the Shapley interaction index. Each line of this 5-column matrix
    corresponds to one constraint of the type "the Shapley interaction
    index of the pair \code{ij} of criteria is greater than the Shapley interaction
    index of the pair \code{kl} of criteria with preference threshold \code{delta.I}".
    A line is structured as follows: the first two elements encode
    \code{ij}, the second two \code{kl}, and the fifth element contains
    the preference threshold \code{delta.I}.}
  
  \item{A.interaction.interval}{Object of class \code{matrix}
    containing the constraints relative to the type and the magnitude of
    the Shapley interaction index for pairs of criteria. Each line of
    this 4-column matrix corresponds to one constraint of the type
    "the Shapley interaction index of the pair \code{ij} of criteria
    lies in the interval \code{[a,b]}". The interval \code{[a,b]} has to
    be included in \code{[-1,1]}. A line is structured as follows: the first two elements encode
    \code{ij}, the third element encodes \code{a}, and the fourth element
    encodes \code{b}.}

  \item{A.inter.additive.partition}{Object of class \code{numeric}
    encoding a partition of the set of criteria imposing that there be
    no interactions among criteria belonging to different classes 
    of the partition. The partition is to be given under the form of a
    vector of integers from \code{{1,\dots,n}} of length \code{n} such
    that two criteria belonging to the same class are "marked" by the
    same integer. For instance, the partition \code{{{1,3},{2,4},{5}}} can
    be encoded as \code{c(1,2,1,2,3)}. See Fujimoto and Murofushi (2000)
    for more details on the concept of mu-inter-additive partition.}

  \item{sigf}{Precision (default: 7 significant figures). Parameter to
    be passed to the \code{ipop} function (quadratic programming)
    of the \pkg{kernlab} package.}

  \item{maxiter}{Maximum number of iterations. Parameter to
    be passed to the \code{ipop} function (quadratic programming)
    of the \pkg{kernlab} package.}
  
  \item{epsilon}{Object of class \code{numeric} containing the
    threshold value for the monotonicity constraints, i.e. the
    difference between the "weights" of two subsets whose cardinals
    differ exactly by 1 must be greater than \code{epsilon}.}
}

\details{
   The quadratic program is solved using the \code{ipop} function of
   the \pkg{kernlab} package.
 }
 
\value{
  The function returns a list structured as follows:

  \item{solution}{Object of class
    \code{Mobius.capacity} containing the \enc{Möbius}{Mobius} transform of the
    \code{k}-additive solution, if any.}
  
  \item{dual}{The dual solution of the problem.}

  \item{how}{Character string describing the type of convergence.}

  \item{residuals}{Differences between the provided global evaluations
    and those returned by the obtained model.}
}

\references{
  K. Fujimoto and T. Murofushi (2000) \emph{Hierarchical decomposition of the
  Choquet integral}, in: Fuzzy Measures and Integrals: Theory and
  Applications, M. Grabisch, T. Murofushi, and M. Sugeno Eds, Physica
  Verlag, pages 95-103.
  
  M. Grabisch, H.T. Nguyen and E.A. Walker (1995), \emph{Fundamentals of
    uncertainty calculi with applications to fuzzy inference}, Kluwer
  Academic, Dordrecht.

  M. Grabisch and M. Roubens (2000), \emph{Application of the Choquet Integral in
  Multicriteria Decision Making}, in: Fuzzy Measures and Integrals: Theory and
  Applications, M. Grabisch, T. Murofushi, and M. Sugeno Eds, Physica Verlag,
  pages 415-434.

  P. Miranda and M. Grabisch (1999), \emph{Optimization issues for fuzzy measures},
  International Journal of Fuzziness and Knowledge-based Systems 7:6, pages
  545-560.
}

\seealso{
  \code{\link{Mobius.capacity-class}},
  \cr \code{\link{heuristic.ls.capa.ident}},
  \cr \code{\link{lin.prog.capa.ident}},
  \cr \code{\link{mini.var.capa.ident}},
  \cr \code{\link{mini.dist.capa.ident}},
  \cr \code{\link{ls.sorting.capa.ident}},
  \cr \code{\link{entropy.capa.ident}}.  
}

\examples{

## the number of data
n.d <- 20

## a randomly generated 5-criteria matrix
C <- matrix(rnorm(5*n.d,10,2),n.d,5)


## the corresponding global scores
g <- numeric(n.d)
mu <- capacity(c(0:29,29,29)/29)
for (i in 1:n.d)
  g[i] <- Choquet.integral(mu,C[i,])

\dontrun{
## the full solution 
lsc <- least.squares.capa.ident(5,5,C,g)
a <- lsc$solution
a
mu.sol <- zeta(a)

## the difference between mu and mu.sol
mu@data - mu.sol@data

## the residuals
lsc$residuals

## the mean square error
mean(lsc$residuals^2)

## a 3-additive solution 
lsc <- least.squares.capa.ident(5,3,C,g)
a <- lsc$solution
mu.sol <- zeta(a)
mu@data - mu.sol@data
lsc$residuals
}



## a similar example based on the Sipos integral

## a randomly generated 5-criteria matrix
C <- matrix(rnorm(5*n.d,0,2),n.d,5)

## the corresponding global scores
g <- numeric(n.d)
mu <- capacity(c(0:29,29,29)/29)
for (i in 1:n.d)
  g[i] <- Sipos.integral(mu,C[i,])

\dontrun{
## the full solution 
lsc <- least.squares.capa.ident(5,5,C,g,Integral = "Sipos")
a <- lsc$solution
mu.sol <- zeta(a)
mu@data - mu.sol@data
lsc$residuals

## a 3-additive solution 
lsc <- least.squares.capa.ident(5,3,C,g,Integral = "Sipos")
a <- lsc$solution
mu.sol <- zeta(a)
mu@data - mu.sol@data
lsc$residuals
}



## additional constraints

## a Shapley preorder constraint matrix
## Sh(1) - Sh(2) >= -delta.S
## Sh(2) - Sh(1) >= -delta.S
## Sh(3) - Sh(4) >= -delta.S
## Sh(4) - Sh(3) >= -delta.S
## i.e. criteria 1,2 and criteria 3,4
## should have the same global importances
delta.S <- 0.01    
Asp <- rbind(c(1,2,-delta.S),
             c(2,1,-delta.S),
             c(3,4,-delta.S),
             c(4,3,-delta.S)
            )

## a Shapley interval constraint matrix
## 0.3 <= Sh(1) <= 0.9 
Asi <- rbind(c(1,0.3,0.9))


## an interaction preorder constraint matrix
## such that I(12) = I(45)
delta.I <- 0.01
Aip <- rbind(c(1,2,4,5,-delta.I),
             c(4,5,1,2,-delta.I))

## an interaction interval constraint matrix
## i.e. -0.20 <= I(12) <= -0.15 
delta.I <- 0.01
Aii <- rbind(c(1,2,-0.2,-0.15))

## an inter-additive partition constraint
## criteria 1,2,3 and criteria 4,5 are independent 
Aiap <- c(1,1,1,2,2)





## a more constrained solution

lsc <- least.squares.capa.ident(5,5,C,g,Integral = "Sipos",
                                 A.Shapley.preorder = Asp,
                                 A.Shapley.interval = Asi,
                                 A.interaction.preorder = Aip,
                                 A.interaction.interval = Aii,
                                 A.inter.additive.partition = Aiap,
                                 sigf = 5)

a <- lsc$solution
mu.sol <- zeta(a)
mu@data - mu.sol@data
lsc$residuals
summary(a)
}

\keyword{math}