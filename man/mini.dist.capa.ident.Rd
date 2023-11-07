\name{mini.dist.capa.ident}
\alias{mini.dist.capa.ident}

\title{Minimum distance capacity identification}

\description{Creates an object of class \code{Mobius.capacity} using a
  minimum distance principle. More precisely, this function determines,
  if it exists, the closest capacity to a user-given game compatible
  with a set of linear constraints. The distance can be chosen among
  three quadratic distances (see help and references hereafter). The
  problem is solved using strictly convex quadratic programming.}

\usage{
mini.dist.capa.ident(a, k, distance = "Choquet.coefficients",
A.Choquet.preorder = NULL, A.Shapley.preorder = NULL,
A.Shapley.interval = NULL, A.interaction.preorder = NULL,
A.interaction.interval = NULL, A.inter.additive.partition = NULL,
epsilon = 1e-6)
}

\arguments{
  
  \item{a}{Object of class \code{Mobius.game} containing the
    \enc{Möbius}{Mobius} transform of the game to be approached.}
  
  \item{k}{Object of class \code{numeric} imposing that the solution is at
    most a k-additive capacity (the \enc{Möbius}{Mobius} transform of subsets whose cardinal is
    superior to \code{k} vanishes).}

  \item{distance}{Object of class \code{character} indicating which
    quadratic distance is to be used in the objective function. The choice
    \code{distance = "Choquet.coefficients"} will result in the
    minimization of the average distance  between Choquet integral
    coefficients, the choice \code{distance = "binary.alternatives"}
    will result in the minimization of the average distance between binary
    alternatives, and the choice \code{distance = "global.scores"} will
    result in the minimization of the average distance between global scores.}
  
  \item{A.Choquet.preorder}{Object of class \code{matrix} containing the
    constraints relative to the preorder of the alternatives. Each line
    of the matrix corresponds to one constraint of the type "alternative
    \code{a} is preferred to alternative \code{b} with preference threshold
    \code{delta.C}". A line is structured as follows: the first \code{n}
    elements encode alternative \code{a}, the next \code{n} elements
    encode alternative \code{b}, and the last element contains the
    preference threshold \code{delta.C}.}
  
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
  
  \item{epsilon}{Object of class \code{numeric} containing the
    thresold value for the monotonicity constraints, i.e. the
    difference between the "weights" of two subsets whose cardinals
    differ exactly by 1 must be greater than \code{epsilon}.}
 
}
\details{
   The quadratic program is solved using the \code{solve.QP} function of
   the \pkg{quadprog} package.
 }

\value{
  The function returns a list structured as follows:

  \item{solution}{Object of class \code{Mobius.capacity} containing the
    \enc{Möbius}{Mobius} transform of the \code{k}-additive solution, if any.}
  
  \item{value}{Value of the objective function.}

  \item{iterations}{Information returned by \code{solve.QP}.}

  \item{iact}{Information returned by \code{solve.QP}.}
}
\references{
  K. Fujimoto and T. Murofushi (2000) \emph{Hierarchical decomposition of the
  Choquet integral}, in: Fuzzy Measures and Integrals: Theory and
  Applications, M. Grabisch, T. Murofushi, and M. Sugeno Eds, Physica
  Verlag, pages 95-103.
  
  I. Kojadinovic (2006), \emph{Quadratic objective functions for
    capacity and bi-capacity identification and approximation}, A
  Quarterly Journal of Operations Research (40R), in press.
  
}

\seealso{
  \code{\link{Mobius.capacity-class}},
  \cr \code{\link{lin.prog.capa.ident}},
  \cr \code{\link{mini.var.capa.ident}},
  \cr \code{\link{least.squares.capa.ident}},
  \cr \code{\link{heuristic.ls.capa.ident}},
  \cr \code{\link{ls.sorting.capa.ident}},
  \cr \code{\link{entropy.capa.ident}}.  
}

\examples{
## some alternatives
a <- c(18,11,18,11,11)
b <- c(18,18,11,11,11)
c <- c(11,11,18,18,11)
d <- c(18,11,11,11,18)
e <- c(11,11,18,11,18)
    
## preference threshold relative
## to the preorder of the alternatives
delta.C <- 1

## corresponding Choquet preorder constraint matrix 
Acp <- rbind(c(d,a,delta.C),
             c(a,e,delta.C),
             c(e,b,delta.C),
             c(b,c,delta.C)
            )

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
## such that I(12) = I(34)
delta.I <- 0.01
Aip <- rbind(c(1,2,3,4,-delta.I),
             c(3,4,1,2,-delta.I))

## an interaction interval constraint matrix
## i.e. -0.20 <= I(12) <= -0.15 
Aii <- rbind(c(1,2,-0.2,-0.15))

## the capacity that we want to approach
x <- runif(31)
for (i in 2:31)
    x[i] <- x[i] + x[i-1]
mu <- normalize(capacity(c(0,x)))
## and its Mobius transform
a.mu <- Mobius(mu)

## some basic checks
\dontrun{
mini.dist.capa.ident(a.mu,5)
mini.dist.capa.ident(a.mu,5,"binary.alternatives")
mini.dist.capa.ident(a.mu,5,"global.scores")
mini.dist.capa.ident(a.mu,3)
mini.dist.capa.ident(a.mu,3,"binary.alternatives")
mini.dist.capa.ident(a.mu,3,"global.scores")
}

## a minimum distance 2-additive solution
min.dist <- mini.dist.capa.ident(a.mu,2,"binary.alternatives",
                                 A.Choquet.preorder = Acp)              
m <- min.dist$solution
m

## a minimum distance 3-additive more constrained solution
min.dist2 <- mini.dist.capa.ident(a.mu,3,"global.scores",
                                   A.Choquet.preorder = Acp,
                                   A.Shapley.preorder = Asp)
m <- min.dist2$solution
m
rbind(c(a,mean(a),Choquet.integral(m,a)),
      c(b,mean(b),Choquet.integral(m,b)),
      c(c,mean(c),Choquet.integral(m,c)),
      c(d,mean(d),Choquet.integral(m,d)),
      c(e,mean(e),Choquet.integral(m,e)))
Shapley.value(m)

\dontrun{
## a minimum distance 5-additive more constrained solution
min.dist3 <- mini.dist.capa.ident(a.mu,5,
                                   A.Choquet.preorder = Acp,
                                   A.Shapley.preorder = Asp,
                                   A.Shapley.interval = Asi,
                                   A.interaction.preorder = Aip,
                                   A.interaction.interval = Aii)

m <- min.dist3$solution
m
rbind(c(a,mean(a),Choquet.integral(m,a)),
      c(b,mean(b),Choquet.integral(m,b)),
      c(c,mean(c),Choquet.integral(m,c)),
      c(d,mean(d),Choquet.integral(m,d)),
      c(e,mean(e),Choquet.integral(m,e)))
summary(m)
}

}

\keyword{math}