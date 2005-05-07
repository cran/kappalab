\name{ls.sorting.capa.ident}
\alias{ls.sorting.capa.ident}

\title{Least squares capacity identification in the framework of a sorting procedure}

\description{
  Sorting alternatives means assigning each alternative to a predefined
  ordered class. The aim of the implemented method is to model a given classification (sorting)
  of the alternatives by means of a Choquet integral. The result of the function
  is an object of class \code{Mobius.capacity}. This function (in
  combination with \code{ls.sorting.treatment}) is an implementation of
  the TOMASO method; see Meyer and Roubens (2005). The input data are
  given under the form of a set of alternatives and associated classes,
  each alternative being described according to a set of criteria. These
  well-known alternatives are called "prototypes". They represent
  alternatives for which the decision maker has an a priori knowledge
  and that he/she is able to assign to one of the
  ordered classes. If the provided classification of the prototypes cannot be
  described by a Choquet integral, an approximative solution, which minimizes
  the "gap" between the given classification and the one derived from
  the Choquet integral, is proposed. The problem is solved by quadratic
  programming. This function should be used in combination with
  \code{ls.sorting.treatment} which allows to evaluate the model which
  has been built and to assign other alternatives to the ordered classes.
}

\usage{
ls.sorting.capa.ident(n, k, C, cl, d, A.Shapley.preorder = NULL,
A.Shapley.interval = NULL, A.interaction.preorder = NULL,
A.interaction.interval = NULL, A.inter.additive.partition = NULL,
sigf = 5, maxiter = 20, epsilon = 1e-6)
}

\arguments{
  \item{n}{Object of class \code{numeric} containing the
    number of elements of the set on which the object of class
    \code{Mobius.capacity} is to be defined (in short, the number of criteria).}
 
  \item{k}{Object of class \code{numeric} imposing that the solution is at
    most a k-additive capacity (the \enc{Möbius}{Mobius} transform of subsets whose cardinal is
    superior to \code{k} vanishes).}

  \item{C}{Object of class \code{matrix} containing the
    \code{n}-column criteria matrix. Each line of this matrix
    corresponds to a prototype.}

  \item{cl}{Object of class \code{numeric} containing the indexes of the
    classes the alternatives are belonging to (the greater the class
    index, the better the prototype is considered by the decision
    maker). Each class index between min(cl) and max(cl) must be present.}

  \item{d}{Object of class \code{numeric} containing the threshold value
    for the classes, i.e. the minimal "distance" between two neighbor
    classes (e.g. the difference in terms of the Choquet integral of the
    worst prototype of class 3 and the best prototype of
    class 2 should be at least \code{d}).}

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

  \item{sigf}{Precision (default: 5 significant figures). Parameter to
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
    \code{k}-additive solution.}
  
  \item{glob.eval}{The global evaluations satisfying the given classification.}

  \item{how}{Information returned by \code{ipop} (cf. \pkg{kernlab}) on
    the convergence of the solver.}
}

\references{
  K. Fujimoto and T. Murofushi (2000) \emph{Hierarchical decomposition of the
  Choquet integral}, in: Fuzzy Measures and Integrals: Theory and
  Applications, M. Grabisch, T. Murofushi, and M. Sugeno Eds, Physica
  Verlag, pages 95-103.
  
  P. Meyer, M. Roubens (2005), \emph{Choice, Ranking and Sorting in Fuzzy Multiple
  Criteria Decision Aid}, in: J. Figueira, S. Greco, and M. Ehrgott,
  Eds, Multiple Criteria Decision Analysis: State of the Art
  Surveys, volume 78 of International Series in Operations Research and
  Management Science, chapter 12, pages 471-506. Springer Science +
  Business Media, Inc., New York.
}

\seealso{
  \code{\link{Mobius.capacity-class}},
  \cr \code{\link{mini.var.capa.ident}},
  \cr \code{\link{least.squares.capa.ident}},
  \cr \code{\link{ls.sorting.treatment}},
  \cr \code{\link{entropy.capa.ident}}.  
}

\examples{

## generate a random problem with 10 prototypes and 4 criteria
n.proto <- 10 ## prototypes
n <- 4  ## criteria
k <- 4  
d <- 0.1
	
## generating random data for the prototypes
C <- matrix(runif(n.proto*n,0,1),n.proto,n)
cl <- numeric(n.proto)

## the corresponding global evaluations
glob.eval <- numeric(n.proto)
a <- capacity(c(0:(2^n-3),(2^n-3),(2^n-3))/(2^n-3))
for (i in 1:n.proto)
  glob.eval[i] <- Choquet.integral(a,C[i,])

## and the classes for the prototypes
cl[glob.eval <= 0.33] <- 1
cl[glob.eval > 0.33 & glob.eval <= 0.66] <-2
cl[glob.eval > 0.66] <- 3

cl

\dontrun{
# starting the calculations
# search for a capacity which satisfies the constraints
lsc <- ls.sorting.capa.ident(n ,k, C, cl, d)

## output of the quadratic program (ipop, package kernlab)
lsc$how

## the capacity satisfying the constraints
lsc$solution
summary(lsc$solution)
## the global evaluations satisfying the constraints
lsc$glob.eval
}

## let us now add some constraints        

## a Shapley preorder constraint matrix
## Sh(1) > Sh(2)
## Sh(3) > Sh(4)
delta.S <-0.01
Asp <- rbind(c(1,2,delta.S), c(3,4,delta.S))

## a Shapley interval constraint matrix
## 0.1 <= Sh(1) <= 0.2 
Asi <- rbind(c(1,0.1,0.2))
        
## an interaction preorder constraint matrix
## such that I(12) > I(34)
delta.I <- 0.01
Aip <- rbind(c(1,2,3,4,delta.I))
        
## an interaction interval constraint matrix
## i.e. 0.2 <= I(12) <= 0.4 
## and 0 < I(34) <= 1
Aii <- rbind(c(1,2,0.2,0.4), c(3,4,delta.I,1))
        
## an inter-additive partition constraint
## criteria 1,2 and criteria 3,4 are independent 
Aiap <- c(1,1,2,2)


## starting the calculations
## search for a capacity which satisfies the constraints
lsc <- ls.sorting.capa.ident(n ,k, C, cl, d,
                                  A.Shapley.preorder = Asp,
                                  A.Shapley.interval = Asi,
                                  A.interaction.preorder = Aip,
                                  A.interaction.interval = Aii,
                                  A.inter.additive.partition = Aiap)

## output of ipop
lsc$how
## the capacity satisfying the constraints
lsc$solution
summary(lsc$solution)
## the global evaluations satisfying the constraints
lsc$glob.eval

}

\keyword{math}
