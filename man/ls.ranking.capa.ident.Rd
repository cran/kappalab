\name{ls.ranking.capa.ident}
\alias{ls.ranking.capa.ident}

\title{Least squares capacity identification in the framework of a ranking procedure}

\description{
	Ranking alternatives means ordering them from the best to the worst alternative. The aim of the implemented method is to model a given ranking by means of a Choquet integral. The result of the function is an object of class \code{Mobius.capacity}. This function is an implementation of the TOMASO method (see Meyer and Roubens (2005)) in the particular ranking framework. The input data are given under the form of a set of alternatives and a partial weak order, each alternative being described according to a set of criteria. These well-known alternatives are called "prototypes". They represent alternatives for which the decision maker has an a priori knowledge and for which he/she is able to build a ranking. If the provided ranking (partial weak order) of the prototypes cannot be described by a Choquet integral, an approximative solution, which minimizes the "gap" between the given ranking and the one derived from the Choquet integral, is proposed. The problem is solved by quadratic programming. 
}

\usage{
ls.ranking.capa.ident(n, k, C, rk, d, A.Shapley.preorder = NULL,
A.Shapley.interval = NULL, A.interaction.preorder = NULL,
A.interaction.interval = NULL, A.inter.additive.partition = NULL,
sigf = 5, maxiter = 20, epsilon = 1e-6)
}

\arguments{
  \item{n}{Object of class \code{numeric} containing the
    number of elements of the set on which the object of class
    \code{Mobius.capacity} is to be defined (in short, the number of criteria).}
 
  \item{k}{Object of class \code{numeric} imposing that the solution is at
    most a k-additive capacity (the Moebius transform of subsets whose cardinal is
    superior to \code{k} vanishes).}

  \item{C}{Object of class \code{matrix} containing the
    \code{n}-column criteria matrix. Each line of this matrix
    corresponds to a prototype.}

  \item{rk}{Object of class \code{matrix} containing the
    constraints relative to the preorder of the prototypes. Each line
    of this 2-column matrix corresponds to one constraint of the type
    "the alternative \code{i} is prefered to the alternative \code{j}". A line is structured as follows: the first element
    encodes \code{i}, the second \code{j}.}

  \item{d}{Object of class \code{numeric} containing the threshold value, i.e. the minimal "distance" between two neighbor alternatives in the given ranking (e.g. the difference in terms of the Choquet integral of the
    a prototype with rank 3 and a prototype with rank 4 should be at least \code{d}).}

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
    \code{Mobius.capacity} containing the Moebius transform of the
    \code{k}-additive solution.}
  
  \item{glob.eval}{The global evaluations satisfying the given ranking.}

  \item{how}{Information returned by \code{ipop} (cf. \pkg{kernlab}) on
    the convergence of the solver.}

	\item{rk.C}{The ranks of the prototypes}
	\item{Choquet.C}{The Choquet integral of the prototypes}

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
  \cr \code{\link{lin.prog.capa.ident}},
  \cr \code{\link{mini.var.capa.ident}},
  \cr \code{\link{mini.dist.capa.ident}},
  \cr \code{\link{least.squares.capa.ident}},
  \cr \code{\link{heuristic.ls.capa.ident}},
  \cr \code{\link{entropy.capa.ident}}.  
}

\examples{

arthur <- c(1, 1, 0.75, 0.25)
lancelot <- c(0.75, 0.75, 0.75, 0.75)
yvain <- c(1, 0.625, 0.5, 1)
perceval <- c(0.25, 0.5, 0.75, 0.75)
erec <- c(0.375, 1, 0.5 , 0.75)

P <- rbind(arthur, lancelot, yvain, perceval, erec)

# lancelot > erec; yvain > erec, erec > perceval, erec > arthur
rk.proto <- rbind(c(2,5), c(3,5), c(5,4), c(5,1))

n<-4
k<-2
d<-0.1

## search for a capacity which satisfies the constraints
lrc <- ls.ranking.capa.ident(n ,k, P, rk.proto, d)

lrc

}

\keyword{math}
