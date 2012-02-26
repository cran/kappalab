\name{ls.sorting.treatment}
\alias{ls.sorting.treatment}

\title{Least squares capacity identification in the framework of a
  sorting procedure: evaluation of the determined capacity}

\description{
  This function assigns alternatives to classes and optionally compares
  the obtained classification to a given one. The classes are
  described by a set of prototypes (well-known alternatives for the
  decision maker) and a   capacity (which can, for instance, be determined by
  \code{ls.sorting.capa.ident}).   This function (in combination with
  \code{ls.sorting.capa.ident}) is an
  implementation of the TOMASO method; see Meyer and Roubens (2005).
}

\usage{
ls.sorting.treatment(P, cl.proto, a, A, cl.orig.A = NULL)
}

\arguments{
  \item{P}{Object of class \code{matrix} containing the
    \code{n}-column criteria matrix. Each line of this matrix
    corresponds to a prototype. A prototype is an alternative for which
    the class is known beforehand.}

  \item{cl.proto}{Object of class \code{numeric} containing the indexes of the
    classes the prototypes \code{P} are belonging to (the greater the class
    index, the better the prototype is considered by the decision maker).}

  \item{a}{Object of class \code{Mobius.capacity} used to model the
    classes determined by \code{P} by a Choquet integral.}

  \item{A}{Object of class \code{matrix} containing the \code{n}-column
    criteria matrix representing the alternatives which have to be classified.} 

  \item{cl.orig.A}{Object of class \code{numeric}
    containing the "true" classes of the alternatives of \code{A}. This
    can be used for the evaluation of the quality of the model.}
}
 
\value{
  The function returns a list structured as follows:
	
  \item{correct.A}{Object of class \code{matrix} which contains the
    different types of assignments of the elements of \code{A}. In
    columns the alternatives. First line: correct assignment to a single class
    (assignment of degree 0). Second line: correct ambiguous assignment to two
    classes (assignment of degree 1), etc. Last line: bad assignment. In
  case no \code{orig.class.A} is given, \code{correct.A} is \code{NULL}.}

  \item{class.A}{Object of class \code{matrix} which contains the
    assignments of the elements of \code{A}. In columns, the
    alternatives. First line, lower class. Second line, upper class. }

  \item{eval.correct}{Object of class \code{numeric} which contains the
    ratio assignments over number of elements of \code{A}(for each type
    of assignment; the first element is the ratio of correct
    assignments). In
  case no \code{orig.class.A} is given, \code{eval.correct} is \code{NULL}.}

  \item{minmax.P}{Object of class \code{matrix} which contains the min
    and the max of each class, according to the prototypes. In columns,
    the classes, first line: the minimum, second line: the maximum.}

  \item{Choquet.A}{Object of class \code{numeric} which contains the
    Choquet integral of the evaluations of the alternatives of \code{A}.}
}

\references{
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
  \cr \code{\link{mini.dist.capa.ident}},
  \cr \code{\link{ls.sorting.capa.ident}},
  \cr \code{\link{least.squares.capa.ident}},
  \cr \code{\link{heuristic.ls.capa.ident}},
  \cr \code{\link{entropy.capa.ident}}.  
}

\examples{

## generate a random problem with 10 prototypes and 4 criteria
n.proto <- 10 ## prototypes
n <- 4  ## criteria
P <- matrix(runif(n.proto*n,0,1),n.proto,n)

## the corresponding global scores, based on a randomly generated
## capacity a
glob.eval <- numeric(n.proto)
a <- capacity(c(0:(2^n-3),(2^n-3),(2^n-3))/(2^n-3))
for (i in 1:n.proto)
  glob.eval[i] <- Choquet.integral(a,P[i,])

## based on these global scores, let us create a classification (3 classes)
cl.proto<-numeric(n.proto)
cl.proto[glob.eval <= 0.33] <- 1
cl.proto[glob.eval > 0.33 & glob.eval<=0.66] <-2
cl.proto[glob.eval > 0.66] <- 3

## search for a capacity which satisfies the constraints
lsc <- ls.sorting.capa.ident(n ,4, P, cl.proto, 0.1)

## output of the QP
lsc$how

## analyse the quality of the model (classify the prototypes by the
## model and compare both assignments)
lst <- ls.sorting.treatment(P,cl.proto,lsc$solution,P,cl.proto)

## assignments of the prototypes
lst$class.A
## assignment types
lst$correct.A
## evaluation
lst$eval.correct

## generate a second set of random alternatives (A)
## their "correct" class is determined as beforehand with the
## randomly generated capacity a
## the goal is to see if we can reproduce this classification
## by the capacity learnt from the prototypes

## a randomly generated criteria matrix of 10 alternatives
A <- matrix(runif(10*n,0,1),10,n)
cl.orig.A <-numeric(10)
	
## the corresponding global scores
glob.eval.A <- numeric(10)
for (i in 1:10)
  glob.eval.A[i] <- Choquet.integral(a,A[i,])

## based on these global scores, let us determine a classification
cl.orig.A[glob.eval.A <= 0.33] <- 1
cl.orig.A[glob.eval.A>0.33 & glob.eval.A<=0.66] <-2
cl.orig.A[glob.eval.A > 0.66] <- 3

## let us now classify the alternatives of A according to the model
## built on P

lst <- ls.sorting.treatment(P,cl.proto,lsc$solution,A,cl.orig.A)

## assignment of the alternatives of A
lst$class.A
## type of assignments
lst$correct.A
## evaluation
lst$eval.correct

## show the learnt capacity
## x11()
## barplot(Shapley.value(lsc$solution), main="Learnt capacity", sub="Shapley")
## summary of the learnt capacity
lsc$solution
summary(lsc$solution)

}

\keyword{math}