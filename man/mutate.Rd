% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mutate.R
\name{mutate}
\alias{mutate}
\title{Mutate a Bacterial Population}
\usage{
mutate(bCell, Rm)
}
\arguments{
\item{bCell}{A vector of 0's and 1's, used to represent the gnome of a bacterial cell.}

\item{Rm}{Mutation rate of the bacteria.}
}
\value{
An object of class "numeric" with the same length as the vector for \code{bCell} that what was inputed.
}
\description{
Mutate is a function that is used to simulate the adaptive mutations that occur within a
bactrial populaton as a result of the evolutionary stress that an antimicrobial substance imposes.
}
\details{
Mutate receives a vector, \code{bCell}, to represent the bacterial genome, along with poulation
mutation rate, \code{Rm}. Each locus within \code{bCell} is either represented as 1 or 0, depending on the
mutation status of that gene, where 0 represnts an unmutated locus, and 1 represnts a mutated locus.
The function then changes 0's to 1's with a probability = \code{Rm}. No reverse mutations occur, and there
is no effect on mutant genes.
}
\examples{
To apply the mutate function to a single vector:
bg <- c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)
mutate(bCell = bg, Rm = 0.2)

To apply mutate to each row of a bacterial population matrix:
P2 <- t(apply(P, 1, function(x) mutate(x, Rm = 0.05)))
}
\author{
Ahmad Abdel-Azim \email{agabdel.azim@gmail.com}

Salma Abdel-Azim \email{salma.abdelazim@gmail.com}

Gamal Abdel-Azim \email{gamal.azim@gmail.com}
}
\keyword{SAR}
\keyword{antibiotic}
\keyword{resistance}
\keyword{simulation}
