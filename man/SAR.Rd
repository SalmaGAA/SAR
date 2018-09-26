\name{runSimulation}
\alias{runSimulation}

\title{Simulation of Ampicillin Resistance in Staphylococcus aureus Bacteria Over Multiple Exposures}

\description{
description to be added later
}

\usage{ runSimulation(gen.interval = 60, Rm = 0.001, Ng = c(6, 10, 12), Nl = c(1, 6, 10),
                      Psize = 300, startingFitness = 0.51, thr = .51, nDays = 15,
                      maxPsize = 2000)
}

\arguments{
  \item{gen.interval}{Generation interval in minutes. Defaults to 60 minutes.}
  \item{Rm}{Rate of mutation in bacterial population. Defaults to 0.001.}
  \item{Ng}{Number of genes per section. One-dimensional array of length equal to number of sections. Defaults to c(6, 10, 12).}
  \item{Nl}{Number of loci per gene per section. One-dimensional array of length equal to number of sections. Defaults to c(1, 6, 10).}
  \item{Psize}{Starting popultaion size. Defaults to 300.}
  \item{startingFitness}{Starting population fitness. Defaults to 0.51.}
  \item{thr}{Bacteria with an overall fitness below this specified threshold do not reproduce to the next generation. Defaults to 0.51.}
  \item{nDays}{Number of days that the bacterial population will be treated with ampicillin. Defaults to 15.}
  \item{maxPsize}{The maximum population size that the bacterial population will reach within a day. Defaults to 2000.}
}

\details{
add details later
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Ahmad Abdel-Azim, Salma Abdel-Azim, and Gamal Abdel-Azim
}
\note{
add later
}

\examples{
x <- runSimulation(gen.interval = 60, Rm = 0.001, Ng = c(6, 10, 12), Nl = c(1, 6, 10),
                      Psize = 100, startingFitness = 0.51, thr = .51, nDays = 2,
                      maxPsize = 2000)
x
}

\keyword{"simulation"}
\keyword{"resistance"}
\keyword{"SAR"}
