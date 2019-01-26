#' @title Simulation of Antibiotic Resistance in Bacteria over Multiple Exposures
#'
#'
#' @author Ahmad Abdel-Azim \email{agabdel.azim@gmail.com}
#' @author Salma Abdel-Azim \email{salma.abdelazim@gmail.com}
#' @author Gamal Abdel-Azim \email{gamal.azim@gmail.com}
#'
#' @name runSimulation
#'
#' @usage runSimulation(Ng, Nl, gen.interval = 60, Rm = 0.001, Psize = 300, startingFitness = 0.51, thr = 0.51, nDays = 15, maxPsize = 2000)
#'
#' @description \code{runSimulation} is used to stochastically simulate the growth of bacteria while under antibiotic stress in an effort to investigate the pattern of resistance in bacteria. The levels of complexity within each region of the bacterial genome can be varied to examine the progression of fitness within a population as well as the trend of the overall growth of bacteria in a culture.
#'
#' @param gen.interval generation interval time, in minutes. Defaults to 60 minutes.
#' @param Rm rate of mutation in bacterial population matrix. Defaults to 0.001.
#' @param Ng one-dimensional array specifying the number of genes within each section of bacterial genome. Length of array should equal number of sections in bacterial genome.
#' @param Nl one-dimensional array specifying the number of mutation sites on each gene along the bacterial genome. Length of array should equal number of sections in bacterial genome.
#' @param Psize starting popultaion size. Defaults to 300.
#' @param startingFitness starting population fitness. Defaults to 0.51.
#' @param thr represents the threshold of fitness which bacteria must meet in order to survive antibiotic stress and proliferate. Bacteria with an overall fitness below this specified threshold do not reproduce to the next generation. Defaults to 0.51.
#' @param nDays number of culture days to be simulated. Defaults to 15.
#' @param maxPsize maximum population size that the bacterial population will be able to reach within a day. Defaults to 2000.
#'
#'
#' @details The \code{runSimulation} function was developed to allow for the investigation of the pattern of fitness development and population growth of bacteria cultured over multiple exposures of antibiotic stress. As mutagenic compounds, antibiotics induce mutations within a bacterial population, which leads to a rapid increase in the fitness of a bacteria. By varing the \code{Nl} and \code{Ng} arguments, the complexity of a bacterial genome can be varied, and the effect such a change will have on a bacterial population can be examined. \cr\cr
#' Bacteria are simulated to be cultured for 24 hours under environmental stress. The argument \code{thr} represents the antibiotic stress placed on a bacterial culture. Bacteria that do not meet this threshold (i.e. bacteria that have not mutated enough in each generation) are barred from proliferation, to represent antibitic stress. In order to avoid population crashes, \code{startingFitness} should be greater than \code{thr}. \cr\cr
#' The argument \code{maxPsize} is set to adhere to the RAM limitations on standard computer machines used to run this simulation.

#'@return \code{runSimulation} returns an object of class matrix with the number of columns equal to 2 plus the number of regions included in the bacterial genome, and the number of rows equal to the number of generations simulated. With 24 hours in a culture day, the number of generations simulated can be calculated as 24 hours multiplied by 60 minutes and divided by \code{gen.interval}. The first column of the returned matrix includes the population size following each generation. The second column provides the overall fitness of the bacterial population following each generation. Subsequent columns provide the fitness of each region from 1 to n regions.

#'@example
#'runSimulation(Ng = c(6, 10, 12), Nl = c(1, 6, 10))
#'

#'
#'@keywords SAR
#'@keywords resistance
#'@keywords antibiotic
#'@keywords simulation
#'


#' @export
runSimulation <- function(gen.interval = 60,
                          Rm = 0.001,
                          Ng,
                          Nl,
                          Psize = 300,
                          startingFitness = 0.51,
                          thr = .51,
                          nDays = 15,
                          maxPsize = 2000) {

  if(length(Nl) != length(Ng))
    stop("Number of sections in regions and loci are different!")

  startPsize = Psize
  # Construct Population --
  P <- constructPopulation(Psize, Nl, Ng, startingFitness)

  sequencedGenes <- rep(1:sum(Ng), rep(Nl, Ng))
  seqSec <- rep(1:length(Ng), Ng)
  gen.number = round((24*60)/gen.interval)
  sections = length(Nl)
  gen.results = matrix(0, nrow=nDays*gen.number, ncol=(2+sections))
  for(i in 1:nDays) {
    gen.number.prev = (i-1)*gen.number
    for(j in 1:gen.number) {
      # Mutate every cell (row) in the population (P)
      P <- t(apply(P, 1, function(x) mutate(x, Rm = Rm)))

      # Get resistance for each cell (row) in the population (P)
      #Rd = apply(P, 1, function(x) getResistance(x, Nl=Nl, Ng=Ng))
      Rd = apply(P, 1, function(x) getResistance.internal(x, sequencedGenes, Ng = Ng))

      # Get resistance of each section
      #sectionRd = t(apply(P, 1, function(x) getSectionResistance(x, Nl=Nl, Ng=Ng)))
      sectionRd = t(apply(P, 1, function(x) getSectionResistance.internal(x, sequencedGenes, seqSec, Ng = Ng)))

      # Reproduce the population P based on a vector of fitness values in Rd
      P = cellDivide(P, Rd, maxPsize, thr=thr)

      # Report generational statistics
      gen.results[gen.number.prev + j,] = reportGeneration(P, Rd, sectionRd)
      cat("Finished gen: ", j, "of day", i, ":", gen.results[gen.number.prev + j,], "\n")
    }
    cat("\n", "Finished day", i, "\n")
    #P <- nextDaySample(P, Rd, startPsize)
    P <- P[sample(1:nrow(P), size=startPsize),]
  }

  return(gen.results)

}


sumSection <- function(gtypeSec, ng, nl) {
  if(nl == 1) sum(gtypeSec)
  else sum(tapply(gtypeSec, rep(1:ng, each=nl), prod))
}


reportGeneration <- function(P, Rd, sectionRd) {
  Psize <- nrow(P)
  avgFit <- mean(Rd)
  avgSectionFit = c(apply(sectionRd, 2, mean))
  return(c(Psize, avgFit, avgSectionFit))
}



nextDaySample <- function(P, Rd, Psize) {
  Rd <- Rd/sum(Rd)
  include <- sample.int(length(Rd), size=Psize, prob=Rd)
  return(P[include,])
}
