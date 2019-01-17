#' @title Simulation of Antibiotic Resistance in Bacteria over Multiple Exposures
#'
#'
#' @author Ahmad Abdel-Azim \email{agabdel.azim@gmail.com}
#' @author Salma Abdel-Azim \email{salma.abdelazim@gmail.com}
#' @author Gamal Abdel-Azim \email{gamal.azim@gmail.com}
#'
#' @name runSimulation
#'
#' @usage runSimulation(gen.interval, Rm, Ng, Nl, Psize, startingFitness, thr, nDays, maxPsize)
#'
#' @description TO BE ADDED LATER
#'
#' @param gen.interval Generation interval in minutes. Defaults to 60 minutes.
#' @param Rm Rate of mutation in bacterial population. Defaults to 0.001.
#' @param Ng Number of genes per section. One-dimensional array of length equal to number of sections. Defaults to c(6, 10, 12).
#' @param Nl Number of loci per gene per section. One-dimensional array of length equal to number of sections. Defaults to c(1, 6, 10).
#' @param Psize Starting popultaion size. Defaults to 300.
#' @param startingFitness Starting population fitness. Defaults to 0.51.
#' @param thr Bacteria with an overall fitness below this specified threshold do not reproduce to the next generation. Defaults to 0.51.
#' @param nDays Number of days that the bacterial population will be treated with ampicillin. Defaults to 15.
#' @param maxPsize The maximum population size that the bacterial population will reach within a day. Defaults to 2000.
#'
#'


#'@example
#'x <- runSimulation(Ng = c(6, 10, 12), Nl = c(1, 6, 10))
#'x

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
