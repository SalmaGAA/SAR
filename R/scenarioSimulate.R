#' @title Simulate an Entire Genome Scenario

#' @author Ahmad Abdel-Azim \email {agabdel.azim@gmail.com}
#' @author Salma Abdel-Azim \email {salma.abdelazim@gmail.com}
#' @author Gamal Abdel-Azim \email {gamal.azim@gmail.com}
#'
#' @name scenarioSimulate
#'
#' @usage scenarioSimulate(Psize, Nl, Ng, runs, nDays, gen.interval, Rm, startingFitness, thr, maxPsize)
#'
#' @param Psize Starting popultaion size. Defaults to 300.
#' @param Nl Number of loci per gene per section. One-dimensional array of length equal to number of sections. Defaults to c(1, 6, 10).
#' @param Ng Number of genes per section. One-dimensional array of length equal to number of sections. Defaults to c(6, 10, 12).
#' @param runs ADD LATER
#' @param nDays ADD LATER
#' @param gen.interval ADD LATER
#' @param Rm ADD LATER
#' @param startingFitness Starting population fitness, represnted by the proportion of mutated bacterial cells in the population. Defaults to 0.51.
#' @param thr ADD LATER
#' @param maxPsize ADD LATER
#'
#' @description ADD LATER
#'
#'
#' @details ADD LATER
#'
#'
#' @return ADD LATER
#'
#'
#' @examples x <- scenarioSimulate()
#'
#'
#' @keywords SAR
#' @keywords resistance
#' @keywords antibiotic
#' @keywords simulation

#' @export
scenarioSimulate <- function(Psize = 300, Nl = c(1, 6, 10), Ng = c(6, 10, 12), runs = 10, nDays = 20,
                             gen.interval = 60, Rm = 0.01, startingFitness = 0.51, thr = .51, maxPsize = 2000) {

  if(length(Nl) != length(Ng))
    stop("Number of sections in regions and loci are different!")

  allData <- array(0, dim = c(round((24*60)/gen.interval)*nDays, (length(Nl)+2), runs))

  for(i in 1:runs) {
    x <- try(runSimulation(Psize = Psize, Nl = Nl, Ng = Ng, nDays = nDays,
                           gen.interval = gen.interval, Rm = Rm, startingFitness = startingFitness,
                           thr = thr, maxPsize = maxPsize), TRUE)
    if(isTRUE(class(x) == "try-error")) {next}
    else {allData[,,i] <- x}
  }
  return(allData)
}
