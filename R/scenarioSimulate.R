#' @title Simulate an Entire Genome Scenario

#' @author Ahmad Abdel-Azim \email{agabdel.azim@gmail.com}
#' @author Salma Abdel-Azim \email{salma.abdelazim@gmail.com}
#' @author Gamal Abdel-Azim \email{gamal.azim@gmail.com}
#'
#' @name scenarioSimulate
#'
#' @usage scenarioSimulate(Nl, Ng, Psize, runs, nDays, gen.interval, Rm, startingFitness, thr, maxPsize)
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
#' @param runs number of times that the simulation will be looped. Defaults to 10.
#'
#' @description \code{scenarioSimulate} loops the \code{runSimulation} function for a specified amount of runs. The purpose fo such a function is to provide accuracy to SAR through replication.
#'
#'
#'@return As simulation scenarios progress through bacterial generations, the status of the simulation will be printed. The population size, the fitness of each region, and the overall fitness of the population will be printed following each generation. After completion, \code{scenarioSimulate} returns an object of class array and mode numeric. The resulting array will be 3-dimensional with the number of columns equal to 2 plus the number of regions included in the bacterial genome, and the number of rows equal to the number of generations simulated. With 24 hours in a culture day, the number of generations simulated can be calculated as 24 hours multiplied by 60 minutes and divided by \code{gen.interval}. The first column of the returned array includes the population size following each generation. The second column provides the overall fitness of the bacterial population following each generation. Subsequent columns provide the fitness of each region from 1 to n regions. The third dimension consists of sheets for every simulation run, with the number of sheets in the third dimension equal to the number of runs specified in the \code{runs} parameter.
#'
#'
#' @examples x <- scenarioSimulate(Nl = c(1, 6, 10), Ng = c(6, 10, 12), runs = 10)
#'
#'
#' @keywords SAR
#' @keywords resistance
#' @keywords antibiotic
#' @keywords simulation
#' @keywords scenario

#' @export
scenarioSimulate <- function(Psize = 300, runs, nDays = 15, Nl, Ng,
                             gen.interval = 60, Rm = 0.01, startingFitness = 0.51, thr = 0.51, maxPsize = 2000) {

  if(length(Nl) != length(Ng))
    stop("Number of sections in regions and loci are different!")

  allData <- array(0, dim = c(round((24*60)/gen.interval)*nDays, (length(Nl)+2), runs))

  for(i in 1:runs) {
    x <- try(runSimulation(Psize = Psize, Nl = Nl, Ng = Ng, nDays = nDays,
                           gen.interval = gen.interval, Rm = Rm, startingFitness = startingFitness,
                           thr = thr, maxPsize = maxPsize), TRUE)
    if(isTRUE(class(x) == "try-error")) {
      cat("\n", "FAILED RUN", i, "of total", runs, "runs", "\n")
      next
      }
    else {
      allData[,,i] <- x
      cat("Finished run", i, "of total", runs, "runs", "\n")
    }
  }
  return(allData)
}
