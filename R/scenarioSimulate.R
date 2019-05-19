#' @title Simulation of Antibiotic Resistance in Bacteria over Multiple Exposures

#' @author Ahmad Abdel-Azim \email{agabdel.azim@gmail.com}
#' @author Salma Abdel-Azim \email{salma.abdelazim@gmail.com}
#' @author Gamal Abdel-Azim \email{gamal.azim@gmail.com}
#'
#' @name scenarioSimulate
#'
#' @usage scenarioSimulate(Nl, Ng, Psize = 300, runs = 10, nDays = 15, gen.interval = 60, Rm = 0.001, startingFitness = 0.51, thr = 0.51, maxPsize = 2000)
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
#' @param runs number of times that the simulation will be repeated. Defaults to 10.
#'
#' @description \code{scenarioSimulate} is used to stochastically simulate the growth of bacteria while under antibiotic stress in an effort to investigate the pattern of resistance in bacteria. The levels of complexity within each region of the bacterial genome can be varied to examine the progression of fitness within a population as well as the trend of the overall growth of bacteria in a culture.
#'
#'
#'@details The \code{scenarioSimulate} function was developed to allow for the investigation of the pattern of fitness development and population growth of bacteria cultured over multiple exposures of antibiotic stress. As mutagenic compounds, antibiotics induce mutations within a bacterial population, which leads to a rapid increase in the fitness of a bacteria. By varing the \code{Nl} and \code{Ng} arguments, the complexity of a bacterial genome can be varied, and the effect such a change will have on a bacterial population can be examined. \cr\cr
#' Bacteria are simulated to be cultured for 24 hours under environmental stress. The argument \code{thr} represents the antibiotic stress placed on a bacterial culture. Bacteria that do not meet this threshold (i.e. bacteria that have not mutated enough in each generation) are barred from proliferation, to represent antibitic stress. In order to avoid population crashes, \code{startingFitness} should be greater than \code{thr}. \cr\cr
#' The argument \code{maxPsize} represents the carrying capacity of the bacterial population. It is also set to adhere to the RAM limitations on standard computer machines used to run this simulation. The argument \code{runs} can be set to indicate the number of times the simulation will be run; more simulation runs provide greater accuracy in the depiction of the simulated population behavior.
#'
#'
#'@return As simulation scenarios progress through bacterial generations, the status of the simulation will be printed. The population size, the fitness of each region, and the overall fitness of the population will be printed following each generation. After completion, \code{scenarioSimulate} returns an object of class array and mode numeric. The resulting array will be 3-dimensional with the number of columns equal to 2 plus the number of regions included in the bacterial genome, and the number of rows equal to the number of generations simulated. With 24 hours in a culture day, the number of generations simulated can be calculated as 24 hours multiplied by 60 minutes and divided by \code{gen.interval}. The first column of the returned array includes the population size following each generation. The second column provides the overall fitness of the bacterial population following each generation. Subsequent columns provide the fitness of each region from 1 to n regions. If the argument \code{runs} is set equal to a value greater than 1, the third dimension will consist of sheets for every simulation run, with the number of sheets equal to the number of runs specified in the \code{runs} parameter.
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
