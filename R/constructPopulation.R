#' @title Construction of Bacterial Population
#'
#' @author Ahmad Abdel-Azim \email {agabdel.azim@gmail.com}
#' @author Salma Abdel-Azim \email {salma.abdelazim@gmail.com}
#' @author Gamal Abdel-Azim \email {gamal.azim@gmail.com}
#'
#' @name constructPopulation
#'
#' @usage constructPopulation(Psize, Nl, Ng, startingFitness)
#'
#' @param Ng Number of genes per section. One-dimensional array of length equal to number of sections. Defaults to c(6, 10, 12).
#' @param Nl Number of loci per gene per section. One-dimensional array of length equal to number of sections. Defaults to c(1, 6, 10).
#' @param Psize Starting popultaion size. Defaults to 300.
#' @param startingFitness Starting population fitness, represnted by the proportion of mutated bacterial cells in the population. Defaults to 0.51.
#'
#' @description TO BE ADDED LATER
#'
#' @export
# Function to construct initial population
constructPopulation <- function(Psize = 300, Nl = c(1, 6, 10), Ng = c(6, 10, 12), startingFitness = 0.51) {
  # Start with a population that survives the antimicrobial substance with a probability = startingFitness
  # A startingFitness greater than 0 is an indicator of previous spontaneous mutations!
  # starting fitness sampled from a random normal with mean = startingFitness and std dev of 1%

  fitness <- rnorm(mean = startingFitness, sd = 0.01, n = Psize)
  # truncate values to be between 0 and 1 by having negative values = 0 and
  # greater-than-1 values = 1

  fitness[fitness < 0] <- 0; fitness[fitness > 1] <- 1;
  Nc = sum(Ng * Nl)
  P <- matrix(0, nrow = Psize, ncol = Nc)
  P <- cbind(P, fitness)  # hook fitness to the end of each row
  P <- t(apply(P, 1, function(x) setFitness(x, Nl=Nl, Ng=Ng)))

  return(P)
}


setFitness <- function(bCell, Nl, Ng) {

  if(length(Nl) != length(Ng))
    stop("Number of sections in regions and loci are different!")

  fitness <- bCell[length(bCell)]
  bCell <- bCell[-length(bCell)]
  sections <- length(Nl)
  rows <- cumsum(Ng)  # 10  20  30
  rrows <- rep(1:rows[1], each = Nl[1])
  for(i in 2:sections)
    rrows <- c(rrows, rep((rows[i-1]+1):(rows[i]), each = Nl[i]))

  total.genes <- sum(Ng)
  Z <- diag(total.genes)
  Z <- Z[rrows,]

  setGenes = numeric(total.genes)
  setGenes[sample.int(total.genes, size = round(fitness*total.genes))] = 1

  return(setGenes %*% t(Z))

}


#' @examples
#' P <- constructPopulation(Psize = 100, Nl = c(1, 2, 3), Ng = c(10, 10, 10), startingFitness = 0.51)

#'@keywords SAR
#'@keywords resistance
#'@keywords antibiotic
#'@keywords simulation
