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
#' @description This function constructs an initial bacterial population with designated population parameters, such as starting fitness an genome complexity.
#'
#' @details A binary 2-dimensional matrix is created, where each row represents a distinct bacterial
#' genome in the poulation. The number of rows in the population matrix will be equal to Psize. Within
#' each row, the number of genes will be the cross product of \code{Nl} and \code{Ng}, i.e. \code{crossprod(Nl, Ng)}.
#' For example, if \code{Ng} is set to three sections with 1, 6, and 10 genes, and \code{Nl} is set to 6, 10, and 12
#' loci per gene within each section, the resulting bacterial genome will have 1*6 + 6*10 + 10*12  =
#' 186 loci across all sections. Thus, each row of the resulting matrix will include 186 columns.
#' Each locus in a bacterial genome is assigned either a 1 or 0, depending on the mutation status of
#' that locus, where a zero represnts an unmutated locus, and 1 represnts a mutated locus. \cr\cr
#'
#' The starting fitness of the bacterial population can be set in this function. Starting fitness is
#' defined as the ability of bacteria to resist antimcrobial stress prior to any mutation-inducing
#' treatments. In other words, a  \code{startingFitness} greater than 0 is an indicator of previous spontaneous
#' mutations. The designated value for startingFitness in this function will create a population that
#' survives an antimicrobial substance with a probability equal to \code{startingFitness}. Mutation values are
#' randomly assigned to bacterial cells in the population matrix with mean = \code{startingFitness} and std
#' dev = 1\%. Genes within each bacterial cell are mutated so that the proportion of mutated genes in
#' each bacterial cell is equal to the assigned fitness value.
#'
#' @return A object of class "matrix" with dimensions \code{Psize} by the cross product of \code{Nl} and \code{Ng}
#' is created by this function.
#'

#' @examples
#' P <- constructPopulation(Psize = 100, Nl = c(1, 2, 3), Ng = c(10, 10, 10), startingFitness = 0.51)

#'@keywords SAR
#'@keywords resistance
#'@keywords antibiotic
#'@keywords simulation


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


