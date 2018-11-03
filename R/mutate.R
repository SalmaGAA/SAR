#' @title Mutate a Bacterial Population
#'
#' @author Ahmad Abdel-Azim \email {agabdel.azim@gmail.com}
#' @author Salma Abdel-Azim \email {salma.abdelazim@gmail.com}
#' @author Gamal Abdel-Azim \email {gamal.azim@gmail.com}
#'
#' @name mutate
#'
#' @usage mutate(bCell, Rm)
#'
#' @param Rm Rate of mutation in bacterial population.
#' @param bCell One bacterial cell in the population, represented as a row in the population matrix, P.
#'
#' @description TO BE ADDED LATER


## Mutate is a function that receives a row in P (bCell) along with poulation mutation rate (Rm)
## The function then changes a 0 to 1 with a probability = Rm. No reverse mutations and no effect on mutant loci.
#mutate0 <- function(bCell, Rm) {
#  n = length(bCell)
#  x = sample(0:1, size = n, prob =  c((1-Rm), Rm), replace = T)
#  bCell = bCell + x
#  bCell[bCell > 1] <- 1
#  return(bCell)
#}


#' @export
# Mutate only wild type loci
mutate <- function(bCell, Rm) {
  zeros <- which(bCell == 0)
  n = length(zeros)
  x = sample(0:1, size = n, prob =  c((1-Rm), Rm), replace = T)
  bCell[zeros] = bCell[zeros] + x
  return(bCell)
}


#' @example
#' P <- t(apply(P, 1, function(x) mutate(x, Rm = 0.05)))
#'
#'
#'@keywords SAR
#'@keywords resistance
#'@keywords antibiotic
#'@keywords simulation


