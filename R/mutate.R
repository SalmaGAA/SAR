#' @title Mutate a Bacterial Population
#'
#' @author Ahmad Abdel-Azim \email{agabdel.azim@gmail.com}
#' @author Salma Abdel-Azim \email{salma.abdelazim@gmail.com}
#' @author Gamal Abdel-Azim \email{gamal.azim@gmail.com}
#'
#' @name mutate
#'
#' @usage mutate(bCell, Rm)
#'
#' @param Rm Mutation rate of the bacteria.
#' @param bCell A vector of 0's and 1's, used to represent the gnome of a bacterial cell.
#'
#' @description Mutate is a function that is used to simulate the adaptive mutations that occur within a
#' bactrial populaton as a result of the evolutionary stress that an antimicrobial substance imposes.
#'
#' @details Mutate receives a vector, \code{bCell}, to represent the bacterial genome, along with poulation
#' mutation rate, \code{Rm}. Each locus within \code{bCell} is either represented as 1 or 0, depending on the
#' mutation status of that gene, where 0 represnts an unmutated locus, and 1 represnts a mutated locus.
#' The function then changes 0's to 1's with a probability = \code{Rm}. No reverse mutations occur, and there
#' is no effect on mutant genes.
#'
#' @return An object of class "numeric" with the same length as the vector for \code{bCell} that what was inputed.



## Mutate is a function that receives a row in P (bCell) along with poulation mutation rate (Rm)
## The function then changes a 0 to 1 with a probability = Rm. No reverse mutations and no effect on mutant loci.
#mutate0 <- function(bCell, Rm) {
#  n = length(bCell)
#  x = sample(0:1, size = n, prob =  c((1-Rm), Rm), replace = T)
#  bCell = bCell + x
#  bCell[bCell > 1] <- 1
#  return(bCell)
#}


#' @examples To apply the mutate function to a single vector:
#' bg <- c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0)
#' mutate(bCell = bg, Rm = 0.2)
#'
#' To apply mutate to each row of a bacterial population matrix:
#' P2 <- t(apply(P, 1, function(x) mutate(x, Rm = 0.05)))


#'@keywords SAR
#'@keywords resistance
#'@keywords antibiotic
#'@keywords simulation


#' @export
# Mutate only wild type loci
mutate <- function(bCell, Rm) {
  zeros <- which(bCell == 0)
  n = length(zeros)
  x = sample(0:1, size = n, prob =  c((1-Rm), Rm), replace = T)
  bCell[zeros] = bCell[zeros] + x
  return(bCell)
}


