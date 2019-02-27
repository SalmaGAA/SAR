#' @title Division of Bacterial Population
#'
#' @author Ahmad Abdel-Azim \email{agabdel.azim@gmail.com}
#' @author Salma Abdel-Azim \email{salma.abdelazim@gmail.com}
#' @author Gamal Abdel-Azim \email{gamal.azim@gmail.com}
#'
#' @name cellDivide
#'
#' @description Functions to manipulate individual operations of division within a cell population.
#'
#' @usage cellDivide(P, Rd, maxPsize, thr)
#'
#' @param thr Bacteria with an overall fitness below this specified threshold do not reproduce to the next generation. Defaults to 0.51.
#' @param maxPsize The maximum population size that the bacterial population will reach within a day. Defaults to 2000.
#' @param P The population of bacteria arranged in a two-dimensional matrix of 0 and 1 values, where 0 represents naive loci and 1 represents mutated loci.
#' @param Rd ADD LATER
#'
#'
#' @details 


#' @examples
#' P2 <- cellDivide(P, Rd = 0.2, maxPsize = 500, thr = 0.51)
#'
#' @keywords SAR
#' @keywords resistance
#' @keywords antibiotic
#' @keywords simulation



#' @export
cellDivide <- function(P, Rd, maxPsize, thr) {

  if(length(Rd) != nrow(P))
    stop("Number of elements in Rd not equal to number of rows in P")

  included <- (Rd >= thr)
  Rd <- Rd[included]
  P <- P[included,]
  n = length(Rd)
  x <- runif(min = 0, max = 1, n=n)
  y <- x <= Rd
  P <- P[y,]; Rd <- Rd[y]
  P <- rbind(P, P); Rd <- c(Rd, Rd)

  # shuffle if nrows <= maxPsize, otherwise pick 'the best' maxPsize!
  if(nrow(P) <= maxPsize)
    P <- P[sample(1:nrow(P), size = nrow(P)),]
  else P <- P[sample(1:nrow(P), size = maxPsize, prob=(Rd/sum(Rd))),]
  P
}


#' # A function to reproduce the population by cell fission/division
#' cellDivide0 <- function(P, Rd) {
#'
#'   if(length(Rd) != nrow(P))
#'     stop("Number of elements in Rd != number of rows in P")
#'
#'   n = length(Rd)
#'   x <- runif(min = 0, max = 1, n=n)
#'   y <- x <= Rd
#'   P <- P[y,]
#'   P <- rbind(P, P)
#'   P <- P[sample(1:nrow(P), size = nrow(P)),]
#'   P
#' }
#'
#'
#' cellDivide1 <- function(P, Rd, sel.intensity = .5) {
#'
#'   if(length(Rd) != nrow(P))
#'     stop("Number of elements in Rd != number of rows in P")
#'
#'   # Use order()
#'   x <- order(Rd)
#'   at = round((1-sel.intensity)*length(Rd)) + 1
#'   y <- Rd >= Rd[x[at]]
#'
#'   P <- P[y,]
#'   P <- rbind(P, P)
#'   P <- P[sample(1:nrow(P), size = nrow(P)),]
#'   P
#' }
