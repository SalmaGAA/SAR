#' @title Set Fitness for a Bacterial Population
#'
#' @author Ahmad Abdel-Azim \email {agabdel.azim@gmail.com}
#' @author Salma Abdel-Azim \email {salma.abdelazim@gmail.com}
#' @author Gamal Abdel-Azim \email {gamal.azim@gmail.com}
#'
#' @name setFitness
#'
#' @usage setFitness(bCell, Nl, Ng)
#'
#' @param Ng Number of genes per section. One-dimensional array of length equal to number of sections. Defaults to c(6, 10, 12).
#' @param Nl Number of loci per gene per section. One-dimensional array of length equal to number of sections. Defaults to c(1, 6, 10).
#' @param bCell One bacterial cell in the population, represented as a row in the population matrix, P.
#'
#' @description TO BE ADDED LATER
#'
#' @export
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
#' P2 <- t(apply(P, 1, function(x) setFitness(x, Nl = c(1, 2, 3), Ng = c(10, 10, 10))))

#'@keywords SAR
#'@keywords resistance
#'@keywords antibiotic
#'@keywords simulation
