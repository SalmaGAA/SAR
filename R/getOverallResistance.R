#' @title Measure Overall Resistance of Bacterial Population
#'
#' @author Ahmad Abdel-Azim \email{agabdel.azim@gmail.com}
#' @author Salma Abdel-Azim \email{salma.abdelazim@gmail.com}
#' @author Gamal Abdel-Azim \email{gamal.azim@gmail.com}
#'
#' @name getOverallResistance
#'
#' @usage getResistance(gtype, Nl, Ng)
#'
#' @param Ng Number of genes per section. One-dimensional array of length equal to number of sections. Defaults to c(6, 10, 12).
#' @param Nl Number of loci per gene per section. One-dimensional array of length equal to number of sections. Defaults to c(1, 6, 10).
#' @param gtype ADD LATER
#'
#'@description TO BE ADDED LATER


# # Quanify Resistance
# getResistance0 <- function(gtype, Nl, Ng) {
#
#   if(length(Nl) != length(Ng))
#     stop("Number of sections in regions and loci are different!")
#
#   sections <- length(Nl)
#   rows <- cumsum(Ng)  # 10  20  30
#   rrows <- rep(1:rows[1], each = Nl[1])
#   for(i in 2:sections)
#     rrows <- c(rrows, rep((rows[i-1]+1):(rows[i]), each = Nl[i]))
#
#   total.genes <- sum(Ng)
#   Z <- diag(total.genes)
#   Z <- Z[rrows,]
#
#   X <- gtype %*% Z
#   A <- numeric(total.genes)
#   A[1:rows[1]] <- Nl[1] - 1
#   for(i in 2:sections) A[(rows[i - 1] + 1):rows[i]] <- Nl[i] - 1
#   resist = c(X) - A
#   resist <- length(resist[resist == 1])
#
#   # Scale resistance before returning the number of 1's
#   # XXX resist = startingFitness + (resist/total.genes)*(1 - startingFitness)
#   resist = (resist/total.genes)
#   return(resist)
# }
#

#' @export
# Quantify Resistance
getResistance <- function(gtype, Nl, Ng) {

  if(length(Nl) != length(Ng))
    stop("Number of sections in regions and loci are different!")

  sections <- length(Nl)
  Ngl = Ng * Nl
  starts <- c(1, cumsum(Ngl)+1)  # cumulative sum for loci
  mutated.genes = 0
  for(i in 1:sections) {
    mutated.genes = mutated.genes + sumSection(gtype[starts[i]:(starts[i+1]-1)], Ng[i], Nl[i])
    }
  # Scale resistance before returning the number of 1's
  return(mutated.genes/sum(Ng))
}

# Internal version for speed
getResistance.internal <- function(gtype, sequencedGenes) {
  return(sum(tapply(gtype, sequencedGenes, prod))/sum(Ng))
}


#' @examples
#' resist <- apply(P, 1, function(x) getResistance(x, Nl = c(1, 2, 3), Ng= c(10, 10, 10)))
#' resist
#'
#'
#' @keywords SAR
#' @keywords resistance
#' @keywords antibiotic
#' @keywords simulation

