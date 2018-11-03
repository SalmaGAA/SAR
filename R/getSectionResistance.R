#' @title Measure Resistance of Each Section of a Bacterial Genome
#'
#' @author Ahmad Abdel-Azim \email {agabdel.azim@gmail.com}
#' @author Salma Abdel-Azim \email {salma.abdelazim@gmail.com}
#' @author Gamal Abdel-Azim \email {gamal.azim@gmail.com}
#'
#' @name getSectionResistance
#'
#' @usage getSectionResistance(gtype, Nl, Ng)
#'
#' @param Ng Number of genes per section. One-dimensional array of length equal to number of sections. Defaults to c(6, 10, 12).
#' @param Nl Number of loci per gene per section. One-dimensional array of length equal to number of sections. Defaults to c(1, 6, 10).
#' @param gtype ADD LATER
#'
#'@description TO BE ADDED LATER


# # Receive a vector of locus genotypes to calculate section ressiatnce values
# # Receive Nl = loci per gene per section and Ng = genes per section
# # Return resistance gained for each section of 'gtype' = fraction of fully mutated genes per section
# getSectionResistance0 <- function(gtype, Nl, Ng) {
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
#   X <- gtype %*% Z
#   A <- numeric(total.genes)
#
#   A[1:rows[1]] <- Nl[1]-1
#   for(i in 2:sections)
#     A[(rows[i-1]+1):rows[i]] <- Nl[i]-1
#
#   resist <- c(X)-A
#   sectionResist <- numeric(sections)
#
#   part = resist[1:rows[1]]
#   sectionResist[1] <- length(part[part == 1])
#
#   for(i in 2:sections) {
#     part = resist[(rows[i-1]+1) : rows[i]]
#     sectionResist[i] <- length(part[part == 1])
#   }
#
#   # return fraction of fully mutated genes per section
#   return(sectionResist/Ng)
# }
#

#' @export
# More efficient version
getSectionResistance <- function(gtype, Nl, Ng) {

  if(length(Nl) != length(Ng))
    stop("Number of sections in regions and loci are different!")

  sections <- length(Nl)
  Ngl = Ng * Nl
  starts <- c(1, cumsum(Ngl)+1)  # cumulative sum for loci
  sectionResist <- numeric(sections)
  for(i in 1:sections) {
    sectionResist[i] = sumSection(gtype[starts[i]:(starts[i+1]-1)], Ng[i], Nl[i])
  }

  # return fraction of fully mutated genes per section
  return(sectionResist/Ng)
}

#' @examples
#' sectionResist <- t(apply(P, 1, function(x) getSectionResistance(x, Nl = c(1, 2, 3), Ng = c(10, 10, 10))))
#' sectionResist
#'
#' @keywords SAR
#' @keywords resistance
#' @keywords antibiotic
#' @keywords simulation

