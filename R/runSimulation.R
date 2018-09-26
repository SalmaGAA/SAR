#9/25/2018
# Run Simulation of Ampicillin Resistance in Staphylococcus aureus


# gen.interval = Generation interval in minutes. Defaults to 60.
# Rm = Rate of mutation in bacterial population. Defaults to 0.001.
# Ng = Number of genes per section. One-dimensional array of length equal to number of sections. Defaults to c(6, 10, 12)
# Nl = Number of loci per gene per section. One-dimensional array of length equal to number of sections. Defaults to c(1, 6, 10)
# Psize = Starting popultaion size. Defaults to 300.
# startingFitness = Starting population fitness. Defaults to 0.51.
# thr = Bacteria with an overall fitness below the threshold do not reproduce to the next generation. Defaults to 0.51.
# nDays = Number of days that the bacterial population will be treated with ampicillin. Defaults to 15.
# maxPsize = The maximum population size that the bacterial population will reach within a day. Defaults to 2000.


runSimulation <- function(gen.interval = 60,
                          Rm = 0.001,
                          Ng = c(6, 10, 12),
                          Nl = c(1, 6, 10),
                          Psize = 300,
                          startingFitness = 0.51,
                          thr = .51,
                          nDays = 15,
                          maxPsize = 2000) {

  if(length(Nl) != length(Ng))
    stop("Number of sections in regions and loci are different!")

  startPsize = Psize
  # Construct Population --
  P <- constructPopulation(Psize, Nl, Ng, startingFitness)

  gen.number = round((24*60)/gen.interval)
  sections = length(Nl)
  gen.results = matrix(0, nrow=nDays*gen.number, ncol=(2+sections))
  for(i in 1:nDays) {
    gen.number.prev = (i-1)*gen.number
    for(j in 1:gen.number) {
      # Mutate every cell (row) in the population (P)
      P <- t(apply(P, 1, function(x) mutate(x, Rm = Rm)))

      # Get resistance for each cell (row) in the population (P)
      Rd = apply(P, 1, function(x) getResistance(x, Nl=Nl, Ng=Ng))

      # Get resistance of each section
      sectionRd = t(apply(P, 1, function(x) getSectionResistance(x, Nl=Nl, Ng=Ng)))

      # Reproduce the population P based on a vector of ftness values in Rd
      P = cellDivide(P, Rd, maxPsize, thr=thr)

      # Report generational statistics
      gen.results[gen.number.prev + j,] = reportGeneration(P, Rd, sectionRd)
      cat("Finished gen: ", j, "of day", i, ":", gen.results[gen.number.prev + j,], "\n")
    }
    cat("\n", "Finished day", i, "\n")
    #P <- nextDaySample(P, Rd, startPsize)
    P <- P[sample(1:nrow(P), size=startPsize),]
  }

  return(gen.results)

}

#
# Function to construct initial population
# ========================================
constructPopulation <- function(Psize, Nl, Ng, startingFitness) {
  # Start with a population that survives the antimicrobial substance with a probability = startingFitness
  #    A greater-than-0 startingFirness is an indicator of previous spontaneous mutations!

  # starting fitness sampled from a random normal with mean = startingFitness and std dev of 1%
  fitness <- rnorm(mean = startingFitness, sd = 0.01, n = Psize)

  # truncate values to be between 0 and 1 by having negative values = 0 and
  # greater-than-1 values = 1
  fitness[fitness < 0] <- 0; fitness[fitness > 1] <- 1;

  Nc = sum(Ng * Nl)
  P <- matrix(0, nrow = Psize, ncol = Nc)
  P <- cbind(P, fitness)  # hock fitness to the end of each row
  P <- t(apply(P, 1, function(x) setFitness(x, Nl=Nl, Ng=Ng)))

  return(P)
}

#
# Mutate is a function that receives a row in P (bCell) along with poulation mutation rate (Rm)
# The function then changes a 0 to 1 with a probability = Rm. No reverse mutations and no effect on mutant loci.
# ==============================================================================================================
mutate0 <- function(bCell, Rm) {
  n = length(bCell)
  x = sample(0:1, size = n, prob =  c((1-Rm), Rm), replace = T)
  bCell = bCell + x
  bCell[bCell > 1] <- 1
  return(bCell)
}

# Mutate only wild type loci
mutate <- function(bCell, Rm) {
  zeros <- which(bCell == 0)
  n = length(zeros)
  x = sample(0:1, size = n, prob =  c((1-Rm), Rm), replace = T)
  bCell[zeros] = bCell[zeros] + x
  return(bCell)
}

#
# A function to reproduce the population by cell fission/division
# Receives

cellDivide0 <- function(P, Rd) {

  if(length(Rd) != nrow(P))
    stop("Number of elements in Rd != number of rows in P")

  n = length(Rd)
  x <- runif(min = 0, max = 1, n=n)
  y <- x <= Rd
  P <- P[y,]
  P <- rbind(P, P)
  P <- P[sample(1:nrow(P), size = nrow(P)),]
  P
}

cellDivide1 <- function(P, Rd, sel.intensity = .5) {

  if(length(Rd) != nrow(P))
    stop("Number of elements in Rd != number of rows in P")

  # Use order()
  x <- order(Rd)
  at = round((1-sel.intensity)*length(Rd)) + 1
  y <- Rd >= Rd[x[at]]

  P <- P[y,]
  P <- rbind(P, P)
  P <- P[sample(1:nrow(P), size = nrow(P)),]
  P
}

cellDivide <- function(P, Rd, maxPsize, thr) {

  if(length(Rd) != nrow(P))
    stop("Number of elements in Rd != number of rows in P")

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


# Quanify Resistance
getResistance0 <- function(gtype, Nl, Ng) {

  if(length(Nl) != length(Ng))
    stop("Number of sections in regions and loci are different!")

  sections <- length(Nl)
  rows <- cumsum(Ng)  # 10  20  30
  rrows <- rep(1:rows[1], each = Nl[1])
  for(i in 2:sections)
    rrows <- c(rrows, rep((rows[i-1]+1):(rows[i]), each = Nl[i]))

  total.genes <- sum(Ng)
  Z <- diag(total.genes)
  Z <- Z[rrows,]

  X <- gtype %*% Z
  A <- numeric(total.genes)
  A[1:rows[1]] <- Nl[1] - 1
  for(i in 2:sections) A[(rows[i - 1] + 1):rows[i]] <- Nl[i] - 1
  resist = c(X) - A
  resist <- length(resist[resist == 1])

  # Scale resistance before returning the number of 1's
  # XXX resist = startingFitness + (resist/total.genes)*(1 - startingFitness)
  resist = (resist/total.genes)
  return(resist)
}

sumSection <- function(gtypeSec, ng, nl) {
  if(nl == 1) sum(gtypeSec)
  else sum(tapply(gtypeSec, rep(1:ng, each=nl), prod))
}


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

reportGeneration <- function(P, Rd, sectionRd) {
  Psize <- nrow(P)
  avgFit <- mean(Rd)
  avgSectionFit = c(apply(sectionRd, 2, mean))
  return(c(Psize, avgFit, avgSectionFit))

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

# Receive a vector of locus genotypes to calculate section ressiatnce values
# Receive Nl = loci per gene per section and Ng = genes per section
# Return resistance gained for each section of 'gtype' = fraction of fully mutated genes per section
getSectionResistance0 <- function(gtype, Nl, Ng) {

  if(length(Nl) != length(Ng))
    stop("Number of sections in regions and loci are different!")

  sections <- length(Nl)
  rows <- cumsum(Ng)  # 10  20  30
  rrows <- rep(1:rows[1], each = Nl[1])
  for(i in 2:sections)
    rrows <- c(rrows, rep((rows[i-1]+1):(rows[i]), each = Nl[i]))

  total.genes <- sum(Ng)
  Z <- diag(total.genes)
  Z <- Z[rrows,]
  X <- gtype %*% Z
  A <- numeric(total.genes)

  A[1:rows[1]] <- Nl[1]-1
  for(i in 2:sections)
    A[(rows[i-1]+1):rows[i]] <- Nl[i]-1

  resist <- c(X)-A
  sectionResist <- numeric(sections)

  part = resist[1:rows[1]]
  sectionResist[1] <- length(part[part == 1])

  for(i in 2:sections) {
    part = resist[(rows[i-1]+1) : rows[i]]
    sectionResist[i] <- length(part[part == 1])
  }

  # return fraction of fully mutated genes per section
  return(sectionResist/Ng)
}


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

nextDaySample <- function(P, Rd, Psize) {
  Rd <- Rd/sum(Rd)
  include <- sample.int(length(Rd), size=Psize, prob=Rd)
  return(P[include,])
}
