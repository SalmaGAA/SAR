#'
#' @title Plot Scenario for Bacterial Genotype
#'
#' @author Ahmad Abdel-Azim \email {agabdel.azim@gmail.com}
#' @author Salma Abdel-Azim \email {salma.abdelazim@gmail.com}
#' @author Gamal Abdel-Azim \email {gamal.azim@gmail.com}
#'
#' @name plotScenario
#'
#' @usage plotScenario(scenario, gen.pd, nDay, Nl, Ng, ylim)
#'
#' @param scenario ADD LATER
#' @param gen.pd ADD LATER
#' @param nDay ADD LATER
#' @param Nl Number of loci per gene per section. One-dimensional array of length equal to number of sections.
#' @param Ng Number of genes per section. One-dimensional array of length equal to number of sections.
#' @param ylim ADD LATER
#'
#' @description ADD LATER
#'
#' @details ADD LATER
#'
#' @return ADD LATER
#'

#' @examples
#' plotScenario(scenario = scenarioSimulate, gen.pd = 24, nDay = 20, Nl = c(1, 2, 3), Ng = c(10, 10, 10), ylim = 0.5, 0.85)

#'@keywords SAR
#'@keywords resistance
#'@keywords antibiotic
#'@keywords simulation


#' @export
plotScenario <- function(scenario, gen.pd, nDay, Nl, Ng, ylim) {
  ## Cleanup zero simulation runs (lost/crashed populations)
  scenario <- cleanUpZeros(scenario)

  ## Flatten 3D array into 2D array and assign column names: ni, fit, fit.1, fit.2, fit.3, SR, Nl1, Nl2, Nl3, Ng1, Ng2, Ng3, Day
  flat.scenario <- flaten2D(scenario, Nl = Nl, Ng = Ng)
  flatScenario$Day <- rep(rep(1:nDay, each = gen.pd), flat.scenario$SR[1])
  flatScenario.last <- flatScenario[seq(gen.pd, nrow(flatScenario), gen.pd), ] ## Keep only last generation per day

  ## plot averages and trend
  par(mfrow = c(1,2))
  plotAvg(flatScenario.last, tt = "Overall and Section Fitness", ylim = ylim)

  legend("bottom", inset = c(0, -0.45), legend=c("Overall Fitness", "Section 1 Fitness", "Section 2 Fitness", "Section 3 Fitness"),
         col=c("dark green", "red", "orange", "yellow"), lty=c(2,NA, NA, NA),  lwd = 2.5, pch=16, cex = 1.5, xpd = NA, ncol = 4, text.width = c(13,13,13,13), x.intersp = .2)

  plotTrend(flatScenario.last, y.lim = ylim, tt = "Overall Trend")

}

## =====================
## Supporting Functions
## =====================

## Cleanup zero matrices from a 3-dim'l array
cleanUpZeros <- function(ar) {
  # ar: 3-dim'l array
  i0 <- which(apply(ar, 3, function(x) {mean(x[,1])}) == 0)
  if(length(i0)) return(ar[,,-i0])
  else return(ar)
}

## Flatten into 2-dim'l and characterize resulting array
flaten2D <- function(ar, Nl, Ng) {
  successRuns <- dim(ar)[3]
  TwoDimAr <- do.call('rbind', lapply(1:successRuns, function(x) ar[,,x]))
  ## attach number of successful runs, number of loci and number of genes to the 2-dim'l flat array
  TwoDimAr <- as.data.frame(TwoDimAr)
  names(TwoDimAr) <- c('ni', 'fit', 'fit.1', 'fit.2', 'fit.3')
  TwoDimAr$SR <- successRuns
  TwoDimAr$Nl1 <- Nl[1]; TwoDimAr$Nl2 <- Nl[2]; TwoDimAr$Nl3 <- Nl[3]
  TwoDimAr$Ng1 <- Ng[1]; TwoDimAr$Ng2 <- Ng[2]; TwoDimAr$Ng3 <- Ng[3]
  return(TwoDimAr)
}

plotAvg <- function(R, y.lim = c(.4, 1), tt) {
  # average replicates
  days <- sort(unique(R$Day))
  fit <- tapply(R$fit, R$Day, mean)
  fit.1 <- tapply(R$fit.1, R$Day, mean)
  fit.2 <- tapply(R$fit.2, R$Day, mean)
  fit.3 <- tapply(R$fit.3, R$Day, mean)

  plot(days, fit, ylim = y.lim, col='dark green', pch=16, cex = 2, cex.lab = 1.7, cex.axis = 1.5,
       ylab = 'Overall and Section Fitness', xlab = 'Days', type = "b", lwd = 2.5, lty = 2, main=tt, cex.main = 2)
  points(days, fit.1, col='red', pch=16, cex = 1.5)
  points(days, fit.2, col='orange', pch=16, cex = 1.5)
  points(days, fit.3, col='yellow', pch=16, cex = 1.5)
}


plotTrend <- function(R, y.lim = c(.5, .85), tt) {
  # Statistical trend based on a 3rd degree polynomial
  days <- sort(unique(R$Day))

  lm.fitday <- lm(fit ~ poly(Day, 3), R)
  plot(R$Day, R$fit, col='gray', pch=16, cex = 1, ylim = y.lim,  ylab = 'Overall Fitness',
       xlab = 'Days', main = tt, cex.lab = 1.7, cex.axis = 1.5, cex.main = 2)
  lines(days, predict(lm.fitday)[1:length(days)],
        col='navy blue', type = "l", lwd = 4, lty = 1)

  summary(lm.fitday)
}

