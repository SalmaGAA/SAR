#'
#' @title Plot Day Averages for Simulated Scenarios
#'
#' @author Ahmad Abdel-Azim \email{agabdel.azim@gmail.com}
#' @author Salma Abdel-Azim \email{salma.abdelazim@gmail.com}
#' @author Gamal Abdel-Azim \email{gamal.azim@gmail.com}
#'
#' @name plotDayAvg
#'
#' @usage plotDayAvg(scenario, gen.pd, nDay, Nl, Ng, ylim)
#'
#' @param scenario ADD LATER
#' @param gen.pd ADD LATER
#' @param nDay ADD LATER
#' @param Nl Number of loci per gene per section. One-dimensional array of length equal to number of sections.
#' @param Ng Number of genes per section. One-dimensional array of length equal to number of sections.
#' @param ylim ADD LATER
#' @legend.position
#'
#' @description ADD LATER
#'
#' @details ADD LATER
#'
#' @return ADD LATER
#'

#' @examples
#' plotDayAvg(scenario = scenarioSimulate, gen.pd = 24, nDay = 20, Nl = c(1, 2, 3), Ng = c(10, 10, 10), ylim = c(0.5, 0.85))

#'@keywords SAR
#'@keywords resistance
#'@keywords antibiotic
#'@keywords simulation


#' @export
plotDayAvg <- function(scenario, gen.pd, nDay, Nl, Ng, ylim, legend.position = c(0,-0.5)) {
  ## Cleanup zero simulation runs (lost/crashed populations)
  scenario <- cleanUpZeros(scenario)

  ## Flatten 3D array into 2D array and assign column names: ni, fit, fit.1, fit.2, fit.3, SR, Nl1, Nl2, Nl3, Ng1, Ng2, Ng3, Day
  flatScenario <- flaten2D(scenario, Nl = Nl, Ng = Ng)
  flatScenario$Day <- rep(rep(1:nDay, each = gen.pd), flatScenario$SR[1])
  flatScenario.last <- flatScenario[seq(gen.pd, nrow(flatScenario), gen.pd), ] ## Keep only last generation per day

  ## plot averages and trend
  par(mar = c(9, 5, 4, 2))
  plotAvg(flatScenario.last, tt = "Overall and Section Fitness", y.lim = ylim)

  legend("bottom", inset = legend.position, legend=c("Overall Fitness", "Section 1 Fitness", "Section 2 Fitness", "Section 3 Fitness"),
         col=c("dark green", "red", "orange", "yellow"), lty=c(2,NA, NA, NA),  lwd = 2.5, pch=16, cex = 1, xpd = NA, ncol = 2)
}

#'
#'
#'
#'
#'
