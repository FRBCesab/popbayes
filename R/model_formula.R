#' @title Write bugs model code
#'
#' @description
#' ...
#'
#' @param jags [bolean] If TRUE, write model code for JAGS. If FALSE, write model code for OPENBUGS.
#'
#' @author Nicolas CASAJUS, \email{nicolas.casajus@@fondationbiodiversite.com}
#' @author Roger PRADEL, \email{roger.pradel@@cefe.cnrs.fr}
#'
#' @export
#'
#' @return
#' For internal use only. No return, write on hard drive a textfile with bugs model code.
#'
#' @examples
#'
#' model_formula(jags = TRUE)
#'
#' model_formula(jags = FALSE)

model_formula <- function(jags = TRUE) {

  if (!is.logical(jags)) {

    stop("`jags` argument must be a boolean.")
  }

  if (jags) {

    header <- "
  data {
    c1 <- c[1]
  }
  "
  } else {

    header <- ""

  }

  core <- "
  model {
    # initialisations of level of Confidence Intervals (could be allowed in entry for more flexibility)
    # currently 95% CI are provided in entry

    for (i in 1:k) {
      sd[i] <- (h[i]-l[i])/3.93
      prec[i] <- pow(sd[i],-2)
    }

    for (i in 1:(k-1)) {
      lint[i] <- t[i+1]-t[i]
    }

    linttot <- t[k]-t[1]
    minN1 <- c1 / 2
    maxN1 <- c1 * 2

    # Priors and constraints
    N[1] ~ dunif(minN1, maxN1)            # Prior for initial population size (N1 inconnue)

    # Likelihood
    # State process
    logN[1]<-log(N[1])
    rcand[1] ~ dnorm(0, 1)
    r[1] <- min(rcand[1], rmax) # r cannot exceed rmax
    rcum[1] <- lint[1]*r[1]
    r2cum[1] <- rcum[1]*r[1]
    logN[2] <- logN[1]+rcum[1]
    N[2] <- exp(logN[2])

    for (i in 2:(k-1)){
      rcand[i] ~ dnorm(r[i-1],lability) # conditions in year i should resemble those in year i-1
      r[i] <- min(rcand[i], rmax) # but cannot exceed rmax
      rcum[i] <- lint[i]*r[i]  # lint = longueur de l'intervalle : nb annees
      r2cum[i] <- rcum[i]*r[i]
      logN[i+1] <- logN[i]+rcum[i]
      N[i+1] <- exp(logN[i+1])
    }

    # Observation process
    for (i in 1:k) {
      c[i] ~ dnorm(N[i],prec[i])
    }

    # derived interesting quantities
    meanr <- sum(rcum[])/linttot  # overall trend
    sdr <- sqrt(sum(r2cum[])/linttot-meanr*meanr) # temporal variation in r
    # cv <- sdr/meanr
    vrrmax <- sdr/rmax
  }
  "

  txt <- c(header, core)
  txt <- paste0(txt, collapse = "\n")

  sink("model.txt")
  cat(txt, fill = TRUE)
  sink()
}
