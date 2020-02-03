#' @title Write bugs model code
#'
#' @description
#' ...
#'
#' @param jags [boolean] If TRUE, write model bugs code for JAGS, otherwise for OpenBUGS.
#'
#' @author Nicolas CASAJUS, \email{nicolas.casajus@@fondationbiodiversite.fr}
#' @author Roger PRADEL, \email{roger.pradel@@cefe.cnrs.fr}
#'
#' @export
#'
#' @return
#' For internal use only.
#' No return. Write on hard drive a textfile with bugs model code.
#'
#' @examples
#'
#' model_formula(jags = TRUE)   # JAGS formulation
#' model_formula(jags = FALSE)  # OpenBUGS formulation



model_formula <- function(jags = TRUE) {

  if (!is.logical(jags)) { stop("`jags` argument must be a boolean.") }


  time   <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")


  header <- "\n###    BUGS Model Formula    ###\n\n"
  header <- c(header, "# ------------------------------")

  if (jags) {

    header <- c(header, "# Program: JAGS")
    header <- c(header, paste0("# Time   : ", time, ""))
    header <- c(header, "# ------------------------------")
    header <- c(header, "\n\n")
    header <- c(header, "# Data formulation\n")
    header <- c(header, "data {\n")
    header <- c(header, "\tc1 <- c[1]")
    header <- c(header, "}")
    header <- c(header, "\n\n")

  } else {

    header <- c(header, "# Program: OpenBUGS")
    header <- c(header, paste0("# Time   : ", time, ""))
    header <- c(header, "# ------------------------------")
    header <- c(header, "\n")

  }


  core <- ""

  core <- c(core, "# Model formulation\n")
  core <- c(core, "model {\n")

  core <- c(core, "\t# Initialisation of level of CI95%")
  core <- c(core, "\tfor (i in 1:k) {")
  core <- c(core, "\t\tsd[i]   <- (h[i] - l[i]) / 3.93")
  core <- c(core, "\t\tprec[i] <- pow(sd[i], -2)")
  core <- c(core, "\t}\n")

  core <- c(core, "\tfor (i in 1:(k - 1)) {")
  core <- c(core, "\t\tlint[i] <- t[i + 1] - t[i]")
  core <- c(core, "\t}\n")
  core <- c(core, "\tlinttot <- t[k] - t[1]\n")

  if (!jags) {
    
    core <- c(core, "\tc1 <- c[1]")

  } 
    
  core <- c(core, "\tminN1 <- c1 / 2")
  core <- c(core, "\tmaxN1 <- c1 * 2\n")

  core <- c(core, "\t# Priors and constraints")
  core <- c(core, "\tN[1] ~ dunif(minN1, maxN1)\n")

  core <- c(core, "\t# Likelihood - State process")
  core <- c(core, "\tlogN[1]  <- log(N[1])")
  core <- c(core, "\trcand[1] ~ dnorm(0, 1)")
  core <- c(core, "\tr[1]     <- min(rcand[1], rmax)")
  core <- c(core, "\trcum[1]  <- lint[1] * r[1]")
  core <- c(core, "\tr2cum[1] <- rcum[1] * r[1]")
  core <- c(core, "\tlogN[2]  <- logN[1] + rcum[1]")
  core <- c(core, "\tN[2]     <- exp(logN[2])\n")

  core <- c(core, "\tfor (i in 2:(k - 1)) {")
  core <- c(core, "\t\trcand[i] ~ dnorm(r[i - 1], lability)")
  core <- c(core, "\t\tr[i]        <- min(rcand[i], rmax)")
  core <- c(core, "\t\trcum[i]     <- lint[i] * r[i]")
  core <- c(core, "\t\tr2cum[i]    <- rcum[i] * r[i]")
  core <- c(core, "\t\tlogN[i + 1] <- logN[i] + rcum[i]")
  core <- c(core, "\t\tN[i + 1]    <- exp(logN[i + 1])")
  core <- c(core, "\t}\n")

  core <- c(core, "\t# Likelihood - Observation process")
  core <- c(core, "\tfor (i in 1:k) {")
  core <- c(core, "\t\tc[i] ~ dnorm(N[i], prec[i])")
  core <- c(core, "\t}\n")

  core <- c(core, "\t# derived interesting quantities")
  core <- c(core, "\tmeanr  <- sum(rcum[]) / linttot")
  core <- c(core, "\tsdr    <- sqrt(sum(r2cum[]) / linttot - meanr * meanr)")
  core <- c(core, "\tvrrmax <- sdr / rmax\n")

  core <- c(core, "}\n")

  txt <- c(header, core)
  txt <- paste0(txt, collapse = "\n")

  cat(txt, file = "bugs_model.txt", append = FALSE, fill = TRUE)
}

