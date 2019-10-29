#' @title Run BUGS population model
#'
#' @description
#' ...
#'
#' @param dsname [string] the R object name containing data.
#' @param MCMCopt [list] a list containing the number of iteration (ni), thin factor (nt), length of burn in (nb), i.e. number of iterations to discard at the beginning and the number of chains (nc).
#' @param modelopt [list] a list of two vectors: the model smoothing factor and a boolean indicating if parameter r must be limited by rmax.
#' @param jags [bolean] If TRUE, write model bugs code for JAGS, otherwise for OpenBUGS.
#'
#' @author Nicolas CASAJUS, \email{nicolas.casajus@@fondationbiodiversite.fr}
#' @author Roger PRADEL, \email{roger.pradel@@cefe.cnrs.fr}
#'
#' @export
#'
#' @return
#' ...
#'
#' @examples
#'
#' # No example



fit_trend <- function(
  dsname   = NULL,
  rmax     = NULL,
  MCMCopt  = list(ni = 50000, nt = 3, nb = 10000, nc = 2),
  modelopt = list(100, TRUE),
  jags     = TRUE) {

  # remove NA (counts + time)
  # require
  # NULL
  # colnames
  # names list
  # logical

  model_formula(jags = jags)

  data <- get(dsname)
  year <- data$year
  k    <- nrow(data)

  if (modelopt[[2]]) {

    rmax <- log(1 + rmax)

  } else {

    rmax <- 10000
  }

  ni <- MCMCopt[["ni"]]
  nt <- MCMCopt[["nt"]]
  nb <- MCMCopt[["nb"]]
  nc <- MCMCopt[["nc"]]

  data_bugs <- list(
    c        = as.vector(data$counts),
    h        = data$csup,
    l        = data$cinf,
    t        = year,
    k        = k,
    rmax     = rmax,
    lability = modelopt[[1]]
  )

  parameters <- c("N", "r", "meanr", "sdr", "vrrmax")


  if (jags) {

    output_bugs <- R2jags::jags(
      data               = data_bugs,
      inits              = inits_plus(data_bugs, nc),
      parameters.to.save = parameters,
      model.file         = "bugs_model.txt",
      n.chains           = nc,
      n.thin             = nt,
      n.iter             = ni,
      n.burnin           = nb,
      working.directory  = getwd()
    )

  } else {

    output_bugs <- R2OpenBUGS::bugs(
      data              = data_bugs,
      inits             = inits_plus(data_bugs, nc),
      parameters        = parameters,
      model.file        = "bugs_model.txt",
      n.chains          = nc,
      n.thin            = nt,
      n.iter            = ni,
      debug             = TRUE,
      n.burnin          = nb,
      working.directory = getwd()
    )

  }

  print(output_bugs, digits = 3)

  xxx <- as.data.frame(output_bugs$BUGSoutput$summary)
  xxx <- data.frame(
    param = rownames(xxx),
    xxx
  )

  save(
    output_bugs,
    file = paste0(
      paste(
        dsname,
        modelopt[[1]],
        modelopt[[2]],
        sep = "_"
      ),
      ".RData"
    )
  )

  return(output_bugs)
}
