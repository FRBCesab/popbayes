#' @title Run BUGS population model
#'
#' @description
#' This is the main function of the package. It applies a Bayesian model to a count series in order to infer 
#' the population trend over time. 
#' The series of counts, extracted from a tibble or dataframe, the name of which is provided as an argument 
#' (parameter dsname), must be accompanied by a 95% confidence interval in the form of 2 parallel series of the  
#' lower and upper bounds, and the series of the times when each count was taken (in the same tibble or dataframe).
#' 
#' There are two types of options: 
#' - model options (parameter modelopt)
#' - MCMC options (parameter MCMCopt)
#' 
#' A. model options are 2: 
#' - a smoothing parameter, 
#' - a logical indicating whether the population growth must remain limited by the species demographic potential. 
#' 
#' The smoothing parameter is actually the precision (the inverse of variance) of a normal distribution centered on 
#' the current relative rate of increase r from which the next candidate relative rate of increase (see below) is drawn. 
#' The highest this number, the tighter the link between successive relative rates of increase. The default 100
#' corresponds to a moderate link.
#' 
#' The relative rate of increase is the change in log population size between two dates. The quantity actually 
#' being modeled is the relative rate of increase per unit of time (usually one year). This quantity reflects 
#' more directly the prevailing conditions than the population size itself, which is the reason why it has been chosen.
#' When the second model option is set to True, the candidate rate of increase is compared to the maximum relative rate 
#' of increase (parameter rmax) and replaced by rmax if greater.
#' 
#' B. MCMC options are the classical MCMC settings (see param MCMCopt below).
#'
#' @param dsname [string] the R object name containing data.
#' @param rmax [numeric] the maximum relative rate of increase of the species i.e. maximum yearly change in log pop size.
#' @param MCMCopt [list] a list containing the number of iteration (ni), thin factor (nt), length of burn in (nb), i.e. number of iterations to discard at the beginning and the number of chains (nc).
#' @param modelopt [list] a list of two vectors: the model smoothing factor and a boolean indicating if parameter r must be limited by rmax.
#' @param jags [boolean] If TRUE, write model bugs code for JAGS, otherwise for OpenBUGS.
#'
#' @author Nicolas CASAJUS, \email{nicolas.casajus@@fondationbiodiversite.fr}
#' @author Roger PRADEL, \email{roger.pradel@@cefe.cnrs.fr}
#'
#' @export
#'
#' @return
#' the ouput of the Bayesian model as provided by Jags or Openbugs. Additionally, the function prints the results 
#' to the console with the method of the package used to call JAGS or OpenBUGS.
#'
#' @examples
#'
#' # No example



fit_trend <- function(dsname = NULL, rmax = NULL, 
                      MCMCopt = list(ni = 50000, nt = 3, nb = 10000, nc = 2),
                      modelopt = list(100, TRUE),
                      jags = TRUE) {

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
    
    if (is.null(rmax)) stop("rmax is required")

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
