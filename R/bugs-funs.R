#' Generate initial values of relative rate of increase (r)
#'
#' @description
#' This function extracts the counts and the corresponding dates from the list 
#' of objects (`data`) to be analyzed by the Bayesian model. It calculates the 
#' observed annual relative rate of increase (_r\_cand\_obs_) for each > 0
#' interval as the ratio of the change in log count divided by the interval 
#' length. It then draws a random initial value of the candidate relative rate 
#' of increase (_rcand_) to be passed to the bugs model as initial value by 
#' drawing in a normal distribution centered on the observed rate  
#' (_r\_cand\_obs_) with a standard deviation of 1. For each _r\_cand\_obs_, 
#' `nc` random values are drawn where `nc` (second argument) is the number of 
#' Markov chains.
#'
#' @param data a `list` created by [fit_trend()]
#' 
#' @param nc a positive `integer`. The number of MCMC chains.
#'
#' @return A `list` of `nc` lists (one per chain). Each list has exactly one 
#' element `rcand`, i.e. a vector of random initial values of candidate 
#' relative rates of increase, one per > 0 interval.    
#'
#' @noRd

inits_plus <- function(data, nc) {
  
  if (!is.list(data)) {
    stop("Argument 'data' must be an output of format_data().")
  }
  
  if (missing(nc)) {
    stop("Argument 'nc' is required.")
  }
  
  if (!is.numeric(nc) || length(nc) != 1) {
    stop("Argument 'nc' must be an integer of length 1.")
  }
  
  if (nc < 1) {
    stop("Argument 'nc' must be a positive integer.")
  }
  
  
  if (!("c" %in% names(data))) {
    stop("Missing 'c' element (count) in 'data'.")
  }
  
  if (!("t" %in% names(data))) {
    stop("Missing 't' element (time) in 'data'.")
  }
  
  if (!("k" %in% names(data))) {
    stop("Missing 'k' element (number of points) in 'data'.")
  }
  
  
  if (length(data$"c") < 4) {
    stop("Series must have at least 4 points.")
  }
  
  if (sum(which(data$"c" < 0)) > 0) {
    stop("Count must be positive or zero.")
  }
  
  if (sum(which(data$"t" <= 0)) > 0) {
    stop("Time must be strictly positive (numerical date).")
  }
  
  
  r_cand_obs <- log(data$"c"[2:data$"k"]) - log(data$"c"[1:(data$"k" - 1)])
  l_interval <- data$"t"[2:data$"k"] - data$"t"[1:(data$"k" - 1)]
  
  r_cand_obs[l_interval != 0] <- r_cand_obs[l_interval != 0] / 
    l_interval[l_interval != 0]
  
  list_start <- list()
  
  for (i in 1:nc) {
    list_start[[i]] <- list("rcand" = stats::rnorm(length(r_cand_obs), 
                                                   r_cand_obs, 1))
  }
  
  list_start
}



#' Formulation of the Bayesian model to be used by JAGS
#'
#' @description
#' The code describes a Bayesian model of the counts of one species/site over 
#' time. The model attempts to reconstitute the trend of the population size 
#' while accounting for the precision of the individual counts and the 
#' demographic potential of the species. It also implements a smoothing rule of 
#' the relative rates of increase in successive dates 
#' 
#' The code generated by this function is saved in a text file later read by 
#' the JAGS program.
#'
#' @param path a `character` string. The directory to save the file. This 
#'   directory must exist and can be an absolute or a relative path.
#'
#' @return No return value.
#' 
#' @noRd

model_formula <- function(path = ".") {
  
  if (!dir.exists(path)) {
    stop("The directory '", path, "' does not exist.")
  }
  
  
  time   <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  
  header <- "\n###    Bayesian Model Formula    ###\n\n"
  header <- c(header, "# ----------------------------------")
  
  header <- c(header, "# Program: JAGS")
  header <- c(header, paste0("# Time   : ", time, ""))
  header <- c(header, "# ------------------------------")
  header <- c(header, "\n\n")
  header <- c(header, "# Data formulation\n")
  header <- c(header, "data {\n")
  header <- c(header, "\tc1 <- c[1]")
  header <- c(header, "}")
  header <- c(header, "\n\n")
  
  
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
  
  cat(txt, file = file.path(path, "model.bug"), append = FALSE, 
      fill = TRUE)
  
  invisible(NULL)
}



#' Fit a Bayesian model to estimate population size trend
#'
#' @description
#' This function applies a Bayesian model to count series in order to infer 
#' the population trend over time. This function only works on the output of 
#' [format_data()] or [filter_series()].
#' 
#' **Important:** This function uses [R2jags::jags()] and the 
#' freeware **JAGS** (\url{https://mcmc-jags.sourceforge.io/}) must be 
#' installed.
#' 
#' There are two types of options: model options (argument `model_opts`) and 
#' MCMC options (argument `mcmc_opts`).
#' 
#' 
#' **A. Model options**
#' 
#' - a smoothing factor: the precision (the inverse of variance) of a normal 
#' distribution centered on the current relative rate of increase r from which 
#' the next candidate relative rate of increase (see below) is drawn. 
#' The highest this number, the tighter the link between successive relative 
#' rates of increase. The default `100` corresponds to a moderate link.
#' 
#' - a logical indicating whether the population growth must remain limited by 
#' the species demographic potential (provided by the argument `rmax` in 
#' [format_data()]).
#' 
#' The relative rate of increase is the change in log population size between 
#' two dates. The quantity actually being modeled is the relative rate of 
#' increase per unit of time (usually one date). This quantity reflects more 
#' directly the prevailing conditions than the population size itself, which is 
#' the reason why it has been chosen.
#' 
#' When the second model option is set to `TRUE`, the candidate rate of 
#' increase is compared to the maximum relative rate of increase 
#' (obtained when using [format_data()]) and replaced by `rmax` if greater.
#' 
#' 
#' **B. MCMC options**
#' 
#' Classical Markov chain Monte Carlo (MCMC) settings (see argument `mcmc_opts`
#' below).
#'
#' @param data a named `list`. The output of [format_data()] or 
#'   [filter_series()].
#' 
#' @param model_opts a `list` of two vectors. The model smoothing factor 
#'   (`numeric`) and a `logical` indicating if the parameter **r** must be 
#'   limited by `rmax` (`TRUE`) or not (`FALSE`). If this second parameter is 
#'   `TRUE`, the argument `rmax` cannot be `NULL` unless species are listed in 
#'   `popbayes` (in `species_info`).
#' 
#' @param mcmc_opts a `list` containing the number of iterations (`ni`), the  
#'   thinning factor (`nt`), the length of burn in (`nb`), i.e. the number of 
#'   iterations to discard at the beginning, and the number of chains (`nc`).
#'   
#' @param path a `character` string. The directory to save BUGS outputs (the 
#'   same as in [format_data()]). 
#'   
#' @export
#'
#' @return An n-element `list` (where `n` is the number of count series). Each
#'   element of the list is a BUGS output as provided by JAGS (also written 
#'   in the folder `path`).
#'
#' @examples
#' ## Load Garamba raw dataset ----
#' file_path <- system.file("extdata", "garamba_survey.csv", 
#'                          package = "popbayes")
#'                          
#' garamba <- read.csv(file = file_path)
#' 
#' ## Create temporary folder ----
#' temp_path <- tempdir()
#' 
#' ## Format dataset ----
#' garamba_formatted <- popbayes::format_data(
#'   data              = garamba, 
#'   path              = temp_path,
#'   field_method      = "field_method",
#'   pref_field_method = "pref_field_method",
#'   conversion_A2G    = "conversion_A2G",
#'   rmax              = "rmax")
#' 
#' ## Get data for Alcelaphus buselaphus at Garamba only ----
#' a_buselaphus <- popbayes::filter_series(garamba_formatted, 
#'                                         location = "Garamba",
#'                                         species  = "Alcelaphus buselaphus")
#'                                         
#' \donttest{
#' ## Fit population trend (requires JAGS) ----
#' a_buselaphus_mod <- popbayes::fit_trend(a_buselaphus, path = temp_path)
#' 
#' ## Check for convergence ----
#' popbayes::diagnostic(a_buselaphus_mod, threshold = 1.1)
#' 
#' ## Plot estimated population trend ----
#' popbayes::plot_trend(series = "garamba__alcelaphus_buselaphus", 
#'                      path   = temp_path)
#' 
#' ## Plot MCMC traceplot ----
#' R2jags::traceplot(a_buselaphus_mod[[1]], ask = TRUE)
#' }

fit_trend <- function(data, model_opts = list(100, TRUE), 
                      mcmc_opts = list(ni = 50000, nt = 3, nb = 10000, nc = 2),
                      path = ".") {
  
  
  ## Check datasets ----
  
  if (!is.list(data)) {
    stop("Argument 'data' must be an output of format_data().")
  }
  
  if (!("data_converted" %in% names(data[[1]]))) {
    stop("Argument 'data' must be an output of format_data().")
  }
  
  
  ## Check path ----
  
  if (!dir.exists(path)) {
    stop("The directory '", path, "' does not exist.")
  }
  
  
  ## Check model options ----
  
  if (!is.list(model_opts)) {
    stop("Argument 'model_opts' must be a list.")
  }
  
  if (length(model_opts) != 2) {
    stop("Argument 'model_opts' must be a list of length 2.")
  }
  
  if (!is.numeric(model_opts[[1]]) || length(model_opts[[1]]) != 1) {
    stop("The first element of the list 'model_opts' must be numeric of ", 
         "length 1.")
  }
  
  if (!is.logical(model_opts[[2]]) || length(model_opts[[2]]) != 1) {
    stop("The second element of the list 'model_opts' must be TRUE or FALSE.")
  }
  
  
  ## Check MCMC options ----
  
  if (!is.list(mcmc_opts)) {
    stop("Argument 'model_opts' must be a list.")
  }
  
  if (length(mcmc_opts) != 4) {
    stop("Argument 'model_opts' must be a list of length 4.")
  }
  
  valid_mcmc_opts_names <- c("ni", "nt", "nb", "nc")
  valid_mcmc_opts_names_msg <- paste0(valid_mcmc_opts_names, 
                                      collapse = "' and '")
  valid_mcmc_opts_names_msg <- paste0("'", valid_mcmc_opts_names_msg, "'")
  
  if (any(!(names(mcmc_opts) %in% valid_mcmc_opts_names))) {
    stop("Invalid parameters in 'mcmc_opts'. ",
         "Required parameters are: ", valid_mcmc_opts_names_msg, ".")
  }
  
  if (!is.numeric(mcmc_opts[[1]]) || length(mcmc_opts[[1]]) != 1) {
    stop("The first element of the list 'mcmc_opts' must be numeric of ", 
         "length 1.")
  }
  
  if (!is.numeric(mcmc_opts[[2]]) || length(mcmc_opts[[2]]) != 1) {
    stop("The second element of the list 'mcmc_opts' must be numeric of ", 
         "length 1.")
  }
  
  if (!is.numeric(mcmc_opts[[3]]) || length(mcmc_opts[[3]]) != 1) {
    stop("The third element of the list 'mcmc_opts' must be numeric of ", 
         "length 1.")
  }
  
  if (!is.numeric(mcmc_opts[[4]]) || length(mcmc_opts[[4]]) != 1) {
    stop("The fourth element of the list 'mcmc_opts' must be numeric of ", 
         "length 1.")
  }
  
  
  params_to_infer <- c("N", "r", "meanr", "sdr", "vrrmax")
  
  bugs_list <- list()
  
  
  ## Loop on count series ----
  
  for (i in seq_len(length(data))) {
    
    
    ## Create sub-directory ----
    
    series_name  <- names(data)[i]
    species_path <- file.path(path, series_name)
    dir.create(species_path, showWarnings = FALSE)
    
    
    ## Model formulation ----
    
    model_formula(species_path)
    
    
    ## Prepare data ----
    
    data_sub <- data[[i]]
    
    species_rmax <- ifelse(model_opts[[2]], log(1 + data_sub$"rmax"), 10000)
    
    data_for_jags <- list(
      c        = data_sub$"data_converted"$"count_conv",
      h        = data_sub$"data_converted"$"upper_ci_conv",
      l        = data_sub$"data_converted"$"lower_ci_conv",
      t        = data_sub$"data_converted"$"date",
      k        = nrow(data_sub$"data_converted"),
      rmax     = species_rmax,
      lability = model_opts[[1]]
    )
    
    
    ## Fit population trend model ----
    
    bugs_outputs <- R2jags::jags(
      data               = data_for_jags,
      inits              = inits_plus(data_for_jags, mcmc_opts[["nc"]]),
      parameters.to.save = params_to_infer,
      model.file         = "model.bug",
      n.chains           = mcmc_opts[["nc"]],
      n.thin             = mcmc_opts[["nt"]],
      n.iter             = mcmc_opts[["ni"]],
      n.burnin           = mcmc_opts[["nb"]],
      working.directory  = species_path)
    
    invisible(file.remove(file.path(species_path, "model.bug")))
    
    bugs_list[[series_name]] <- bugs_outputs
    bugs_outputs <- bugs_list[series_name]
    
    save(bugs_outputs, file = file.path(species_path, paste0(series_name, 
                                                             "_bugs.RData")))
  }
  
  bugs_list
}



#' Import a list of BUGS outputs previously exported
#'
#' @description 
#' This function imports a list of BUGS outputs previously exported by 
#' [fit_trend()]. Users can import one, several, or all models.
#'
#' @param series a vector of `character` strings. One or several count series 
#'   names. If `NULL` (default) BUGS outputs for all count series will be 
#'   imported. Users can run [list_series()] to get the correct spelling of 
#'   count series names.
#'    
#' @param path a `character` string. The directory in which BUGS outputs have 
#'   been saved by the function [fit_trend()].
#'
#' @return An n-element `list` (where `n` is the number of count series). See
#'   [fit_trend()] for further information.
#' 
#' @export
#'
#' @examples
#' ## Load Garamba raw dataset ----
#' file_path <- system.file("extdata", "garamba_survey.csv", 
#'                          package = "popbayes")
#'                          
#' garamba <- read.csv(file = file_path)
#' 
#' ## Create temporary folder ----
#' temp_path <- tempdir()
#' 
#' ## Format dataset ----
#' garamba_formatted <- popbayes::format_data(
#'   data              = garamba, 
#'   path              = temp_path,
#'   field_method      = "field_method",
#'   pref_field_method = "pref_field_method",
#'   conversion_A2G    = "conversion_A2G",
#'   rmax              = "rmax")
#' 
#' ## Select one serie ----
#' a_buselaphus <- popbayes::filter_series(garamba_formatted, 
#'                                         location = "Garamba",
#'                                         species  = "Alcelaphus buselaphus")
#' \donttest{
#' ## Fit population trends (requires JAGS) ----
#' a_buselaphus_mod <- popbayes::fit_trend(a_buselaphus, path = temp_path)
#' 
#' ## Import BUGS outputs for one count series ----
#' popbayes::read_bugs(series = "garamba__alcelaphus_buselaphus", 
#'                     path   = temp_path)
#' 
#' ## Import BUGS outputs for all count series ----
#' popbayes::read_bugs(path = temp_path)
#' }

read_bugs <- function(series = NULL, path = ".") {
    
    
  if (!dir.exists(path)) {
    stop("The directory '", path, "' does not exist.")
  }
  
  filenames <- list.files(path, recursive = TRUE, pattern = "_bugs\\.RData")
  
  if (length(filenames) == 0) {
    stop("No BUGS outputs can be found.")  
  }
  
  
  ## All available series names ----
  
  series_names <- strsplit(filenames, .Platform$"file.sep")
  series_names <- unlist(lapply(series_names, function(x) x[1]))
  
  
  if (!is.null(series)) {
    
    if (!is.character(series)) {
      stop("Argument 'series' must be a character (series name(s)).")
    }
    
    if (any(!(series %in% series_names))) {
      stop("Some count series cannot be found.")
    } 
    
    series_names <- series
  }
  
  data_bugs <- list()
  
  for (series in series_names) {
    
    data_bugs <- c(data_bugs, get(load(file.path(path, series, 
                                                 paste0(series, 
                                                        "_bugs.RData")))))
  }
  
  data_bugs
}



#' Extract estimated parameters from a list of BUGS outputs
#'
#' @description
#' From the output of the function [fit_trend()] (or [read_bugs()]), this 
#' function extracts estimated parameters into a `data.frame`.
#' 
#' The resulting `data.frame` has no particular use in `popbayes` but it can be
#' useful for users.
#'
#' @param data a named `list` of BUGS outputs. The output of [fit_trend()] or
#'   [read_bugs()]
#'
#' @return A `data.frame`.
#' 
#' @export
#'
#' @examples
#' ## Load Garamba raw dataset ----
#' file_path <- system.file("extdata", "garamba_survey.csv", 
#'                          package = "popbayes")
#'                          
#' garamba <- read.csv(file = file_path)
#' 
#' ## Create temporary folder ----
#' temp_path <- tempdir()
#' 
#' ## Format dataset ----
#' garamba_formatted <- popbayes::format_data(
#'   data              = garamba, 
#'   path              = temp_path,
#'   field_method      = "field_method",
#'   pref_field_method = "pref_field_method",
#'   conversion_A2G    = "conversion_A2G",
#'   rmax              = "rmax")
#'                                         
#' ## Select one serie ----
#' a_buselaphus <- popbayes::filter_series(garamba_formatted, 
#'                                         location = "Garamba",
#'                                         species  = "Alcelaphus buselaphus")
#' \donttest{
#' ## Fit population trends (requires JAGS) ----
#' a_buselaphus_mod <- popbayes::fit_trend(a_buselaphus, path = temp_path)
#' 
#' ## Import BUGS outputs for one count series ----
#' bugs <- popbayes::read_bugs(series = "garamba__alcelaphus_buselaphus", 
#'                             path   = temp_path)
#' 
#' ## Extract estimated parameters ----
#' popbayes::bugs_to_df(bugs)
#' }

bugs_to_df <- function(data) {
  
  if (!is.list(data)) {
    stop("Argument 'data' must be an output of fit_trend().")
  }
  
  if (!("BUGSoutput" %in% names(data[[1]]))) {
    stop("Argument 'data' must be an output of fit_trend().")
  }
  
  
  series_names <- names(data)
  
  bugs <- lapply(data, function(x) x[["BUGSoutput"]]$"summary")
  bugs <- lapply(seq_len(length(series_names)), function(x, y, z) {
    bug <- cbind(z[x], rownames(y[[x]]), as.data.frame(y[[x]]))
    colnames(bug)[1:2] <- c("series", "parameter")
    rownames(bug) <- NULL
    bug
  }, y = bugs, z = series_names)
  
  bugs <- do.call(rbind.data.frame, bugs)
  rownames(bugs) <- NULL
  
  bugs
}



#' Check if a BUGS model has converged
#'
#' @description
#' From the output of the function [fit_trend()] (or [read_bugs()]), this 
#' function checks if the estimation of all parameters of one (or several) 
#' BUGS model has converged. This diagnostic is performed by comparing the 
#' `Rhat` value of each parameter to a `threshold` (default is `1.1`). If some 
#' `Rhat` values are greater than this threshold (no convergence), a message
#' listing problematic models is displayed.
#'
#' @param data a named `list` of BUGS outputs. The output of [fit_trend()] or 
#'   [read_bugs()].
#'   
#' @param threshold a `numeric`.
#'
#' @return No return value.
#' 
#' @export
#'
#' @examples
#' ## Load Garamba raw dataset ----
#' file_path <- system.file("extdata", "garamba_survey.csv", 
#'                          package = "popbayes")
#'                          
#' garamba <- read.csv(file = file_path)
#' 
#' ## Create temporary folder ----
#' temp_path <- tempdir()
#' 
#' ## Format dataset ----
#' garamba_formatted <- popbayes::format_data(
#'   data              = garamba, 
#'   path              = temp_path,
#'   field_method      = "field_method",
#'   pref_field_method = "pref_field_method",
#'   conversion_A2G    = "conversion_A2G",
#'   rmax              = "rmax")
#'                                         
#' ## Select one serie ----
#' a_buselaphus <- popbayes::filter_series(garamba_formatted, 
#'                                         location = "Garamba",
#'                                         species  = "Alcelaphus buselaphus")
#' \donttest{
#' ## Fit population trends (requires JAGS) ----
#' a_buselaphus_mod <- popbayes::fit_trend(a_buselaphus, path = temp_path)
#' 
#' ## Check for convergence ----
#' popbayes::diagnostic(a_buselaphus_mod)
#' }

diagnostic <- function(data, threshold = 1.1) {
  
  if (!is.list(data)) {
    stop("Argument 'data' must be an output of fit_trend().")
  }
  
  if (!("BUGSoutput" %in% names(data[[1]]))) {
    stop("Argument 'data' must be an output of fit_trend().")
  }
  
  
  series_names <- names(data)
  
  bugs <- lapply(data, function(x) x[["BUGSoutput"]]$"summary")
  bugs <- lapply(bugs, function(x) {
    length(which(x[ , "Rhat"] > threshold))
  })
  
  bugs <- unlist(bugs)
  
  pos <- which(bugs > 0)
  
  if (length(pos)) {
    
    no_converged <- names(bugs)[pos]
    cli::cli_alert_danger(c(
      "Some parameters have not converged for the following count series: ",
      "{.val {no_converged}}."
    ))
    
  } else {
    
    cli::cli_alert_success("All models have converged.")
  }
  
  invisible(NULL)
}
