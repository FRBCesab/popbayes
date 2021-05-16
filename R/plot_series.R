#' Plot original vs converted counts (with confident interval)
#' 
#' @description
#' From the output of the function [format_data()], this function plots a
#' panel of two figures for one counts series:
#' - on the left side, scatter plot overlapping original and converted counts;
#' - on the right side, scatter plot of converted counts with boundaries of the
#' 95% confident interval.
#' 
#' @param data a named list. The output of the function [format_data()].
#' 
#' @param series a numeric or a character. Either a series name or a series 
#'   index.
#'   
#' @param title a logical. If `TRUE` (default) a title is added at the top-right
#'    of the plot.
#'   
#' @param path a character. The directory to save the plot (if `save = TRUE`). 
#'   This directory must exist and can be an absolute or a relative path.
#'   
#' @param save a logical. If `TRUE` (default is `FALSE`) the plot is saved in 
#'   `path`.
#'
#' @return NULL
#' 
#' @export
#'
#' @examples
#' data(garamba)
#' 
#' ## Format dataset ----
#' garamba_formatted <- format_data(data = garamba)
#' 
#' ## Get series names ----
#' names(garamba_formatted)
#' 
#' ## Plot for Alcelaphus buselaphus at Garamba ----
#' plot_series(data = garamba_formatted,
#'             series = "garamba__alcelaphus_buselaphus", save = FALSE)
#'             
#' ## alternatively...
#' plot_series(data = garamba_formatted, series = 1, save = FALSE)

plot_series <- function(data, series = NULL, title = TRUE, path = ".", 
                        save = FALSE) {
  
  
  if (!is.list(data)) {
    stop("Argument 'data' must be an output of format_data().")
  }
  
  if (!("data_converted" %in% names(data[[1]]))) {
    stop("Argument 'data' must be an output of format_data().")
  }
  
  if (is.null(series)) {
    stop("No series name provided. Please use the argument 'series'.")
  }
  
  if (!is.character(series) && !is.numeric(series)) {
    stop("Argument 'series' must be either a series name or a series index.")
  }
  
  if (length(series) != 1) {
    stop("Argument 'series' must be of length 1.")
  }
    
  if (is.numeric(series)) {
    if (!(series %in% 1:length(data))) {
      stop("Argument 'series' is out of series length.")
    }
  }

  if (is.character(series)) {
    if (!(series %in% names(data))) {
      stop("The series '", series, "' cannot be found in 'data'.")
    }
  }
  
  
  ## Subset series ----
  
  data_copy <- data
  data <- data[[series]]
  
  
  ## Get plot extent ----
  
  counts_range <- c(0, max(c(data$"data_converted"$"upper_ci_conv", 
                             data$"data_original"$"counts_orig")))
  # years_range <- range(data$"years")
  years_range <- range(unlist(lapply(data_copy, function(x) x$years)))
  
  
  ## Export as PNG (if required) ----
  
  if (save) {
    
    if (!dir.exists(path)) {
      stop("The directory '", path, "' does not exist.")
    }
    
    filename <- paste0(names(data_copy[series]), "_data.png")
    grDevices::png(filename = file.path(path, filename), width = 12, 
                   height = 6, pointsize = 14, res = 600, units = "in")
  }
    
  
  ## Graphical parameters ----
  
  opar <- par(no.readonly = TRUE)
  
  par(family = "serif", mgp = c(1.5, 0.1, 0), mar = c(1.5, 3.1, 1.9, 1.1),
      cex.axis = 0.75)
  par(mfrow = c(1, 2))
  
  
  ###
  ### Left plot - Original vs. Converted
  ###
  
  ## Empty plot ----
  
  plot(0, type = "n", xlim = years_range, ylim = counts_range, axes = FALSE,
       xlab = "", ylab = "Original (black) and converted (grey) counts",
       font.lab = 1)
  
  
  ## Background ----
  
  grid(lwd = 0.25, lty = 1)
  axis(1, at = axTicks(1), lwd = 0)
  axis(2, at = axTicks(2), lwd = 0)
  box(lwd = 0.5)
  
  
  ## Deviation between original and converted data ----
  
  for (i in 1:nrow(data$"data_original")) {
    
    lines(x = rep(data$"data_original"[i, "year"], 2),
          y = c(data$"data_original"[i, "counts_orig"], 
                data$"data_converted"[i, "counts_conv"]),
          lty = 4, lwd = 0.75)
  }
  
  
  ## Add data ----
  
  points(x = data$"data_original"$"year", 
         y = data$"data_original"$"counts_orig", pch = 19, cex = 1.2)
  
  points(x = data$"data_converted"$"year", 
         y = data$"data_converted"$"counts_conv", pch = 19, cex = 1.2)
  
  points(x = data$"data_converted"$"year", 
         y = data$"data_converted"$"counts_conv", pch = 19, col = "#aaaaaa")
  
  
  ###
  ### Right plot - Converted with 95% CI
  ###
  
  
  ## Empty plot ----
  
  plot(0, type = "n", xlim = years_range, ylim = counts_range, axes = FALSE,
       xlab = "", ylab = "Converted counts with 95% CI",
       font.lab = 1)
  
  
  ## Background ----
  
  grid(lwd = 0.25, lty = 1)
  axis(1, at = axTicks(1), lwd = 0)
  axis(2, at = axTicks(2), lwd = 0)
  box(lwd = 0.5)

  
  ## Add 95% CI ----
  
  for (i in 1:nrow(data$"data_converted")) {
    
    lines(x = rep(data$"data_converted"[i, "year"], 2),
          y = c(data$"data_converted"[i, "lower_ci_conv"], 
                data$"data_converted"[i, "upper_ci_conv"]))
    
    
    lines(x = c(data$"data_converted"[i, "year"] - 
                  data$"data_converted"[i, "year"] * 0.0002,
                data$"data_converted"[i, "year"] + 
                  data$"data_converted"[i, "year"] * 0.0002),
          y = rep(data$"data_converted"[i, "lower_ci_conv"], 2))
    
    lines(x = c(data$"data_converted"[i, "year"] - 
                  data$"data_converted"[i, "year"] * 0.0002,
                data$"data_converted"[i, "year"] + 
                  data$"data_converted"[i, "year"] * 0.0002),
          y = rep(data$"data_converted"[i, "upper_ci_conv"], 2))
  }
  
  
  ## Add data ----

  points(x = data$"data_converted"$"year", 
         y = data$"data_converted"$"counts_conv", pch = 19, cex = 1.2)
  
  points(x = data$"data_converted"$"year", 
         y = data$"data_converted"$"counts_conv", pch = 19, col = "#aaaaaa")
  
  
  ## Add brand ----
  
  if (title)
    mtext(paste(data$"location", "-", data$"species"), side = 3, line = 0.5, 
          adj = 1, at = par()$usr[2], cex = 1.2, font = 4)
    
  
  ## Restore user graphical parameters ----
  
  if (save)
    dev.off()
  
  par(opar)
  
  invisible(NULL) 
}
