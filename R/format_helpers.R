#' Remove rows or return error if NA counts detected
#'
#' @param data a data frame
#' @param col a character of length 1 (column to inspect)
#' @param na_rm a logical. If TRUE delete rows. Otherwise return an error
#'
#' @return A data frame (same as `data`).
#' 
#' @noRd

is_na_counts <- function(data, col, na_rm) {
  
  if (any(is.na(data[ , col]))) {
    
    if (!na_rm) {
      
      stop("The column '", col, "' cannot contain NA. If you want to ", 
           "remove missing counts, please use 'na_rm = TRUE'.")
      
    } else {
      
      pos <- which(is.na(data[ , col]))
      usethis::ui_info(paste0("Removing {usethis::ui_value(length(pos))} ",
                              "rows with NA values in 'counts' field."))
      
      data <- data[-pos, ]
    }
  }
  
  data
}



#' Check precision data for sampling counts
#'
#' - Check for missing precision values: if neither ci, sd, var, cv are found, 
#'   returns an error (`na_rm = FALSE`) or deletes rows (`na_rm = TRUE`).
#' - Check for CI boundaries: lower and upper boundaries are requires. If not, 
#'   returns an error.
#' - Check CI boundaries values: if lower_ci > estimate or upper_ci < estimate,
#'   returns an error.
#'
#' @param data a data frame
#' @param precision_cols a character vector (column names of precision 
#'   information, i.e. `'lower_ci_orig'`, `'sd_orig'`, etc.)
#' @param na_rm a logical. If `TRUE` delete rows. Otherwise return an error.
#'
#' @return A data frame (same as `data`).
#' 
#' @noRd

is_na_precision <- function(data, precision_cols, na_rm) {
  
  
  ## Check for missing precision values ----
  
  sampling_rows <- which(data[ , "stat_method"] == "S")
  
  if (length(sampling_rows)) {
    
    is_na_precision <- apply(data[sampling_rows, precision_cols], 1, 
                             function(x) {
                               x <- sum(ifelse(is.na(x), 0, 1))
                               ifelse(x == 0, TRUE, FALSE)
                             })
  }
  
  
  ## Remove or stop is missing precision values for S ----
  
  if (sum(is_na_precision)) {
    
    if (!na_rm) {
      
      stop("Precision column(s) cannot contain NA for sampling counts. If you ", 
           "want to remove missing values please use 'na_rm = TRUE'.")
      
    } else {
      
      usethis::ui_info(paste0("Removing {usethis::ui_value(", 
                              "sum(is_na_precision))} rows with NA values ", 
                              "in 'precision' field(s) for (S)amplings."))
      
      data <- data[-sampling_rows[which(is_na_precision)], ]
    }
  }
  
  
  ## Check for CI bounds (require both) for S ----
  
  sampling_rows <- which(data[ , "stat_method"] == "S")
  
  if (length(sampling_rows)) {
    
    if ("lower_ci_orig" %in% colnames(data)) {
      
      is_na_lower <- is.na(data[sampling_rows, "lower_ci_orig"])
      is_na_upper <- is.na(data[sampling_rows, "upper_ci_orig"])
      
      pos <- which((is_na_upper + is_na_lower) == 1)
      
      if (length(pos)) {
        
        ci_cols <- which(precision_cols %in% c("lower_ci_orig", 
                                               "upper_ci_orig"))
        
        tmp <- as.data.frame(data[sampling_rows[pos], precision_cols[-ci_cols]])
        
        if (ncol(tmp)) {
          
          is_na_precision <- apply(tmp, 1, 
                                   function(x) {
                                     x <- sum(ifelse(is.na(x), 0, 1))
                                     ifelse(x == 0, TRUE, FALSE)
                                   })
          
          if (sum(is_na_precision)) {
            stop("Both lower and upper CI are required if others precision ", 
                 "information are not provided.")
          }
          
        } else {
          
          stop("Both lower and upper CI are required if others precision ", 
               "information are not provided.")
        }
      }
    }
  }
  
  
  ## Check CI bounds values for S ----
  
  if (length(sampling_rows)) {
    
    if ("lower_ci_orig" %in% colnames(data)) {
      
      pos <- which(data[sampling_rows, "lower_ci_orig"] > 
                   data[sampling_rows, "counts_orig"])
      
      if (length(pos)) {
        stop("Some lower CI values are greater than estimation.")
      }
      
      
      pos <- which(data[sampling_rows, "lower_ci_orig"] < 0)
      
      if (length(pos)) {
        stop("Lower CI values must be positive.")
      }
    }
    
    
    if ("upper_ci_orig" %in% colnames(data)) {
      
      pos <- which(data[sampling_rows, "upper_ci_orig"] < 
                   data[sampling_rows, "counts_orig"])
      
      if (length(pos)) {
        stop("Some upper CI values are lesser than estimation.")
      }
      
      
      pos <- which(data[sampling_rows, "upper_ci_orig"] < 0)
      
      if (length(pos)) {
        stop("Upper CI values must be positive.")
      }
    }
  }
  
  data
}



#' Compute 95% confident interval boundaries
#'
#' - Always computes CI boundaries for T and G
#' - Derives CI boundaries for S from SD, VAR, or CV, unless CI boundaries are
#'   provided
#'
#' @param data a data frame
#' @param precision_cols a character vector (column names of precision 
#'   information, i.e. `'lower_ci_orig'`, `'sd_orig'`, etc.)
#'
#' @return A data frame identical to `data` with two additional columns: 
#'   `lower_ci_conv` and `upper_ci_conv`.
#' 
#' @noRd

compute_ci <- function(data, precision_cols) {
  
  data$"counts_conv"   <- data$"counts_orig"
  data$"lower_ci_conv" <- NA
  data$"upper_ci_conv" <- NA
  
  
  ## Compute CI boundaries for Total counts ----
  
  pos <- which(data$"stat_method" == "T")
  
  if (length(pos)) {
    data[pos, "lower_ci_conv"] <- data[pos, "counts_orig"] * 0.95
    data[pos, "upper_ci_conv"] <- data[pos, "counts_orig"] * 1.20
  }
  
  
  ## Compute CI boundaries for Guesstimates ----
  
  pos <- which(data$"stat_method" == "G")
  
  if (length(pos)) {
    data[pos, "lower_ci_conv"] <- data[pos, "counts_orig"] * 0.80
    data[pos, "upper_ci_conv"] <- data[pos, "counts_orig"] * 1.20
  }
  
  
  ## Compute CI boundaries for Sampling counts ----
  
  pos <- which(data$"stat_method" == "S")
  
  if (length(pos)) {
    
    for (i in pos) {
      
      found <- 0
      
      if (found == 0) {
        if ("lower_ci_orig" %in% precision_cols) {
          if (!is.na(data[i, "lower_ci_orig"])) {          # CI bounds
            data[i, "lower_ci_conv"] <- data[i, "lower_ci_orig"]
            data[i, "upper_ci_conv"] <- data[i, "upper_ci_orig"]
            found <- 1
          }
        }
      }
      
      if (found == 0) {
        if ("sd_orig" %in% precision_cols) {
          if (!is.na(data[i, "sd_orig"])) {                # Standard deviation
            data[i, "lower_ci_conv"] <- 
              data[i, "counts_orig"] - 1.96 * data[i, "sd_orig"]
            data[i, "upper_ci_conv"] <- 
              data[i, "counts_orig"] + 1.96 * data[i, "sd_orig"]
            found <- 1
          }
        }
      }
      
      if (found == 0) {
        if ("var_orig" %in% precision_cols) {
          if (!is.na(data[i, "var_orig"])) {               # Variance
            data[i, "lower_ci_conv"] <- 
              data[i, "counts_orig"] - 1.96 * sqrt(data[i, "var_orig"])
            data[i, "upper_ci_conv"] <- 
              data[i, "counts_orig"] + 1.96 * sqrt(data[i, "var_orig"])
            found <- 1
          }
        }
      }
      
      if (found == 0) {
        if ("cv_orig" %in% precision_cols) {
          if (!is.na(data[i, "cv_orig"])) {                # Coeff of variation
            data[i, "lower_ci_conv"] <- 
              data[i, "counts_orig"] * (1 - 1.96 * data[i, "cv_orig"])
            data[i, "upper_ci_conv"] <- 
              data[i, "counts_orig"] * (1 + 1.96 * data[i, "cv_orig"])
            found <- 1
          }
        }
      }
    }
  }
  
  data
}



#' Detect unique counts series and create an unique identifier
#'
#' @param data a data frame
#' @param quiet a logical. If `TRUE`, suppress message
#'
#' @return A data frame with three columns:
#'   - `id` (an unique identifier of the series);
#'   - `location` (the site of the series);
#'   - `species` (the species of the series).
#' 
#' @noRd

get_series <- function(data, quiet = TRUE) {
  
  series_id <- unique(paste(data[ , "location"], data[ , "species"], 
                            sep = "__"))
  
  series_infos <- data.frame("id" = series_id)
  
  series_infos$"location" <- unlist(lapply(strsplit(series_id, "__"), 
                                           function(x) x[1]))
  
  series_infos$"species"  <- unlist(lapply(strsplit(series_id, "__"), 
                                           function(x) x[2]))
  
  series_infos <- series_infos[order(series_infos$"id", decreasing = FALSE), ]
  
  series_infos$"id" <- tolower(series_infos$"id")
  series_infos$"id" <- gsub("\\s{2,}", "", series_infos$"id")
  series_infos$"id" <- gsub("\\s", "_", series_infos$"id")
  series_infos$"id" <- gsub("[[:punct:]]", "_", series_infos$"id")
  
  if (!quiet)
    usethis::ui_done(paste0("Detecting {usethis::ui_value(nrow(series_infos))}",
                            " series (with {usethis::ui_value(length(unique(",
                            "series_infos$location)))} locations and ",
                            "{usethis::ui_value(length(unique(",
                            "series_infos$species)))} species)."))
  
  series_infos
}



convert_counts <- function(data, field_method, conversion_data) {
  
  if (!is.null(field_method)) {
    
    series_infos <- get_series(data, quiet = TRUE)
    
    for (i in 1:nrow(series_infos)) {
      
      species_name <- series_infos[i, "species"]
      
      series_rows <- which(data$"location" == series_infos[i, "location"] & 
                           data$"species" == series_infos[i, "species"])
      
      methods_used <- data[series_rows, "field_method"]
      
      conv_row    <- which(conversion_data$"species" == species_name)
      method_pref <- conversion_data[conv_row, "pref_field_method"]
      conv_fact   <- conversion_data[conv_row, "conversion_fact"]
      
      conv_fact   <- ifelse(method_pref == "G", 1 / conv_fact, conv_fact)
      conv_fact   <- ifelse(methods_used == method_pref, 1, conv_fact)
      
      data[series_rows, "counts_conv"]   <- 
        data[series_rows, "counts_conv"]   * conv_fact
      data[series_rows, "lower_ci_conv"] <- 
        data[series_rows, "lower_ci_conv"] * conv_fact
      data[series_rows, "upper_ci_conv"] <- 
        data[series_rows, "upper_ci_conv"] * conv_fact
      
      data[series_rows, "field_method_conv"] <- method_pref
    }
  }
  
  data
}

zero_counts <- function() {}
