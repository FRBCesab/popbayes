#' Extract original/converted data from a list
#'
#' @description
#' From the output of the function [format_data()], this function extracts
#' data frames (for all counts series) containing converted counts 
#' (`converted = TRUE`) or original counts (`converted = FALSE`). This function
#' can also be applied on the output of [find_series()].
#' 
#' The resulting data frame has no particular use in `popbayes` but it can be 
#' useful for users.
#'
#' @param data a named list. The output of the function [format_data()].
#' 
#' @param converted a logical. If `TRUE` (default) extracts converted data, 
#'   otherwise returns original data.
#'
#' @return A data frame.
#' 
#' @export
#'
#' @examples
#' data("garamba")
#' 
#' ## Format dataset ----
#' garamba_formatted <- format_data(garamba)
#' 
#' class(garamba_formatted)
#' 
#' ## Extract converted counts data ----
#' converted_data <- get_series_data(garamba_formatted, converted = TRUE)
#' 
#' class(converted_data)
#' 
#' ## Extract original counts data ----
#' original_data <- get_series_data(garamba_formatted, converted = FALSE)
#' 
#' dim(converted_data)
#' dim(original_data)
#' dim(garamba)

get_series_data <- function(data, converted = TRUE) {
  
  if (!is.list(data)) {
    stop("Argument 'data' must be an output of format_data().")
  }
  
  if (!("data_converted" %in% names(data[[1]]))) {
    stop("Argument 'data' must be an output of format_data().")
  }
  
  if (!is.logical(converted)) {
    stop("Argument 'converted' must be TRUE or FALSE.")
  }
  
  if (length(converted) != 1) {
    stop("Argument 'converted' must be TRUE or FALSE.")
  }
  
  element <- ifelse(converted, "data_converted", "data_original")
  
  data <- lapply(data, function(x) x[[element]])
  data <- do.call(rbind.data.frame, data)
  rownames(data) <- NULL
  
  data
}