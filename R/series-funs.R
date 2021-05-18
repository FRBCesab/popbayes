#' Detect counts series and create an unique identifier by series
#'
#' For internal use only.
#'
#' @param data a data frame. 
#' 
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



#' Subset a counts series list based on location and/or species
#' 
#' @description
#' From the output of the function [format_data()], this function subsets the
#' list of counts series based on a species name and/or a site. If both the 
#' species and the location are provided, the resulting list will have a length 
#' of 1 (data for the species at the location), otherwise the list will have a 
#' length >= 1.
#' 
#' @param data a named list. The output of the function [format_data()].
#' 
#' @param species a character of length 1. A species name.
#' 
#' @param location a character of length 1. A site name.
#'
#' @return A subset of `data`, i.e. a named list.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' data("garamba")
#' 
#' ## Format dataset ----
#' garamba_formatted <- popbayes::format_data(garamba)
#' 
#' ## Get series names ----
#' names(garamba_formatted)
#' 
#' ## Or...
#' popbayes::list_series()
#' 
#' ## Get data for Alcelaphus buselaphus (at all sites) ----
#' popbayes::filter_series(garamba_formatted, species = "Alcelaphus buselaphus")
#' 
#' ## Get data at Garamba (for all species) ----
#' popbayes::filter_series(garamba_formatted, location = "Garamba")
#' 
#' ## Get data for Alcelaphus buselaphus at Garamba only ----
#' popbayes::filter_series(garamba_formatted, location = "Garamba",
#'                         species = "Alcelaphus buselaphus")
#' }

filter_series <- function(data, species = NULL, location = NULL) {
  
  if (!is.list(data)) {
    stop("Argument 'data' must be an output of format_data().")
  }
  
  if (!("data_converted" %in% names(data[[1]]))) {
    stop("Argument 'data' must be an output of format_data().")
  }
  
  
  ## No filter methods ----
  
  if (is.null(species) && is.null(location)) {
    
    usethis::ui_oops(paste0("No species nor location provided to filter ", 
                            "series."))
    
    return(NULL)
  }
  
  
  ## Find series by species ----
  
  if (!is.null(species)) {
    
    if (!is.character(species)) {
      stop("Argument 'species' must be a character.")
    }
    
    if (length(species) > 1) {
      stop("Argument 'species' must be a character of length 1.")
    }
    
    species_detected <- unlist(lapply(data, function(x, species) 
      ifelse(x$"species" == species, TRUE, FALSE), species = species))
    
    if (sum(species_detected) == 0) {
      stop("Wrong species spelling.")
    }
  }
  
  
  ## Find series by locations ----
  
  if (!is.null(location)) {
    
    if (!is.character(location)) {
      stop("Argument 'location' must be a character.")
    }
    
    if (length(location) > 1) {
      stop("Argument 'location' must be a character of length 1.")
    }
    
    location_detected <- unlist(lapply(data, function(x, location) 
      ifelse(x$"location" == location, TRUE, FALSE), location = location))
    
    if (sum(location_detected) == 0) {
      stop("Wrong location spelling.")
    }
  }
  
  
  ## Find intersection ----
  
  if (!is.null(species) && !is.null(location)) {
    
    series_match <- which(species_detected & location_detected)
    
    if (length(series_match)) {
      usethis::ui_done(paste0("Found {usethis::ui_value(length(", 
                              "series_match))} series with ", 
                              "{usethis::ui_value(species)} and ", 
                              "{usethis::ui_value(location)}."))
    } else {
      
      usethis::ui_oops(paste0("No series found with ", 
                              "{usethis::ui_value(species)} and ", 
                              "{usethis::ui_value(location)}."))
      
      return(NULL)
    }
  }
  
  
  ## Otherwise ----
  
  if (!is.null(species) && is.null(location)) {
    
    series_match <- species_detected
    
    usethis::ui_done(paste0("Found {usethis::ui_value(",
                            "sum(species_detected))} ",
                            "series with {usethis::ui_value(species)}."))
  }
  
  if (is.null(species) && !is.null(location)) {
    
    series_match <- location_detected
    
    usethis::ui_done(paste0("Found {usethis::ui_value(",
                            "sum(location_detected))} ",
                            "series with {usethis::ui_value(location)}."))
  }
  
  data[series_match]
}



#' Extract original/converted counts series data from a list
#'
#' @description
#' From the output of the function [format_data()] (or [filter_series()]), this 
#' function extracts data frames containing converted counts 
#' (`converted = TRUE`) or original counts (`converted = FALSE`) for one, 
#' several, or all counts series.
#' 
#' The resulting data frame has no particular use in `popbayes` but it can be 
#' useful for users.
#'
#' @param data a named list. The output of [format_data()] or [filter_series()].
#' 
#' @param converted a logical. If `TRUE` (default) extracts converted counts, 
#'   otherwise returns original counts
#'
#' @return A data frame.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' data("garamba")
#' 
#' ## Format dataset ----
#' garamba_formatted <- popbayes::format_data(garamba)
#' 
#' class(garamba_formatted)
#' 
#' ## Extract converted counts data ----
#' converted_data <- popbayes::series_to_df(garamba_formatted, converted = TRUE)
#' 
#' class(converted_data)
#' 
#' ## Extract original counts data ----
#' original_data <- popbayes::series_to_df(garamba_formatted, converted = FALSE)
#' 
#' dim(converted_data)
#' dim(original_data)
#' dim(garamba)
#' }

series_to_df <- function(data, converted = TRUE) {
  
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



#' Import a list of counts series previously exported
#'
#' @description 
#' This function imports a list of counts series data previously exported by 
#' [format_data()]. Users can import one, several, or all counts series data.
#'
#' @param series a vector of characters. One or several counts series names to
#'   be imported. If `NULL` (default) all available counts series will be 
#'   imported.
#'    
#' @param path a character. The directory in which counts series have been 
#'   saved by the function [format_data()].
#'
#' @return A n-elements list (where n is the number of counts series). See
#'   [format_data()] for further information.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' ## Import all counts series ----
#' counts_series <- popbayes::read_series()
#' 
#' ## Import one counts series ----
#' a_buselaphus <- popbayes::read_series("garamba__alcelaphus_buselaphus")
#' }

read_series <- function(series = NULL, path = ".") {
  
  
  if (!dir.exists(path)) {
    stop("The directory '", path, "' does not exist.")
  }
  
  filenames <- list.files(path, recursive = TRUE, pattern = "_data\\.RData")
  
  if (length(filenames) == 0) {
    stop("No counts series can be found.")  
  }
  
  
  ## All available series names ----
  
  series_names <- strsplit(filenames, .Platform$"file.sep")
  series_names <- unlist(lapply(series_names, function(x) x[1]))

  
  if (!is.null(series)) {
    
    if (!is.character(series)) {
      stop("Argument 'series' must be a character (series name(s)).")
    }
    
    if (any(!(series %in% series_names))) {
      stop("Some counts series cannot be found.")
    } 
    
    series_names <- series
  }

  data_series <- list()
  
  for (series in series_names) {
    
    data_series <- c(data_series, get(load(file.path(path, series, 
                                                     paste0(series, 
                                                            "_data.RData")))))
  }
  
  data_series
}



#' Retrieve the counts series names
#'
#' @description This function retrieves the counts series names generated by
#' the function [format_data()].
#'
#' @param path a character. The directory in which counts series have been 
#'   saved by the function [format_data()].
#'
#' @return A vector of counts series names (character).
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' data("garamba")
#' 
#' ## Format dataset ----
#' garamba_formatted <- format_data(garamba)
#' 
#' ## Retrieve counts series names ----
#' popbayes::list_series()
#' }

list_series <- function(path) {
  
  if (!dir.exists(path)) {
    stop("The directory '", path, "' does not exist.")
  }
  
  filenames <- list.files(path, recursive = TRUE, pattern = "_data\\.RData")
  
  if (length(filenames) == 0) {
    stop("No counts series can be found.")  
  }
  
  
  ## All available series names ----
  
  series_names <- strsplit(filenames, .Platform$"file.sep")
  series_names <- unlist(lapply(series_names, function(x) x[1]))
  
  series_names
}
