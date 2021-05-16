#' Subset a series based on a location and/or a species
#' 
#' @description
#' From the output of the function [format_data()], this function subsets the
#' list of series based on a species name and/or a site. It both the species and
#' the location are provided, the resulting list will have a length of 1 (data 
#' for the species at the location), otherwise the list will have a length >= 1.
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
#' data(garamba)
#' 
#' ## Format dataset ----
#' garamba_formatted <- format_data(data = garamba)
#' 
#' ## Get series names ----
#' names(garamba_formatted)
#' 
#' ## Get data for Alcelaphus buselaphus (at all sites) ----
#' find_series(data = garamba_formatted, species = "Alcelaphus buselaphus")
#' 
#' ## Get data for at Garamba (for all species) ----
#' find_series(data = garamba_formatted, location = "Garamba")
#' 
#' ## Get data for Alcelaphus buselaphus at Garamba only ----
#' find_series(data = garamba_formatted, species = "Alcelaphus buselaphus", 
#'             location = "Garamba")

find_series <- function(data, species = NULL, location = NULL) {
  
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
