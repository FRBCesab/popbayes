#' Format individual counts series
#'
#' @description
#' This function provides an easy way to get individual counts series ready to
#' be analyzed by the package `popbayes`. It must be used prior to all other 
#' functions.
#' 
#' This function formats individual counts series (passed through the argument 
#' `data`) by selecting and renaming columns, checking columns format and 
#' content, and removing missing data (if `na_rm = TRUE`). It converts the 
#' original data frame into a list of counts series that will be analyzed later
#' by the function [fit_trend()] to estimate population trend.
#' 
#' To be usable for the estimation of population trend, counts data must be 
#' accompanied by information on precision. The population trend model requires 
#' a 95% confident interval (CI).
#' If estimates are total counts or guesstimates, this function will construct 
#' boundaries of the 95% CI by applying the rules set out in **_???_**.
#' If counts were estimated by a sampling method user needs to specify a 
#' measure of precision. Precision is preferably provided in the form of a 95% 
#' CI by means of two fields: `lower_ci` and `upper_ci`. It may also be given 
#' in the form of a standard deviation (`sd`), a variance (`var`), or a 
#' coefficient of variation (`cv`). If the fields `lower_ci` and `upper_ci` are 
#' both absent (or `NA`), fields `sd`, `var`, and `cv` are examined in this 
#' order. When one is found valid (no missing value), a 95% CI is derived 
#' assuming a normal distribution. 
#' The field `stat_method` must be present in the data frame `data` to indicate
#' if counts are **total counts** (`'T'`), **sampling** (`'S'`), or 
#' **guesstimate** (`'G'`).
#' 
#' If the series mixes aerial and ground counts, a field `field_method` must 
#' also be present and must contain either `'A'` (aerial counts), or `'G'` 
#' (ground counts). As all counts must refer to the same field method, a 
#' conversion will be performed to homogenize counts. This conversion is based 
#' on a **preferred field method** and a **conversion factor** both specific to 
#' a species/category. This conversion factor (a multiplicative factor) 
#' will be apply to an aerial count to get an equivalent ground count (if the
#' preferred field method is `'G'`) or to a ground count to get an equivalent 
#' aerial count (if the preferred field method is `'A'`).
#' These two parameters, named `pref_field_method` and `conversion_fact`, can 
#' be present in the data frame `data` or in the data frame `info`.
#' Alternatively, the package `popbayes` provides their values for some 
#' African large mammals.
#' If the field `field_method` is absent in data, counts are assumed to be 
#' obtained with one field method.
#'
#'
#'
#' @param data a data frame with at least five columns: `location`, `species`, 
#'   `year`, `counts`, and `stat_method`.
#'   
#'   The `stat_method` field indicates the method used to estimate counts. It 
#'   can contain: `T` (total counts), `G` (guesstimate), and/or `S` (sampling). 
#' 
#'   If individuals counts were estimated by **sampling**, additional column(s) 
#'   providing a measure of precision is also required (e.g. `lower_ci` and 
#'   `upper_ci`, or `sd`, `cv`, `var`). Precision metrics can be different 
#'   between counts. For instance, some sampling counts can have a `sd` value 
#'   and others `lower_ci` and `upper_ci`. In that case three columns are 
#'   required (`lower_ci`, `upper_ci`, and `sd`). See above section 
#'   **Description** for further information on the computation of the 95% 
#'   confident interval of estimates. 
#' 
#'   If the individuals were counted by different methods, an additional field 
#'   `field_method` is also required. It can contain: `G` (ground counts) 
#'   and/or `A` (aerial counts). See above section **Description** for further 
#'   information on the counts conversion.
#'   
#'   Others fields can be present either in `data` or `info` (see below).
#' 
#' @param info (optional) a data frame with species in rows and the following
#'   columns: `species` (species name), `pref_field_method`,and 
#'   `conversion_fact`. See above section **Description** for further 
#'   information on these fields.
#'   Default is `NULL` (i.e. these information must be present in `data` 
#'   if not available in `popbayes`).
#' 
#' @param location a character of length 1. The column name in `data` of the
#'   site. This field is used to distinguish counts series from different sites
#'   (if required) and to create an unique series name.
#'   Default is `'location'`.
#'   
#' @param species a character of length 1. The column name in `data` (and 
#'   in `info` if provided) of the species. This field is used to distinguish 
#'   counts series for different species (if required) and to create an unique 
#'   series name.
#'   Default is `'species'`.
#'   
#' @param year a character of length 1. The column name in `data` of the year.
#'   This column `year` must be in the form 1999, 2000, etc. (numeric).
#'   Default is `'year'`.
#'   
#' @param counts a character of length 1. The column name in `data` of the
#'   number of individuals. This column must be numerical.
#'   Default is `'counts'`.
#'  
#' @param stat_method a character of length 1. The column name in `data` of 
#'   the method used to estimate individuals counts. It can contain `'T'` 
#'   (total counts), `'G'` (guesstimate), and/or `'S'` (sampling). If some 
#'   counts are coded as `'S'`, precision column(s) must also be provided (see 
#'   below).
#'   Default is `'stat_method'`. 
#' 
#' @param lower_ci (optional) a character of length 1. The column name in `data`
#'   of the lower boundary of the 95% CI of the estimate (i.e. `counts`). If 
#'   provided the upper boundary of the 95% CI (argument `upper_ci`) must be 
#'   also provided. This argument is only required if some counts have been 
#'   estimated by a sampling method. But user may prefer use other precision 
#'   measures, e.g. standard deviation (argument `sd`), variance (argument 
#'   `var`), or coefficient of variation (argument `cv`). 
#'   Default is `'lower_ci'`.
#'   
#' @param upper_ci (optional) a character of length 1. The column name in `data`
#'   of the upper boundary of the 95% CI of the estimate (i.e. `counts`). If 
#'   provided the lower boundary of the 95% CI (argument `lower_ci`) must be 
#'   also provided.
#'   Default is `'upper_ci'`.
#'   
#' @param sd (optional) a character of length 1. The column name in `data` of 
#'   the standard deviation of the estimate.
#'   Default is `NULL`.
#'   
#' @param var (optional) a character of length 1. The column name in `data` of 
#'   the variance of the estimate.
#'   Default is `NULL`.
#'    
#' @param cv (optional) a character of length 1. The column name in `data` of 
#'   the coefficient of variation of the estimate.
#'   Default is `NULL`.
#'   
#' @param field_method (optional) a character of length 1. The column name in 
#'   `data` of the field method used to count individuals. Counts can be ground 
#'   counts (coded as `'G'`) or aerial counts (coded as `'A'`). This argument 
#'   is optional if individuals have been counted by the same method. See above 
#'   section **Description** for further information on the counts conversion.
#'   Default is `'field_method'`.
#'   
#' @param pref_field_method (optional) a character of length 1. The column name
#'   in `data` of the preferred field method of the species. This argument is
#'   only required is `field_method` is not NULL (i.e. individuals have been 
#'   counted by different methods). Alternatively this value can be passed in
#'   `info`(or internally retrieved if the species is listed in the package). 
#'   See above section **Description** for further information on the counts 
#'   conversion.
#'   Default is `NULL`.
#' 
#' @param conversion_fact (optional) a character of length 1. The column name
#'   in `data` of the counts conversion factor of the species. This argument is
#'   only required is `field_method` is not NULL (i.e. individuals have been 
#'   counted by different methods). Alternatively this value can be passed in
#'   `info` (or internally retrieved if the species is listed in the package).
#'   See above section **Description** for further information on the counts 
#'   conversion.
#'   Default is `NULL`.
#'   
#' @param path a character. The directory to save formatted data. 
#'   This directory must exist and can be an absolute or a relative path.
#' 
#' @param na_rm a logical. If `TRUE` counts and precision measures (for 
#'   sampling counts) with `NA` values will be removed.
#'   Default is `FALSE` (returns an error to inform user if `NA` are detected).
#'
#' @return A n-elements list (where n is the number of counts series). The name 
#'   of each element of this list is a combination of location and species. 
#'   Each element of the list is a list with the following content:
#'   \itemize{
#'   \item \code{location} a character. The name of the series site.
#'   \item \code{species} a character. The name of the series species.
#'   \item \code{years} a numerical vector. The sequence of years of the series.
#'   \item \code{n_years} an integer. The number of unique years.
#'   \item \code{stat_methods} a character vector. The different stat methods 
#'     of the series.
#'   \item \code{field_methods} (optional) a character vector. The different 
#'     field methods of the series.
#'   \item \code{pref_field_method} (optional) a character. The preferred 
#'     field method of the species (`'A'` or `'G'`).
#'   \item \code{conversion_fact} (optional) a numeric. The conversion factor 
#'     of the species used to convert counts from a field method to its 
#'     preferred field method.
#'   \item \code{data_original} a data frame. Original data of the series with 
#'     renamed columns. Some rows can have been deleted (if `na_rm = TRUE`).
#'   \item \code{data_converted} a data frame. Data containing computed 
#'     boundaries of the 95% CI (`lower_ci_conv` and `upper_ci_conv`). If 
#'     counts have been obtained by different field methods, contains also 
#'     converted counts (`counts_conv`) based on the preferred field method and 
#'     conversion factor of the species. This data frame will be used by the 
#'     function [fit_trend()] to fit population models.
#'   }
#'   
#'   **Note:** Some original series can be discarded if one of these two 
#'   conditions is met: 1) the series contains only zero counts, and 2) the 
#'   series contains only a few years (< 4 years).
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
#' length(garamba_formatted)
#' 
#' ## Retrieve counts series names ----
#' popbayes::list_series()
#' }

format_data <- function(data, info = NULL, year = "year", counts = "counts", 
                        location = "location", species = "species", 
                        stat_method = "stat_method", lower_ci = "lower_ci", 
                        upper_ci = "upper_ci", sd = NULL, var = NULL, cv = NULL, 
                        field_method = "field_method", pref_field_method = NULL,
                        conversion_fact = NULL, path = ".", na_rm = FALSE) {
  
  
  ## Check Data dataset ----
  
  if (missing(data)) {
    stop("Argument 'data' is required.")
  }
  
  if (!is.data.frame(data)) {
    stop("Argument 'data' must be a data frame.")
  }
  
  
  ## Check Location field ----
  
  if (!is.character(location) || length(location) != 1) {
    stop("Argument 'location' must be a column name (character of length 1).")
  }
  
  if (!(location %in% colnames(data))) {
    stop("The column '", location, "' (argument location) is absent from ", 
         "'data'. Please check the spelling.")
  }
  
  data[ , location] <- as.character(data[ , location])
  
  if (any(is.na(data[ , location]))) {
    stop("The column '", location, "' cannot contain NA.")
  }
  
  location_list <- sort(unique(data[ , location]))
  
  
  ## Check Species field ----
  
  if (!is.character(species) || length(species) != 1) {
    stop("Argument 'species' must be a column name (character of length 1).")
  }
  
  if (!(species %in% colnames(data))) {
    stop("The column '", species, "' (argument species) is absent from ", 
         "'data'. Please check the spelling.")
  }
  
  data[ , species] <- as.character(data[ , species])
  
  if (any(is.na(data[ , species]))) {
    stop("The column '", species, "' cannot contain NA.")
  }
  
  species_list <- sort(unique(data[ , species]))
  
  
  ## Check Year field ----
  
  if (!is.character(year) || length(year) != 1) {
    stop("Argument 'year' must be a column name (character of length 1).")
  }
  
  if (!(year %in% colnames(data))) {
    stop("The column '", year, "' (argument year) is absent from 'data'. ", 
         "Please check the spelling.")
  }
  
  if (!is.numeric(data[ , year])) {
    stop("The column '", year, "' must be numeric.")
  }
  
  if (any(is.na(data[ , year]))) {
    stop("The column '", year, "' cannot contain NA.")
  }
  
  
  ## Check Counts field ----
  
  if (!is.character(counts) || length(counts) != 1) {
    stop("Argument 'counts' must be a column name (character of length 1).")
  }
  
  if (!(counts %in% colnames(data))) {
    stop("The column '", counts, "' (argument counts) is absent from 'data'. ", 
         "Please check the spelling.")
  }
  
  if (!is.numeric(data[ , counts])) {
    stop("The column '", counts, "' must be numeric.")
  }
  
  if (length(which(data[ , counts] < 0))) {
    stop("The column '", counts, "' must be positive (or zero).")
  }
  
  
  ## Detect and delete (or error) NA in counts ----
  
  data <- is_na_counts(data, counts, na_rm)
  
  if (nrow(data) == 0) {
    stop("All counts are NA. Please check your data.")
  }
  
  
  ## Check Stat Method field ----
  
  if (!is.character(stat_method) || length(stat_method) != 1) {
    stop("Argument 'stat_method' must be a column name (character of length ",
         "1).")
  }
  
  if (!(stat_method %in% colnames(data))) {
    stop("The column '", stat_method, "' (argument stat_method) is absent ", 
         "from 'data'. Please check the spelling.")
  }
  
  valid_stat_methods <- c("T", "G", "S")
  
  valid_stat_methods_msg <- paste0(valid_stat_methods, collapse = "' or '")
  valid_stat_methods_msg <- paste0("'", valid_stat_methods_msg, "'")
  
  data[ , stat_method] <- as.character(data[ , stat_method])
  
  if (any(is.na(data[ , stat_method]))) {
    stop("The column '", stat_method, "' cannot contain NA.")
  }
  
  if (any(!(data[ , stat_method] %in% valid_stat_methods))) {
    stop("Invalid value(s) for 'stat_method' in 'data'. ",
         "Allowed values are: ", valid_stat_methods_msg, ".")
  }
  
  
  ## Check precision columns ----
  
  if (!is.null(lower_ci)) {
    
    if (!is.character(lower_ci) || length(lower_ci) != 1) {
      stop("Argument 'lower_ci' must be a column name (character of length 1).")
    }
    
    if (!(lower_ci %in% colnames(data))) {
      stop("The column '", lower_ci, "' (argument lower_ci) is absent ", 
           "from 'data'. Please check the spelling.")
    }
    
    if (!is.numeric(data[ , lower_ci])) {
      stop("The column '", lower_ci, "' must be numeric.")
    }
    
    if (length(which(data[ , lower_ci] < 0))) {
      stop("The column '", lower_ci, "' must be positive (or zero).")
    }
  }
  
  if (!is.null(upper_ci)) {
    
    if (!is.character(upper_ci) || length(upper_ci) != 1) {
      stop("Argument 'upper_ci' must be a column name (character of length 1).")
    }
    
    if (!(upper_ci %in% colnames(data))) {
      stop("The column '", upper_ci, "' (argument upper_ci) is absent ", 
           "from 'data'. Please check the spelling.")
    }
    
    if (!is.numeric(data[ , upper_ci])) {
      stop("The column '", upper_ci, "' must be numeric.")
    }
    
    if (length(which(data[ , upper_ci] < 0))) {
      stop("The column '", upper_ci, "' must be positive (or zero).")
    }
  }
  
  if (!is.null(sd)) {
    
    if (!is.character(sd) || length(sd) != 1) {
      stop("Argument 'sd' must be a column name (character of length 1).")
    }
    
    if (!(sd %in% colnames(data))) {
      stop("The column '", sd, "' (argument sd) is absent from 'data'. ",
           "Please check the spelling.")
    }
    
    if (!is.numeric(data[ , sd])) {
      stop("The column '", sd, "' must be numeric.")
    }
    
    if (length(which(data[ , sd] <= 0))) {
      stop("The column '", sd, "' must be strictly positive.")
    }
  }
  
  if (!is.null(var)) {
    
    if (!is.character(var) || length(var) != 1) {
      stop("Argument 'var' must be a column name (character of length 1).")
    }
    
    if (!(var %in% colnames(data))) {
      stop("The column '", var, "' (argument var) is absent from 'data'. ",
           "Please check the spelling.")
    }
    
    if (!is.numeric(data[ , var])) {
      stop("The column '", var, "' must be numeric.")
    }
    
    if (length(which(data[ , var] <= 0))) {
      stop("The column '", var, "' must be strictly positive.")
    }
  } 
  
  if (!is.null(cv)) {
    
    if (!is.character(cv) || length(cv) != 1) {
      stop("Argument 'cv' must be a column name (character of length 1).")
    }
    
    if (!(cv %in% colnames(data))) {
      stop("The column '", cv, "' (argument cv) is absent from 'data'. ",
           "Please check the spelling.")
    }
    
    if (!is.numeric(data[ , cv])) {
      stop("The column '", cv, "' must be numeric.")
    }
    
    if (length(which(data[ , cv] <= 0))) {
      stop("The column '", cv, "' must be strictly positive.")
    }
  } 
  
  
  ## Check Precision Information (if required) ----
  
  if ("S" %in% data[ , stat_method]) {
    
    if (is.null(lower_ci) && is.null(upper_ci) && is.null(sd) && is.null(cv) && 
        is.null(var)) {
      
      stop("No valid measure of precision is available for sampling counts. ", 
           "Add 'lower_ci' and 'upper_ci' and/or 'sd', 'var', 'cv' ", 
           "information.")
    }
    
    if (!is.null(lower_ci) && is.null(upper_ci)) {
      stop("You must provide both lower and upper CI.")
    }
    
    if (is.null(lower_ci) && !is.null(upper_ci)) {
      stop("You must provide both lower and upper CI.")
    }
  }
  
  
  ## Check Field Method field (optional) ----
  
  if (!is.null(field_method)) {
    
    if (!is.character(field_method) || length(field_method) != 1) {
      stop("Argument 'field_method' must be a column name (character of ",
           "length 1).")
    }
    
    if (!(field_method %in% colnames(data))) {
      stop("The column '", field_method, "' (argument field_method) is ", 
           "absent from 'data'. Please check the spelling.")
    }
    
    valid_field_methods <- c("G", "A")
    
    valid_field_methods_msg <- paste0(valid_field_methods, collapse = "' or '")
    valid_field_methods_msg <- paste0("'", valid_field_methods_msg, "'")
    
    data[ , field_method] <- as.character(data[ , field_method])
    
    if (any(is.na(data[ , field_method]))) {
      stop("The column '", field_method, "' cannot contain NA.")
    }  
    
    if (any(!(data[ , field_method] %in% valid_field_methods))) {
      stop("Invalid value(s) for 'field_method' in 'data'. ",
           "Allowed values are: ", valid_field_methods_msg, ".")
    }
  }
  
  
  ## Check Info dataset ----
  
  if (!is.null(info)) {
    
    if (!is.data.frame(info)) {
      stop("Argument 'info' must be a data frame.")
    }
    
    valid_info_colnames <- c("species", "pref_field_method", "conversion_fact")
    
    valid_info_colnames_msg <- paste0(valid_info_colnames, collapse = "' and '")
    valid_info_colnames_msg <- paste0("'", valid_info_colnames_msg, "'")
    
    
    if (any(!(colnames(info) %in% valid_info_colnames))) {
      stop("Invalid columns in 'info'. ",
           "Required variables are: ", valid_info_colnames_msg, ".")
    }
    
    if (any(is.na(info))) {
      stop("The dataset 'info' cannot contain NA.")
    }
  }
  
  
  ## Check Logical ----
  
  if (!is.logical(na_rm)) {
    stop("Argument 'na_rm' must be a TRUE or FALSE.")
  }
  
  
  ## Check conversion data (if required) ----
  
  if (!is.null(field_method)) {
    
    if (!is.null(info)) {           ## Conversion data in info
      
      if (!any(species_list %in% info$"species")) {
        stop("Some species listed in 'data' are missing from 'info'.")
      }
      
      usethis::ui_done("Conversion data found in 'info'.")
      
      conversion_data <- info[info$"species" %in% species_list, ]
      
    } else {                        ## Conversion data in data
      
      if (!is.null(pref_field_method) && !is.null(conversion_fact)) {
        
        if (!is.character(pref_field_method) || 
            length(pref_field_method) != 1) {
          stop("Argument 'pref_field_method' must be a column name ", 
               "(character of length 1).")
        }
        
        if (!(pref_field_method %in% colnames(data))) {
          stop("The column '", pref_field_method, "' (argument ",
               "pref_field_method) is absent from 'data'. ",
               "Please check the spelling.")
        }
        
        if (any(is.na(data[ , pref_field_method]))) {
          stop("The column '", pref_field_method, "' cannot contain NA.")
        }
        
        if (!is.character(conversion_fact) || 
            length(conversion_fact) != 1) {
          stop("Argument 'conversion_fact' must be a column name ", 
               "(character of length 1).")
        }
        
        if (!(conversion_fact %in% colnames(data))) {
          stop("The column '", conversion_fact, "' (argument ",
               "conversion_fact) is absent from 'data'. ",
               "Please check the spelling.")
        }
        
        if (any(is.na(data[ , conversion_fact]))) {
          stop("The column '", conversion_fact, "' cannot contain NA.")
        }
        
        pref_data <- tapply(data[ , pref_field_method], data[ , species], 
                            function(x) unique(x))
        
        pref_data_unique <- unlist(lapply(pref_data, function(x) length(x)))
        
        if (any(pref_data_unique != 1)) {
          stop("Each species must have a unique preferred field method.")
        }
        
        conv_data <- tapply(data[ , conversion_fact], data[ , species], 
                            function(x) unique(x))
        
        conv_data_unique <- unlist(lapply(conv_data, function(x) length(x)))
        
        if (any(conv_data_unique != 1)) {
          stop("Each species must have a unique conversion factor.")
        }
        
        pref_data <- data.frame("species"           = names(pref_data), 
                                "pref_field_method" = unlist(pref_data))
        
        conv_data <- data.frame("species"           = names(conv_data), 
                                "conversion_fact"   = unlist(conv_data))
        
        conversion_data <- merge(pref_data, conv_data, by = "species")
        
        usethis::ui_done("Conversion data found in 'data'.")
        
      } else {                      ## Conversion data in popbayes
        
        if (!any(species_list %in% conversion_data$"species")) {
          stop("Some species listed in 'data' are not available in popbayes. ", 
               "Please use the argument 'info' or add conversion information ",
               "in 'data'.")
        }
        
        conversion_data <- conversion_data[conversion_data$"species" %in% 
                                             species_list, ]
        
        usethis::ui_done("Conversion data found in 'popbayes'.")
      }
    }
    
    if (!is.character(conversion_data[ , "pref_field_method"])) {
      stop("Preferred field method must be a character.")
    }
    
    if (any(!(conversion_data[ , "pref_field_method"] %in% 
              valid_field_methods))) {
      stop("Invalid value(s) for 'pref_field_method'. ",
           "Allowed values are: ", valid_field_methods_msg, ".")
    }
    
    if (!is.numeric(conversion_data[ , "conversion_fact"])) {
      stop("Conversion factor must be a numeric.")
    }
    
    rownames(conversion_data) <- NULL
  }
  
  
  ## Detect Precision(s) Columns ----
  
  precision_cols <- NULL
  
  if ("S" %in% data[ , stat_method]) {
    
    if (!is.null(lower_ci)) {
      
      precision_cols <- c(precision_cols, lower_ci)
      names(precision_cols)[length(precision_cols)] <- "lower_ci"
    }
    
    if (!is.null(upper_ci)) {
      
      precision_cols <- c(precision_cols, upper_ci)
      names(precision_cols)[length(precision_cols)] <- "upper_ci"
    }
    
    if (!is.null(sd)) {
      
      precision_cols <- c(precision_cols, sd)
      names(precision_cols)[length(precision_cols)] <- "sd"
    }
    
    if (!is.null(var)) {
      
      precision_cols <- c(precision_cols, var)
      names(precision_cols)[length(precision_cols)] <- "var"
    }
    
    if (!is.null(cv)) {
      
      precision_cols <- c(precision_cols, cv)
      names(precision_cols)[length(precision_cols)] <- "cv"
    }
  }
  
  
  ## Rename columns ----
  
  data_renamed <- data[ , c(location, species, year, stat_method)]
  colnames(data_renamed) <- c("location", "species", "year", "stat_method")
  
  if (!is.null(field_method)) {
    data_renamed <- data.frame(data_renamed, 
                               "field_method" = data[ , field_method])
  }
  
  data_renamed <- data.frame(data_renamed, 
                             "counts_orig" = data[ , counts])
  
  if (!is.null(precision_cols)) {
    
    precision_data <- data[ , precision_cols]
    precision_cols <- paste0(names(precision_cols), "_orig")
    colnames(precision_data) <- precision_cols
    
    data_renamed <- data.frame(data_renamed, precision_data)
  }
  
  
  ## Check precision measures ----
  
  data_renamed <- is_na_precision(data_renamed, precision_cols, na_rm)
  
  if (nrow(data_renamed) == 0) {
    stop("All your counts are sampling counts without precision measures. ", 
         "Please check your data.")
  }
  
  
  ## Check path ----
  
  if (!dir.exists(path)) {
    stop("The directory '", path, "' does not exist.")
  }
  
  
  
  ##
  ## ... end of checks ----
  ## 
  
  
  
  ## Compute 95% CI boundaries ----
  
  data_renamed <- compute_ci(data_renamed, precision_cols)
  
  
  ## Convert counts to preferred field method ----
  
  data_renamed <- convert_counts(data_renamed, field_method, conversion_data)
  
  
  ## Detect series ----
  
  series_infos <- get_series(data_renamed, quiet = TRUE)
  
  
  ## Split original data by series ----
  
  data_series    <- list()
  series_ignored <- 0
  
  for (i in 1:nrow(series_infos)) {
    
    id <- series_infos[i, "id"]
    
    sel_rows <- which(data_renamed$"location" == series_infos[i, "location"] &
                        data_renamed$"species"  == series_infos[i, "species"])
    
    data_sub <- zero_counts(data_renamed[sel_rows, ], na_rm)
    
    if (nrow(data_sub) < 4) {
      
      if (!na_rm) {
        
        stop("Counts series '", id, "' have not enough data (< 4). Remove ",
             "this counts series or use 'na_rm = TRUE'.")
        
      } else {
        
        data_sub <- data.frame()
      }
    }
    
    if (nrow(data_sub)) {
      
      data_sub <- data_sub[order(data_sub$"year", decreasing = FALSE), ]
      rownames(data_sub) <- NULL
      
      if (!is.null(field_method)) {
        
        species_row <- which(conversion_data$"species" == 
                               series_infos[i, "species"])
        
        pref_field_method <- conversion_data[species_row, "pref_field_method"]
        conversion_fact   <- conversion_data[species_row, "conversion_fact"]
        
        field_methods <- sort(unique(data_sub[ , "field_method"]))
        
      } else {
        
        pref_field_method <- NULL
        conversion_fact   <- NULL
        field_methods     <- NULL
      }
      
      
      data_series[[id]] <- list(
        "location"          = series_infos[i, "location"],
        "species"           = series_infos[i, "species"],
        "years"             = data_sub[ , "year"],
        "n_years"           = length(unique(data_sub[ , "year"])),
        "stat_methods"      = sort(unique(data_sub[ , "stat_method"])),
        "field_methods"     = field_methods,
        "pref_field_method" = pref_field_method,
        "conversion_fact"   = conversion_fact,
        "data_original"     = data_sub[ , -grep("_conv", colnames(data_sub))],
        "data_converted"    = data_sub[ , -grep("_orig", colnames(data_sub))])
      
      
      ## Export sub-list ----
      
      species_path <- file.path(path, id)
      dir.create(species_path, showWarnings = FALSE)
      
      formatted_data <- data_series[id]
      save(formatted_data, file = file.path(species_path,
                                            paste0(id, "_data.RData")))
    }
  }
  
  
  if (length(data_series) == 0) {
    
    stop("No series detected. Check your data.")
    
  } else {
    
    usethis::ui_done(paste0("Detecting {usethis::ui_value(length(", 
                            "data_series))} series"))
  }
  
  data_series
}



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



#' Convert counts
#'
#' This function converts counts (and lower and upper bounds 95% CI) based on
#' the preferred field method and the species conversion factor.
#' 
#' **Important:** if the preferred field method is provided and the 
#' `field_method` column is present in `data` counts are always converted to 
#' the preferred field method.
#'
#' @param data a data frame. Counts dataset.
#' 
#' @param field_method a character of length 1. The column name in `data`.
#' 
#' @param conversion_data a data frame. Conversion data (see `conversion_data`).
#'
#' @return A data frame (same as `data`).
#' 
#' @noRd

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
      
      conv_fact   <- ifelse(method_pref == "A", 1 / conv_fact, conv_fact)
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



#' Special cases: zero counts
#'
#' For a counts series, identify zero counts.
#' 
#' If there is only zero counts, returns an error (`na_rm = FALSE`) or delete 
#' series (`na_rm = TRUE`). 
#' 
#' If there is zero counts and non-zero counts, replaces 0 by the minimum 
#' non-zero counts value (and replaces `lower_ci_conv` and `upper_ci_conv` by 
#' the corresponding values).
#'
#' @param data a data frame
#' 
#' @param na_rm a logical. If TRUE delete rows. Otherwise return an error
#'
#' @return A data frame (same as `data`).
#' 
#' @noRd

zero_counts <- function(data, na_rm) {
  
  pos <- which(data[ , "counts_conv"] == 0)
  
  if (length(pos)) {
    
    if (length(pos) == nrow(data)) {
      
      if (!na_rm) {
        
        stop("Some series has only zero counts. Check your data or ", 
             "'use na_rm = TRUE'.")
        
      } else {
        
        return(data[-pos, ])
      }
      
    } else {
      
      non_zero_counts <- data[-pos, ]
      which_min_counts <- which.min(non_zero_counts[ , "counts_conv"])[1]
      
      data[pos, "counts_conv"]   <- non_zero_counts[which_min_counts, 
                                                    "counts_conv"]
      data[pos, "lower_ci_conv"] <- non_zero_counts[which_min_counts, 
                                                    "lower_ci_conv"]
      data[pos, "upper_ci_conv"] <- non_zero_counts[which_min_counts, 
                                                    "upper_ci_conv"]
    }
  }
  
  data
}
