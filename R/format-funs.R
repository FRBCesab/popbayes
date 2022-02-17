#' Format count series
#'
#' @description
#' This function provides an easy way to get count series ready to be analyzed 
#' by the package `popbayes`. It must be used prior to all other functions.
#' 
#' This function formats the count series (passed through the argument 
#' `data`) by selecting and renaming columns, checking columns format and 
#' content, and removing missing data (if `na_rm = TRUE`). It converts the 
#' original data frame into a list of count series that will be analyzed later
#' by the function [fit_trend()] to estimate population trends.
#' 
#' To be usable for the estimation of population trends, counts must be 
#' accompanied by information on precision. The population trend model requires 
#' a 95% confident interval (CI).
#' If estimates are total counts or guesstimates, this function will construct 
#' boundaries of the 95% CI by applying the rules set out in 
#' \url{https://frbcesab.github.io/popbayes/articles/popbayes.html}.
#' If counts were estimated by a sampling method the user needs to specify a 
#' measure of precision. Precision is preferably provided in the form of a 95% 
#' CI by means of two fields: `lower_ci` and `upper_ci`. It may also be given 
#' in the form of a standard deviation (`sd`), a variance (`var`), or a 
#' coefficient of variation (`cv`). If the fields `lower_ci` and `upper_ci` are 
#' both absent (or `NA`), fields `sd`, `var`, and `cv` are examined in this 
#' order. When one is found valid (no missing value), a 95% CI is derived 
#' assuming a normal distribution. 
#' The field `stat_method` must be present in `data` to indicate
#' if counts are **total counts** (`'T'`), **sampling** (`'S'`), or 
#' **guesstimate** (`'G'`).
#' 
#' If a series mixes aerial and ground counts, a field `field_method` must 
#' also be present and must contain either `'A'` (aerial counts), or `'G'` 
#' (ground counts). As all counts must eventually refer to the same field 
#' method for a correct estimation of trend, a conversion will be performed to  
#' homogenize counts. This conversion is based on a **preferred field method**
#' and a **conversion factor** both specific to a species/category. 
#' The preferred field method specifies the conversion direction. The 
#' conversion factor is the multiplicative factor that must be applied to an 
#' aerial count to get an equivalent ground count (note that if the preferred 
#' field method is `'A'`, ground counts will be divided by the conversion 
#' factor to get the equivalent aerial count).
#' 
#' The argument `rmax` represents the maximum change in log population size 
#' between two dates (i.e. the relative rate of increase). It will be used 
#' by [fit_trend()] but must be provided in this function.
#' 
#' These three parameters, named `pref_field_method`, `conversion_A2G`, and 
#' `rmax` can be present in `data` or in a second `data.frame` 
#' (passed through the argument `info`).
#' Alternatively, the package `popbayes` provides their values for some 
#' African large mammals.
#' 
#' **Note:** If the field `field_method` is absent in `data`, counts are 
#' assumed to be obtained with one field method.
#'
#'
#'
#' @param data a `data.frame` with at least five columns: `location`, 
#'   `species`, `date`, `count`, and `stat_method`.
#'   
#'   The `stat_method` field indicates the method used to estimate counts. It 
#'   can contain: `T` (total counts), `G` (guesstimate), and/or `S` (sampling). 
#' 
#'   If individual counts were estimated by **sampling**, additional column(s) 
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
#' @param info (optional) a `data.frame` with species in rows and the following
#'   columns: `species` (species name), `pref_field_method`, 
#'   `conversion_A2G`, and `rmax`. See above section **Description** for 
#'   further information on these fields.
#'   Default is `NULL` (i.e. these information must be present in `data` 
#'   if not available in `popbayes`).
#' 
#' @param location a `character` string. The column name in `data` of the
#'   site. This field is used to distinguish count series from different sites
#'   (if required) and to create an unique series name.
#'   Default is `'location'`.
#'   
#' @param species a `character` string. The column name in `data` (and 
#'   in `info` if provided) of the species. This field is used to distinguish 
#'   count series for different species (if required) and to create an unique 
#'   series name.
#'   Default is `'species'`.
#'   
#' @param date a `character` string. The column name in `data` of the date.
#'   This column `date` must be in a numerical form with possibly a decimal 
#'   part.
#'   Default is `'date'`.
#'   
#' @param count a `character` string. The column name in `data` of the
#'   number of individuals. This column must be numerical.
#'   Default is `'count'`.
#'  
#' @param stat_method a `character` string. The column name in `data` of 
#'   the method used to estimate individuals counts. It can contain `'T'` 
#'   (total counts), `'G'` (guesstimate), and/or `'S'` (sampling). If some 
#'   counts are coded as `'S'`, precision column(s) must also be provided (see 
#'   below).
#'   Default is `'stat_method'`. 
#' 
#' @param lower_ci (optional) a `character` string. The column name in `data`
#'   of the lower boundary of the 95% CI of the estimate (i.e. `count`). If 
#'   provided, the upper boundary of the 95% CI (argument `upper_ci`) must be 
#'   also provided. This argument is only required if some counts have been 
#'   estimated by a sampling method. But user may prefer use other precision 
#'   measures, e.g. standard deviation (argument `sd`), variance (argument 
#'   `var`), or coefficient of variation (argument `cv`). 
#'   Default is `'lower_ci'`.
#'   
#' @param upper_ci (optional) a `character` string. The column name in `data`
#'   of the upper boundary of the 95% CI of the estimate (i.e. `count`). If 
#'   provided, the lower boundary of the 95% CI (argument `lower_ci`) must be 
#'   also provided.
#'   Default is `'upper_ci'`.
#'   
#' @param sd (optional) a `character` string. The column name in `data` of 
#'   the standard deviation of the estimate.
#'   Default is `NULL`.
#'   
#' @param var (optional) a `character` string. The column name in `data` of 
#'   the variance of the estimate.
#'   Default is `NULL`.
#'    
#' @param cv (optional) a `character` string. The column name in `data` of 
#'   the coefficient of variation of the estimate.
#'   Default is `NULL`.
#'   
#' @param field_method (optional) a `character` string. The column name in 
#'   `data` of the field method used to count individuals. Counts can be ground 
#'   counts (coded as `'G'`) or aerial counts (coded as `'A'`). This argument 
#'   is optional if individuals have been counted by the same method. See above 
#'   section **Description** for further information on the count conversion.
#'   Default is `'field_method'`.
#'   
#' @param pref_field_method (optional) a `character` string. The column name
#'   in `data` of the preferred field method of the species. This argument is
#'   only required is `field_method` is not `NULL` (i.e. individuals have been 
#'   counted by different methods). Alternatively, this value can be passed in
#'   `info` (or internally retrieved if the species is listed in the package). 
#'   See above section **Description** for further information on the count
#'   conversion.
#'   Default is `'pref_field_method'`.
#' 
#' @param conversion_A2G (optional) a `character` string. The column name
#'   in `data` of the count conversion factor of the species. This argument is
#'   only required if `field_method` is not `NULL` (i.e. individuals have been 
#'   counted by different methods). Alternatively this value can be passed in
#'   `info` (or internally retrieved if the species is listed in the package).
#'   See above section **Description** for further information on the count
#'   conversion.
#'   Default is `'conversion_A2G'`.
#'   
#' @param rmax (optional) a `character` string. The column name in `data` of 
#'   the species demographic potential (i.e. the relative rate of increase of 
#'   the population). This is the change in log population size between two 
#'   dates and will be used later by [fit_trend()].
#'   Default is `'rmax'`.
#'   
#' @param path a `character` string. The directory to save formatted data. 
#'   This directory must exist and can be an absolute or a relative path.
#'   Default is the current working directory.
#' 
#' @param na_rm a `logical.` If `TRUE`, counts with `NA` values will be 
#'   removed.
#'   Default is `FALSE` (returns an error to inform user if `NA` are detected).
#'
#'
#'
#' @return An n-elements `list` (where `n` is the number of count series). The
#'   name of each element of this list is a combination of location and 
#'   species. Each element of the list is a `list` with the following content:
#'   \itemize{
#'   \item \code{location} a `character` string. The name of the series site.
#'   \item \code{species} a `character` string. The name of the series species.
#'   \item \code{date} a `numerical` vector. The sequence of dates of the 
#'     series.
#'   \item \code{n_dates} an `integer.` The number of unique dates.
#'   \item \code{stat_methods} a `character` vector. The different stat methods
#'     of the series.
#'   \item \code{field_methods} (optional) a `character` vector. The different
#'     field methods of the series.
#'   \item \code{pref_field_method} (optional) a `character` string. The 
#'     preferred field method of the species (`'A'` or `'G'`).
#'   \item \code{conversion_A2G} (optional) a `numeric`. The conversion factor
#'     of the species used to convert counts to its preferred field method.
#'   \item \code{rmax} a `numeric`. The maximum population growth rate of the
#'     species.
#'   \item \code{data_original} a `data.frame`. Original data of the series 
#'     with renamed columns. Some rows may have been deleted 
#'     (if `na_rm = TRUE`).
#'   \item \code{data_converted} a `data.frame`. Data containing computed 
#'     boundaries of the 95% CI (`lower_ci_conv` and `upper_ci_conv`). If 
#'     counts have been obtained by different field methods, contains also 
#'     converted counts (`count_conv`) based on the preferred field method and
#'     conversion factor of the species. This `data.frame` will be used by the
#'     function [fit_trend()] to fit population models.
#'   }
#'   
#'   **Note:** Some original series can be discarded if one of these two 
#'   conditions is met: 1) the series contains only zero counts, and 2) the 
#'   series contains only a few dates (< 4 dates).
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
#' garamba_formatted <- popbayes::format_data(garamba, path = temp_path)
#' 
#' ## Number of count series ----
#' length(garamba_formatted)
#' 
#' ## Retrieve count series names ----
#' popbayes::list_series(path = temp_path)
#' 
#' ## Print content of the first count series ----
#' names(garamba_formatted[[1]])
#' 
#' ## Print original data ----
#' garamba_formatted[[1]]$"data_original"
#' 
#' ## Print converted data ----
#' garamba_formatted[[1]]$"data_converted"

format_data <- function(data, info = NULL, date = "date", count = "count", 
                        location = "location", species = "species", 
                        stat_method = "stat_method", lower_ci = "lower_ci", 
                        upper_ci = "upper_ci", sd = NULL, var = NULL, 
                        cv = NULL, field_method = "field_method", 
                        pref_field_method = "pref_field_method",
                        conversion_A2G = "conversion_A2G", rmax = "rmax", 
                        path = ".", na_rm = FALSE) {
  
  
  ## Check Data dataset ----
  
  if (missing(data)) {
    stop("Argument 'data' is required.")
  }
  
  if (!is.data.frame(data)) {
    stop("Argument 'data' must be a data frame.")
  }
  
  
  ## Check Location field ----
  
  if (!is.character(location) || length(location) != 1) {
    stop("Argument 'location' must be a column name (character string).")
  }
  
  if (!(location %in% colnames(data))) {
    stop("The column '", location, "' (argument location) is absent from ", 
         "'data'. Please check the spelling.")
  }
  
  data[ , location] <- as.character(data[ , location])
  
  if (any(is.na(data[ , location]))) {
    stop("The column '", location, "' cannot contain NA.")
  }
  
  
  ## Check Species field ----
  
  if (!is.character(species) || length(species) != 1) {
    stop("Argument 'species' must be a column name (character string).")
  }
  
  if (!(species %in% colnames(data))) {
    stop("The column '", species, "' (argument species) is absent from ", 
         "'data'. Please check the spelling.")
  }
  
  data[ , species] <- as.character(data[ , species])
  
  if (any(is.na(data[ , species]))) {
    stop("The column '", species, "' cannot contain NA.")
  }
  
  
  ## Check Date field ----
  
  if (!is.character(date) || length(date) != 1) {
    stop("Argument 'date' must be a column name (character string).")
  }
  
  if (!(date %in% colnames(data))) {
    stop("The column '", date, "' (argument date) is absent from 'data'. ", 
         "Please check the spelling.")
  }
  
  if (!is.numeric(data[ , date])) {
    stop("The column '", date, "' must be numeric.")
  }
  
  if (any(is.na(data[ , date]))) {
    stop("The column '", date, "' cannot contain NA.")
  }
  
  
  ## Check Count field ----
  
  if (!is.character(count) || length(count) != 1) {
    stop("Argument 'count' must be a column name (character string).")
  }
  
  if (!(count %in% colnames(data))) {
    stop("The column '", count, "' (argument count) is absent from 'data'. ", 
         "Please check the spelling.")
  }
  
  if (!is.numeric(data[ , count])) {
    stop("The column '", count, "' must be numeric.")
  }
  
  if (length(which(data[ , count] < 0))) {
    stop("The column '", count, "' must be positive (or zero).")
  }
  
  
  ## Check Logical ----
  
  if (!is.logical(na_rm)) {
    stop("Argument 'na_rm' must be a TRUE or FALSE.")
  }
  
  
  ## Detect and delete (or error) NA in counts ----
  
  data <- is_na_counts(data, count, na_rm)
  
  if (nrow(data) == 0) {
    stop("All counts are NA. Please check your data.")
  }
  
  species_list <- sort(unique(data[ , species]))
  
  
  ## Check Stat Method field ----
  
  if (!is.character(stat_method) || length(stat_method) != 1) {
    stop("Argument 'stat_method' must be a column name (character string).")
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
      stop("Argument 'lower_ci' must be a column name (character string).")
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
      stop("Argument 'upper_ci' must be a column name (character string).")
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
      stop("Argument 'sd' must be a column name (character string).")
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
      stop("Argument 'var' must be a column name (character string).")
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
      stop("Argument 'cv' must be a column name (character string).")
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
    
    valid_info_colnames <- c("species", "pref_field_method", "conversion_A2G",
                             "rmax")
    
    valid_info_colnames_msg <- paste0(valid_info_colnames, collapse = "' and '")
    valid_info_colnames_msg <- paste0("'", valid_info_colnames_msg, "'")
    
    
    if (any(!(valid_info_colnames %in% colnames(info)))) {
      stop("Invalid columns in 'info'. ",
           "Required variables are: ", valid_info_colnames_msg, ".")
    }
    
    if (any(is.na(info[ , valid_info_colnames]))) {
      stop("The dataset 'info' cannot contain NA.")
    }
  }
  
  
  ## Check conversion data (if required) ----
  
  if (!is.null(field_method)) {
    
    if (!is.null(info)) {           ## Conversion data in info
      
      if (!all(species_list %in% info$"species")) {
        stop("Some species listed in 'data' are missing from 'info'.")
      }
      
      usethis::ui_done("Conversion data found in 'info'.")
      
      conversion_data <- info[info$"species" %in% species_list, ]
      
    } else {                        ## Conversion data in data
      
      if (!is.null(pref_field_method) && !is.null(conversion_A2G)) {
        
        if (!is.character(pref_field_method) || 
            length(pref_field_method) != 1) {
          stop("Argument 'pref_field_method' must be a column name ", 
               "(character string).")
        }
        
        if (!(pref_field_method %in% colnames(data))) {
          stop("The column '", pref_field_method, "' (argument ",
               "pref_field_method) is absent from 'data'. ",
               "Please check the spelling.")
        }
        
        if (any(is.na(data[ , pref_field_method]))) {
          stop("The column '", pref_field_method, "' cannot contain NA.")
        }
        
        if (!is.character(conversion_A2G) || 
            length(conversion_A2G) != 1) {
          stop("Argument 'conversion_A2G' must be a column name ", 
               "(character string).")
        }
        
        if (!(conversion_A2G %in% colnames(data))) {
          stop("The column '", conversion_A2G, "' (argument ",
               "conversion_A2G) is absent from 'data'. ",
               "Please check the spelling.")
        }
        
        if (any(is.na(data[ , conversion_A2G]))) {
          stop("The column '", conversion_A2G, "' cannot contain NA.")
        }
        
        pref_data <- tapply(data[ , pref_field_method], data[ , species], 
                            function(x) unique(x))
        
        pref_data_unique <- unlist(lapply(pref_data, function(x) length(x)))
        
        if (any(pref_data_unique != 1)) {
          stop("Each species must have an unique preferred field method.")
        }
        
        conv_data <- tapply(data[ , conversion_A2G], data[ , species], 
                            function(x) unique(x))
        
        conv_data_unique <- unlist(lapply(conv_data, function(x) length(x)))
        
        if (any(conv_data_unique != 1)) {
          stop("Each species must have an unique conversion factor.")
        }
        
        pref_data <- data.frame("species"           = names(pref_data), 
                                "pref_field_method" = unlist(pref_data))
        
        conv_data <- data.frame("species"           = names(conv_data), 
                                "conversion_A2G"    = unlist(conv_data))
        
        conversion_data <- merge(pref_data, conv_data, by = "species")
        
        usethis::ui_done("Conversion data found in 'data'.")
        
      } else {                      ## Conversion data in popbayes
        
        if (!any(species_list %in% species_info$"species")) {
          stop("Some species listed in 'data' are not available in popbayes. ", 
               "Please use the argument 'info' or add conversion information ",
               "in 'data'.")
        }
        
        conversion_data <- species_info[species_info$"species" %in% 
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
    
    if (!is.numeric(conversion_data[ , "conversion_A2G"])) {
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
  
  
  ## Get rmax ----
  
  if (!is.null(info)) {           ## rmax data in info
    
    if (!all(species_list %in% info$"species")) {
      stop("Some species listed in 'data' are missing from 'info'.")
    }
    
    rmax_data <- info[info$"species" %in% species_list, ]
    
  } else {                        ## rmax data in data
    
    if (!is.null(rmax)) {
      
      if (!is.character(rmax) || length(rmax) != 1) {
        stop("Argument 'rmax' must be a column name (character string).")
      }
      
      if (!(rmax %in% colnames(data))) {
        stop("The column '", rmax, "' (argument rmax) is absent from 'data'. ",
             "Please check the spelling.")
      }
      
      if (any(is.na(data[ , rmax]))) {
        stop("The column '", rmax, "' cannot contain NA.")
      }
      
      rmax_data <- tapply(data[ , rmax], data[ , species], 
                          function(x) unique(x))
      
      rmax_data_unique <- unlist(lapply(rmax_data, function(x) length(x)))
      
      if (any(rmax_data_unique != 1)) {
        stop("Each species must have an unique rmax.")
      }
      
      rmax_data <- data.frame("species" = names(rmax_data), 
                              "rmax"    = unlist(rmax_data))
      
    } else {                      ## rmax data in popbayes
      
      if (!any(species_list %in% species_info$"species")) {
        stop("Some species listed in 'data' are not available in popbayes. ", 
             "Please use the argument 'info' or add rmax information ",
             "in 'data'.")
      }
      
      rmax_data <- species_info[species_info$"species" %in% species_list, ]
    }
  }
  
  if (!is.numeric(rmax_data[ , "rmax"])) {
    stop("rmax must be a numeric.")
  }
  
  rownames(rmax_data) <- NULL
  
  
  ## Rename columns ----
  
  data_renamed <- data[ , c(location, species, date, stat_method)]
  colnames(data_renamed) <- c("location", "species", "date", "stat_method")
  
  if (!is.null(field_method)) {
    data_renamed <- data.frame(data_renamed, 
                               "field_method" = data[ , field_method])
  }
  
  data_renamed <- data.frame(data_renamed, 
                             "count_orig" = data[ , count])
  
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
  
  for (i in seq_len(nrow(series_infos))) {
    
    id <- series_infos[i, "id"]
    
    sel_rows <- which(data_renamed$"location" == series_infos[i, "location"] &
                      data_renamed$"species"  == series_infos[i, "species"])
    
    data_sub <- zero_counts(data_renamed[sel_rows, ], na_rm)
    
    if (nrow(data_sub) < 4) {
      
      if (!na_rm) {
        
        stop("Count series '", id, "' have not enough data (< 4). Remove ",
             "this count series or use 'na_rm = TRUE'.")
        
      } else {
        
        data_sub <- data.frame()
      }
    }
    
    if (nrow(data_sub)) {
      
      data_sub <- data_sub[order(data_sub$"date", decreasing = FALSE), ]
      rownames(data_sub) <- NULL
      
      if (!is.null(field_method)) {
        
        species_row <- which(conversion_data$"species" == 
                             series_infos[i, "species"])
        
        pref_field_method <- conversion_data[species_row, "pref_field_method"]
        conversion_A2G   <- conversion_data[species_row, "conversion_A2G"]
        
        field_methods <- sort(unique(data_sub[ , "field_method"]))
        
      } else {
        
        pref_field_method <- NULL
        conversion_A2G    <- NULL
        field_methods     <- NULL
      }
      
      rmax <- rmax_data[rmax_data$"species" == series_infos[i, "species"], 
                        "rmax"]
      
      
      data_series[[id]] <- list(
        "location"          = series_infos[i, "location"],
        "species"           = series_infos[i, "species"],
        "dates"             = data_sub[ , "date"],
        "n_dates"           = length(unique(data_sub[ , "date"])),
        "stat_methods"      = sort(unique(data_sub[ , "stat_method"])),
        "field_methods"     = field_methods,
        "pref_field_method" = pref_field_method,
        "conversion_A2G"    = conversion_A2G,
        "rmax"              = rmax,
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
    
    stop("No count series detected. Check your data.")
    
  } else {
    
    usethis::ui_done(paste0("Detecting {usethis::ui_value(length(", 
                            "data_series))} count series"))
  }
  
  data_series
}



#' Remove rows or return error if NA counts detected
#'
#' @param data a `data.frame`
#' @param col a `character` string (column to inspect)
#' @param na_rm a `logical`. If `TRUE`, deletes rows. Otherwise, returns an 
#'   error.
#'
#' @return A `data.frame` (same as `data` with possibly some rows removed).
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
#' - Check for CI boundaries: lower and upper boundaries are required. If not, 
#'   returns an error.
#' - Check CI boundaries values: if lower_ci > count or upper_ci < count,
#'   returns an error.
#'
#' @param data a `data.frame`
#' @param precision_cols a `character` vector (column names of precision 
#'   information, i.e. `'lower_ci_orig'`, `'sd_orig'`, etc.)
#' @param na_rm a `logical`. If `TRUE` delete rows. Otherwise return an error.
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
  
  
  ## Remove or stop if missing precision values for S ----
  
  if (sum(is_na_precision)) {
    
    if (!na_rm) {
      
      stop("Precision column(s) cannot all be NA for sampling counts. If you ", 
           "want to remove counts missing precision information, please use ", 
           "'na_rm = TRUE'.")
      
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
            stop("Unless another type of precision information is provided, ", 
                 "both lower and upper CI bounds are required.")
          }
          
        } else {
          
          stop("Unless another type of precision information is provided, ", 
               "both lower and upper CI bounds are required.")
        }
      }
    }
  }
  
  
  ## Check CI bounds values for S ----
  
  if (length(sampling_rows)) {
    
    if ("lower_ci_orig" %in% colnames(data)) {
      
      pos <- which(data[sampling_rows, "lower_ci_orig"] > 
                     data[sampling_rows, "count_orig"])
      
      if (length(pos)) {
        stop("At least one CI lower bound is greater than the corresponding ", 
             "count.")
      }
      
      
      pos <- which(data[sampling_rows, "lower_ci_orig"] < 0)
      
      if (length(pos)) {
        stop("CI lower bounds must be positive.")
      }
    }
    
    
    if ("upper_ci_orig" %in% colnames(data)) {
      
      pos <- which(data[sampling_rows, "upper_ci_orig"] < 
                     data[sampling_rows, "count_orig"])
      
      if (length(pos)) {
        stop("At least one CI upper bound is smaller than the corresponding ", 
             "count.")
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
#' - Always computes CI bounds for stat methods T and G
#' - Derives CI bounds for stat method S from SD, VAR, or CV, unless CI bounds 
#'   are provided
#'
#' @param data a `data.frame`
#' @param precision_cols a `character` vector (names of columns with precision 
#'   information, i.e. `'lower_ci_orig'`, `'sd_orig'`, etc.)
#'
#' @return A `data.frame` identical to `data` with two additional columns: 
#'   `lower_ci_conv` and `upper_ci_conv`.
#' 
#' @noRd

compute_ci <- function(data, precision_cols) {
  
  data$"count_conv"    <- data$"count_orig"
  data$"lower_ci_conv" <- NA
  data$"upper_ci_conv" <- NA
  
  
  ## Compute CI boundaries for Total counts ----
  
  pos <- which(data$"stat_method" == "T")
  
  if (length(pos)) {
    data[pos, "lower_ci_conv"] <- data[pos, "count_orig"] * 0.95
    data[pos, "upper_ci_conv"] <- data[pos, "count_orig"] * 1.20
  }
  
  
  ## Compute CI boundaries for Guesstimates ----
  
  pos <- which(data$"stat_method" == "G")
  
  if (length(pos)) {
    data[pos, "lower_ci_conv"] <- data[pos, "count_orig"] * 0.80
    data[pos, "upper_ci_conv"] <- data[pos, "count_orig"] * 1.20
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
              data[i, "count_orig"] - 1.96 * data[i, "sd_orig"]
            data[i, "upper_ci_conv"] <- 
              data[i, "count_orig"] + 1.96 * data[i, "sd_orig"]
            found <- 1
          }
        }
      }
      
      if (found == 0) {
        if ("var_orig" %in% precision_cols) {
          if (!is.na(data[i, "var_orig"])) {               # Variance
            data[i, "lower_ci_conv"] <- 
              data[i, "count_orig"] - 1.96 * sqrt(data[i, "var_orig"])
            data[i, "upper_ci_conv"] <- 
              data[i, "count_orig"] + 1.96 * sqrt(data[i, "var_orig"])
            found <- 1
          }
        }
      }
      
      if (found == 0) {
        if ("cv_orig" %in% precision_cols) {
          if (!is.na(data[i, "cv_orig"])) {                # Coeff of variation
            data[i, "lower_ci_conv"] <- 
              data[i, "count_orig"] * (1 - 1.96 * data[i, "cv_orig"])
            data[i, "upper_ci_conv"] <- 
              data[i, "count_orig"] * (1 + 1.96 * data[i, "cv_orig"])
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
#' This function converts counts (and 95% CI lower and upper bounds) based on
#' the preferred field method and the species conversion factor.
#' 
#' **Important:** if the preferred field method is provided and the 
#' `field_method` column is present in `data`, counts are always converted 
#' toward the preferred field method.
#'
#' @param data a `data.frame`. Counts dataset.
#' 
#' @param field_method a `character` string. The column name in `data`.
#' 
#' @param conversion_data a `data.frame`. Conversion data (see `species_info`).
#'
#' @return A `data.frame` (same as `data`).
#' 
#' @noRd

convert_counts <- function(data, field_method, conversion_data) {
  
  if (!is.null(field_method)) {
    
    series_infos <- get_series(data, quiet = TRUE)
    
    for (i in seq_len(nrow(series_infos))) {
      
      species_name <- series_infos[i, "species"]
      
      series_rows <- which(data$"location" == series_infos[i, "location"] & 
                           data$"species" == series_infos[i, "species"])
      
      methods_used <- data[series_rows, "field_method"]
      
      conv_row    <- which(conversion_data$"species" == species_name)
      method_pref <- conversion_data[conv_row, "pref_field_method"]
      conv_fact   <- conversion_data[conv_row, "conversion_A2G"]
      
      conv_fact   <- ifelse(method_pref == "A", 1 / conv_fact, conv_fact)
      conv_fact   <- ifelse(methods_used == method_pref, 1, conv_fact)
      
      data[series_rows, "count_conv"]   <- 
        data[series_rows, "count_conv"]    * conv_fact
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
#' For a count series, identify zero counts.
#' 
#' If there are only zero counts, returns an error (`na_rm = FALSE`) or delete 
#' series (`na_rm = TRUE`). 
#' 
#' If there are zero and non-zero counts, replaces 0 counts by the smaller 
#' non-zero count (and replaces `lower_ci_conv` and `upper_ci_conv` by 
#' the corresponding values).
#'
#' @param data a `data.frame`
#' 
#' @param na_rm a `logical`. If `TRUE` delete series with all 0 counts. 
#'   Otherwise, return an error.
#'
#' @return A `data.frame` (same as `data`).
#' 
#' @noRd

zero_counts <- function(data, na_rm) {
  
  pos <- which(data[ , "count_conv"] == 0)
  
  if (length(pos)) {
    
    if (length(pos) == nrow(data)) {
      
      if (!na_rm) {
        
        stop("Some series have only zero counts. Check your data or ", 
             "'use na_rm = TRUE'.")
        
      } else {
        
        return(data[-pos, ])
      }
      
    } else {
      
      non_zero_counts  <- data[-pos, ]
      which_min_counts <- which.min(non_zero_counts[ , "count_conv"])[1]
      
      data[pos, "count_conv"]    <- non_zero_counts[which_min_counts, 
                                                    "count_conv"]
      data[pos, "lower_ci_conv"] <- non_zero_counts[which_min_counts, 
                                                    "lower_ci_conv"]
      data[pos, "upper_ci_conv"] <- non_zero_counts[which_min_counts, 
                                                    "upper_ci_conv"]
    }
  }
  
  data
}
