#' @title Extract and prepare the individual count series present in a data frame for use by function fit_trend
#'
#' @description
#' To be usable for the estimation of population trend, count data must be accompanied by dates collected and 
#' information on precision. Two fields are compulsory:
#'     - field 'count': contains the original counts
#'     - field 'year': contains the date each count was taken (given as a year with decimal part e.g. 2019.33).
#'     
#' Precision is preferably provided in the form of a 95% confidence interval (CI) by means of 2 fields:
#'     - field 'lower_bound': the lower boundary of the 95% CI
#'     - field 'upper_bound': the upper boundary of the 95% CI
#' It may also be given in the form of a standard deviation, a variance, or a coefficient of variation.  
#' If the fields 'lower_bound' and 'upper_bound' are both absent, or that at least one contains NA, 
#' fields 'sd' (standard deviation), 'var' (variance), and 'cv' (coefficient of variation) are examined 
#' in this order. When one is found valid, a 95% CI is derived assuming a normal distribution. 
#' 
#' If precision is altogether missing, the count must be a total count or a guesstimate and specified as such 
#' in a field 'stat_method'. Then, a 95% CI interval will be constructed by applying the rules set out 
#' in .... 
#'     - field 'stat_method': either T if a total count, S if a sample count, or G if a guesstimate. 
#'         
#' By applying the preceding rules, the function fills in the fields 'cinf' and 'csup' of the output dataframe.  
#' They contain respectively the lower and upper boundaries of the 95% CI. 
#' 
#' If the series mixes aerial and ground counts, a field 'field_method' must be present. 
#'     - field 'field_method': either A for an aerial count, or G for a ground count
#'     
#' As all counts must refer to the same field method, a reference field method must be chosen. It is specified 
#' via the parameter preferred_field_method. If this parameter is not specified, the reference field method 
#' depends on the species category. The species category is specified via the parameter 'category'. 
#'     - 'category' may be MLB for medium-sized light or brown species (20-150kg), 
#'                         LLB for large light or brown species (>150kg), 
#'                         LD for large dark species (>150kg), 
#'                         Giraffe, 
#'                         Elephant
#' By default, the preferred method is aerial for LD, Giraffe and Elephant species, ground for MLB and LLB species.
#' 
#' Conversion to the preferred method uses a category-specific coefficient of conversion, which is the 
#' multiplicative factor to apply to an aerial count to get an equivalent ground count. The conversion 
#' coefficient may be provided via the parameter 'conversion_fac_A2G'. Alternatively, it is defined internally   
#' in the package.
#'   
#' @param data [data frame] data to be analysed (contains at the minimum counts, dates and count precisions)
#'
#' @author Nicolas CASAJUS, \email{nicolas.casajus@@fondationbiodiversite.fr}
#' @author Roger PRADEL, \email{roger.pradel@@cefe.cnrs.fr}
#'
#' @export
#'
#' @return
#' a dataframe with the original data and additional fields required by function fit_trend (trend analysis)  
#'
#' @examples
#'
#' # No example


sort_out_series <- function(data) {
  
  if (!is.data.frame(data)) {
    
    stop("input argument must be a data.frame")
  }
  
  valid_category <- c(
    "MLB","LLB","LD","Giraffe","Elephant"
  )
  valid_field_meth <- c(
    "G","A"
  )
  
  #### identifying series ####
  
  series_discriminants <- c('species',
                          'location'
                          )
  
  SDF <- series_discriminants[series_discriminants %in% colnames(data)]
  series_name <- apply(data[SDF], 1, 'paste0', collapse='')
  identified_series <- unique(series_name)
  
  # template for output #
  
  series_listing <- vector(mode = "list", length = length(identified_series))
  names(series_listing) <- identified_series
  
  #### identifying series characteristics ####
  
  ## which field are present ##
  
  isrmax <- "rmax" %in% colnames(data)
  isCF_A2G <- "conversion_fac_A2G" %in% colnames(data)
  isPFM <- "preferred_field_method" %in% colnames(data)
  iscategory <- "category" %in% colnames(data)
  
  if (isrmax) {
    
    if (!is.numeric(data$rmax)) stop("rmax must be numeric")
  }
  
  if (isCF_A2G) {
    
    if (!is.numeric(data$conversion_fac_A2G)) stop("conversion factor from aerial to ground count must be numeric")
  }
  
  ## default initialization ##
  
  found_rmax <- vector("numeric", length(identified_series))
  names(found_rmax) <- identified_series
  found_rmax[] <- NA
  
  found_CF_A2G <- found_rmax
  
  found_category <- factor(rep(NA, length(identified_series)), valid_category)
  names(found_category) <- identified_series
  
  found_PFM <- factor(rep(NA, length(identified_series)), valid_field_meth)
  names(found_PFM) <- identified_series
  
  ## examining series one by one ##
  
  for (series in identified_series) {
    
    assign(series, data[which(series_name == series), ])
    cur_series <- get(series)
    
    # rmax #
    
    if (isrmax) {
      
      rmax <- unique(cur_series$rmax[!is.na(cur_series$rmax)]) 
      
      if (length(rmax) > 1) stop("series ", series, ": rmax must be unique")
      
      if (length(rmax) != 0) {
        
        if (rmax <= 0) stop("series ", series, ": rmax must be strictly positive")
        
        found_rmax[series] <- rmax
      }
    }
    
    # conversion factor from aerial to ground #
    
    if (isCF_A2G) {
      
      CF_A2G <- unique(cur_series$conversion_fac_A2G[!is.na(cur_series$conversion_fac_A2G)]) 
      
      if (length(CF_A2G) > 1) stop("series ", series, ": conversion factor from aerial to ground count must be unique")
      
      if (length(CF_A2G) != 0) {
        
        if (CF_A2G <= 0) stop("series ", series, ": conversion factor from aerial to ground must be strictly positive")
        
        found_CF_A2G[series] <- CF_A2G
      }
    }
    
    # preferred field method #
    
    if (isPFM) {
      
      PFM <- unique(cur_series$preferred_field_method[!is.na(cur_series$preferred_field_method)]) 
      
      if (length(PFM) > 1) stop("series ", series, ": preferred field method must be unique")
      
      if (length(PFM) != 0) {
        
        if (!PFM %in% valid_field_meth) stop("series ", series, ": invalid preferred field method. Allowed field methods are ", valid_field_meth)
        
        found_PFM[series] <- PFM
      }
    }
    
    # species category #
    
    if (iscategory) {
      
      category <- unique(cur_series$category[!is.na(cur_series$category)]) 
      
      if (length(category) > 1) stop("series ", series, ": species category must be unique")
      
      if (length(category) != 0) {
        
        if (!category %in% valid_category) stop("series ", series, ": invalid species category. Allowed categories are ", valid_category)
        
        found_category[series] <- category
      }
    }
    
    # trying to complete missing information with internal tables #
    
    if ("species" %in% colnames(cur_series)) {
      
      cur_species <- as.character(cur_series$species[[1]])
      
      if (cur_species %in% row.names(species)) {
        
        if ( is.na(found_category[series]) ) found_category[series] <- species[cur_species,"category"]
        
        if ( is.na(found_rmax[series]) ) found_rmax[series] <- species[cur_species,"rmax"]
      }      
    }
    
    if ( !is.na(found_category[series]) ) {
      
      if ( is.na(found_CF_A2G[series]) ) found_CF_A2G[series] <- conversion_A2G[found_category[series],"A2G"]
      
      if ( is.na(found_PFM[series]) ) found_PFM[series] <- conversion_A2G[found_category[series],"preferred_field_method"]
    }
    
    series_listing[[series]] <- list( 
      data = cur_series, 
      rmax = found_rmax[series],
      conversion_factor_A2G = found_CF_A2G[series],
      category = found_category[series],
      preferred_field_method = found_PFM[series]
    )
  }
  
  return(series_listing)
}  
