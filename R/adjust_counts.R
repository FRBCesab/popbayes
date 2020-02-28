#' @title Prepare 1 count series for use by function fit_trend
#'
#' @description
#' To be usable for the estimation of population trend, count data must be 
#' accompanied by dates collected and information on precision.
#' Two fields are compulsory:
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
#' coefficient may be provided via the parameter 'conversion_factor_A2G'. Alternatively, it is defined internally   
#' in the package.
#'   
#' @param data [data frame] data to be analysed (contains at the minimum counts, dates and count precisions)
#' @param category [string] species category: MLB, LLB, LD, elephant, or Giraffe
#' @param preferred_field_method [character] 'G' for ground or 'A' for aerial.
#' @param conversion_factor_A2G [numeric] multiplicative factor to apply to aerial counts to obtain equivalent ground counts.
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


adjust_counts <- function(
  data,
  category = NULL,
  preferred_field_method = NULL,
  conversion_factor_A2G = NULL,
  ...
) {
  
  if (!is.data.frame(data)) {
    
    stop("data must be a data.frame")
  }
  
  valid_stat_meth <- c(
    "S","T","G"
  )
  valid_field_meth <- c(
    "G","A"
  )
  valid_category <- c(
    "MLB","LLB","LD","Giraffe","Elephant"
  )
  
  n <- nrow(data)
  formatted_data <- data.frame(data, c = data$count, cinf = numeric(n), csup = numeric(n))
  
  #### checking presence of valid counts and dates  #### 
  
  compulsory_col_names <- c(
    "year", "count"
  )

  if (sum(compulsory_col_names %in% colnames(data)) != length(compulsory_col_names)) {
      
    stop("'count' and 'year' fields are required")
  }

  if (!is.numeric(data$count)) {
    
    stop("'count' must be numeric")
  }
  
  if (!is.numeric(data$year)) {
    
    stop("'year' must be numeric")
  }
  
  if (sum(data$count >= 0, na.rm = T) < length(data$count)) {

    stop("Counts must be positive or zero.")
  }

  if (sum(data$year >= 0, na.rm = T) < length(data$count)) {
    
    stop("Dates must be positive or zero.")
  }
  
  #### checking that fields with precision information are valid  #### 
  
  orig_CI_fields <- c(
    "lower_bound", "upper_bound"
  )
  
  isCI <- sum(orig_CI_fields %in% colnames(data)) == length(orig_CI_fields)
  if (isCI) {
    
    if (!(is.numeric(data$lower_bound) && is.numeric(data$upper_bound))) {
      
      stop("'lower_bound' and 'upper_bound' must be numeric")
    }
  }
  
  issd <- 'sd' %in% colnames(data)
  if (issd) {
    
    if (!is.numeric(data$sd)) {
      
      stop("'sd' must be numeric")
    }
  }
  
  isvar <- 'var' %in% colnames(data)
  if (isvar) {
    
    if (!is.numeric(data$var)) {
      
      stop("'var' must be numeric")
    }
  }
  
  iscv <- 'cv' %in% colnames(data)
  if (iscv) {
    
    if (!is.numeric(data$cv)) {
      
      stop("'cv' must be numeric")
    }
  }
  
  noprecision <- !(isCI || issd || isvar || iscv)
  
    if (noprecision) {
    
    msg1 <- "no valid measure of precision is available\n neither CI ('lower_bound', 'upper_bound'), nor 'sd', 'var', or 'cv'"  
  } 
  
  #### checking field stat_method is valid and precision is present for sampling counts ####
  
  isstat_method <- 'stat_method' %in% colnames(data)
      
  if (isstat_method) {
    
    if (!all(data$stat_method %in% valid_stat_meth)) {
      
      stop("invalid statistical method. Allowed stat. methods are ", valid_stat_meth)
    }
  }
  
  if (noprecision && !isstat_method) {
    
    stop(paste(msg1, "statistical method is required (field 'stat_method')", sep = "/n"))
  }
  
  if (noprecision && any(data$stat_method) == 'S') {
    
    stop('measure of precision is required for sampling counts')
  }

  #### building Confidence Intervals #### 

  for (i in 1:n) {
    
    if (isCI &&
        !is.na(data$lower_bound[i]) && 
        !is.na(data$upper_bound[i])
    ) {
      
      if (!(data$lower_bound[i] >= 0 && data$upper_bound[i] >= 0)) {
        
        stop("CI boundaries must be positive or zero (observation #", i,").")
      }
      
      if (data$lower_bound[i] > data$count[i] || 
          data$count[i] > data$upper_bound[i]) {
        
        stop("CI boundaries must verify: lower_bound < count < upper_bound (observation #", i,").")
      }
      
      formatted_data$cinf[i] <- data$lower_bound[i]
      formatted_data$csup[i] <- data$upper_bound[i]
      
    } else {
      
      if (issd &&
          !is.na(data$sd[i])
      ) {
        
        if (data$sd[i] <= 0
        ) {
          
          stop("standard deviation must be strictly positive (observation #", i,").")
        }
        
        formatted_data$cinf[i] <- data$count[i] - 1.96 * data$sd[i]
        formatted_data$csup[i] <- data$count[i] + 1.96 * data$sd[i]
        
      } else {
        
        if (isvar &&
            !is.na(data$var[i])
        ) {
          
          if (data$var[i] <= 0
          ) {
            
            stop("variance must be strictly positive (observation #", i,").")
          }
          
          formatted_data$cinf[i] <- data$count[i] - 1.96 * sqrt(data$sd[i])
          formatted_data$csup[i] <- data$count[i] + 1.96 * sqrt(data$sd[i])

        } else {
          
          if (iscv &&
              !is.na(data$cv[i])
          ) {
            
            if (data$cv[i] <= 0
            ) {
              
#!!!!! it may happen that count = 0 and then cv = 0. What to do in this case?
              
              stop("coefficient of variation must be strictly positive (observation #", i,").")
            }
            
            formatted_data$cinf[i] <- data$count[i] * (1 - 1.96 * data$cv[i])
            formatted_data$csup[i] <- data$count[i] * (1 + 1.96 * data$cv[i])
            
          } else {
            
            if (isstat_method) {
              
              if (data$stat_method[i] == 'T') {
                
                formatted_data$cinf[i] <- formatted_data$count[i] * 0.95
                formatted_data$csup[i] <- formatted_data$count[i] * 1.20
              } else {
                
                if (data$stat_method[i] == 'G') {
                  
                  formatted_data$cinf[i] <- formatted_data$count[i] * 0.80
                  formatted_data$csup[i] <- formatted_data$count[i] * 1.20
                } else {
                  
                  stop("observation #", i,": in the absence of information on precision, statistical method must be 'G' or 'T'")
                }
              }
            } else {
              
              stop("observation #", i,": in the absence of information on precision, the statistical method must be specified")
            } # e_o 
          } # e_o cv              
        } # e_o var      
      } # e_o sd       
    } # e_o CI         
  } # e_o for            
              
  #### reconciling ground and aerial counts ####
  
  # What is available?
  
  isfield_method <- 'field_method' %in% colnames(data)
  ispreferred_field_method <- !is.null(preferred_field_method) && !is.na(preferred_field_method)
  iscategory <- !is.null(category) && !is.na(category)
  isconversion_factor <- !is.null(conversion_factor_A2G) && !is.na(conversion_factor_A2G)
  
  # is this valid information?
  
  if (ispreferred_field_method && !(preferred_field_method %in% valid_field_meth)) {
    
    stop("invalid preferred field method. Allowed field methods are ", valid_field_meth)
  }
  if (iscategory && !(category %in% valid_category)) {
    
    stop("invalid category. Allowed categories are ", valid_category)
  }
  if (isconversion_factor) {
    
    if (!is.numeric(conversion_factor_A2G)) {
      
      stop("conversion factor must be numeric")
    } else {
      
      if (conversion_factor_A2G <= 0) {
        
        stop("conversion factor must be strictly positive")
      }
    }
  }
  if (isfield_method) {
    
    if (!all(data$field_method %in% valid_field_meth)) {
      
      stop("invalid field method. Allowed field methods are ", valid_field_meth)
    } 
  }
  
  # determining preferred field method
    
  if (!ispreferred_field_method) {
    
    if (!iscategory) {
      
      preferred_field_method <- 'A'
      warning(" aerial taken as reference field method by default")  
    } else {
      
      preferred_field_method <- conversion_A2G[ category, "preferred_field_method"]
    }
  }
  
  # determining conversion factor
  
  if (!isconversion_factor) {
    
    if (!iscategory) {
      
      if (all(data$field_method == preferred_field_method)) {
        
        conversion_factor_A2G <- 1
      } else {
        
        stop("at least one of the parameters 'category' and 'conversion_factor_A2G' is required")
      }
    } else {
      
      conversion_factor_A2G <- conversion_A2G[ category, "A2G"]
    }
  }
  
  ## converting towards preferred field method
  
  if (preferred_field_method == 'G') {
    
    conversion_factor <- conversion_factor_A2G
  } else {
    
    conversion_factor <- 1/conversion_factor_A2G
  }
  
  if (!isfield_method) {
    
    warning("As no field method is specified, counts are assumed to be comparable and kept as they are.")
    conversion_factor <- 1
  }
  
  formatted_data[data$field_method != preferred_field_method, c("c","cinf","csup")] <- formatted_data[data$field_method != preferred_field_method, c("c","cinf","csup")] * conversion_factor

  #### special cases ####
  
  ## count is 0
  
  cmin <- min(formatted_data$c[formatted_data$c > 0])
  
  if (cmin == Inf) {
    
    stop("species was never detected. No trend can be derived")
  }
  formatted_data$c[formatted_data$c == 0] <- cmin
  
  ## cinf is negative
  
  formatted_data$cinf[formatted_data$cinf < 0] <- 0
  
  return(formatted_data)

}

  