#' @title check count data
#'
#' @description
#' To be usable for the estimation of population trend, count data must be accompanied by dates collected 
#' and information of precision. This function checks that the data have the necessary information and are 
#' in a suitable format for the function fit_trend.
#'
#' @param data [data frame] data to be analysed (contains at the minimum counts, dates and count precision)
#' @param conversion_fac [numeric] multiplicative factor to apply to ground counts n order to obtain equivalent aerial counts
#'
#' @author Nicolas CASAJUS, \email{nicolas.casajus@@fondationbiodiversite.fr}
#' @author Roger PRADEL, \email{roger.pradel@@cefe.cnrs.fr}
#'
#' @export
#'
#' @return
#' ...
#'
#' @examples
#'
#' # No example


check_data <- function(
  data,
  conversion_fac = NULL
) {


  col_names <- c(
    "cinf", "csup", "cv", "sd", "var",
    "stat_method", "field_method",
    "location", "species", "year", "counts"
  )
# remove 'which' in following line?
  if (sum(which(colnames(data) %in% col_names)) < length(col_names)) {

    stop("...")
  }

# check number of species
# check number of location



  # if (sum(which(data$c < 0)) > 0) {
  #
  #   stop("Counts must be positive or zero.")
  #
  # }

  # not na counts

  if (
    sum(which(is.na(data[ , "stat_method"]))) > 0 ||
    sum(which(is.na(data[ , "field_method"]))) > 0
  ) {

    stop("...")
  }

  if (length(unique(data[ , "field_method"])) > 1) {

    if (is.null(conversion_fac)) {

      stop("`conversion_fac` must be provided.")
    }
  }

  # "T", "S", "G"


  for (i in 1:nrow(data)) {

    if (data[i, "counts"] > 1) {

      if (data[i, "stat_method"] == "T") {

        data[i, "cinf"] <- data[i, "counts"] * 0.95
        data[i, "csup"] <- data[i, "counts"] * 1.20
      }

      if (data[i, "stat_method"] == "G") {

        data[i, "cinf"] <- data[i, "counts"] * 0.80
        data[i, "csup"] <- data[i, "counts"] * 1.20
      }

      if (data[i, "stat_method"] == "S") {

        if (!is.na(data[i, "cinf"]) && is.na(data[i, "csup"])) {

          stop("...")
        }

        if (is.na(data[i, "cinf"]) && !is.na(data[i, "csup"])) {

          stop("...")
        }

        if (!is.na(data[i, "cinf"]) && !is.na(data[i, "csup"])) {

          if (data[i, "cinf"] > data[i, "csup"]) {

            stop("...")
          }

          if (data[i, "cinf"] > data[i, "counts"]) {

            stop("...")
          }

          if (data[i, "counts"] > data[i, "csup"]) {

            stop("...")
          }
        }

        if (is.na(data[i, "cinf"]) && is.na(data[i, "csup"])) {

          if (!is.na(data[i, "cv"])) {

            if (data[i, "cv"] <= 0) {

              stop(paste0("`cv` must be strickly positive (observation #", i, ")."))
            }

            data[i, "cinf"] <- data[i, "counts"] * (1 - 1.96 * data[i, "cv"])
            data[i, "csup"] <- data[i, "counts"] * (1 + 1.96 * data[i, "cv"])

          } else {

            if (!is.na(data[i, "sd"])) {

              if (data[i, "sd"] <= 0) {

                stop(paste0("`sd` must be strickly positive (observation #", i, ")."))
              }

              data[i, "cinf"] <- data[i, "counts"] - 1.96 * data[i, "sd"]
              data[i, "csup"] <- data[i, "counts"] + 1.96 * data[i, "sd"]

            } else {

              if (!is.na(data[i, "var"])) {

                if (data[i, "var"] <= 0) {

                  stop(paste0("`var` must be strickly positive (observation #", i, ")."))
                }

                data[i, "cinf"] <- data[i, "counts"] - 1.96 * sqrt(data[i, "var"])
                data[i, "csup"] <- data[i, "counts"] + 1.96 * sqrt(data[i, "var"])

              } else {

                stop(paste0("No precision given for observation #", i, "."))

              } # e_o var
            } # e_o sd
          } # e_o cv
        } # e_o no_ci
      } # e_o field_method
    } # e_o count <= 1
  } # e_o for


  if (length(unique(data[ , "field_method"])) > 1) {

    data[data[ , "field_method"] == "T", "counts"] <- data[data[ , "field_method"] == "T", "counts"] * conversion_fac
    data[data[ , "field_method"] == "T", "cinf"]   <- data[data[ , "field_method"] == "T", "cinf"] * conversion_fac
    data[data[ , "field_method"] == "T", "csup"]   <- data[data[ , "field_method"] == "T", "csup"] * conversion_fac
  }

  data[data[ , "cinf"] < 0, "cinf"] <- 0.001 * data[data[ , "cinf"] < 0, "counts"]

  data[data[ , "counts"] <= 1, "cinf"] <- 0.01
  data[data[ , "counts"] <= 1, "csup"] <- min(data[data[ , "counts"] > 0, "counts"])


  return(data)
}
