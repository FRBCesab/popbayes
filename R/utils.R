#' Compute rmax from adult female body mass
#' 
#' @description
#' The demographic potential of a species is limited. The intrinsic rate of 
#' increase `rmax` is the maximum increase in log population size that a 
#' species can attain in a year. According to Sinclair (2003), it is related 
#' to the body mass of adult females by the formula:
#' 
#' $$ 1.375*W^{-0.315} $$
#' 
#' @param w a numerical vector. Adult female body mass (in kg).
#' 
#' @export
#'  
#' @return A numerical vector of `rmax` values.
#' 
#' @references 
#' Sinclair (2013) Mammal population regulation, keystone processes and
#'   ecosystem dynamics. _Philosophical Transactions: Biological Sciences_,
#'   **358**, 1729-1740.
#' 
#' @examples
#' body_masses <- c(55, 127)
#' names(body_masses) <- c("Impala", "Tiang")
#' w_to_rmax(body_masses)

w_to_rmax <- function(w) 1.375 * w ^ -0.315
