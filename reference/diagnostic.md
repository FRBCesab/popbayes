# Check if a BUGS model has converged

From the output of the function
[`fit_trend()`](https://frbcesab.github.io/popbayes/reference/fit_trend.md)
(or
[`read_bugs()`](https://frbcesab.github.io/popbayes/reference/read_bugs.md)),
this function checks if the estimation of all parameters of one (or
several) BUGS model has converged. This diagnostic is performed by
comparing the `Rhat` value of each parameter to a `threshold` (default
is `1.1`). If some `Rhat` values are greater than this threshold (no
convergence), a message listing problematic models is displayed.

## Usage

``` r
diagnostic(data, threshold = 1.1)
```

## Arguments

- data:

  a named `list` of BUGS outputs. The output of
  [`fit_trend()`](https://frbcesab.github.io/popbayes/reference/fit_trend.md)
  or
  [`read_bugs()`](https://frbcesab.github.io/popbayes/reference/read_bugs.md).

- threshold:

  a `numeric`.

## Value

No return value.

## Examples

``` r
## Load Garamba raw dataset ----
file_path <- system.file("extdata", "garamba_survey.csv", 
                         package = "popbayes")
                         
garamba <- read.csv(file = file_path)

## Create temporary folder ----
temp_path <- tempdir()

## Format dataset ----
garamba_formatted <- popbayes::format_data(
  data              = garamba, 
  path              = temp_path,
  field_method      = "field_method",
  pref_field_method = "pref_field_method",
  conversion_A2G    = "conversion_A2G",
  rmax              = "rmax")
#> ✔ Detecting 10 count series.
                                        
## Select one serie ----
a_buselaphus <- popbayes::filter_series(garamba_formatted, 
                                        location = "Garamba",
                                        species  = "Alcelaphus buselaphus")
#> ✔ Found 1 series with "Alcelaphus buselaphus" and "Garamba".
# \donttest{
## Fit population trends (requires JAGS) ----
a_buselaphus_mod <- popbayes::fit_trend(a_buselaphus, path = temp_path)
#> Compiling data graph
#>    Resolving undeclared variables
#>    Allocating nodes
#>    Initializing
#>    Reading data back into data table
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 15
#>    Unobserved stochastic nodes: 15
#>    Total graph size: 227
#> 
#> Initializing model
#> 

## Check for convergence ----
popbayes::diagnostic(a_buselaphus_mod)
#> ✔ All models have converged.
# }
```
