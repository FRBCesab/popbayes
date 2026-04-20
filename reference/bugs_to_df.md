# Extract estimated parameters from a list of BUGS outputs

From the output of the function
[`fit_trend()`](https://frbcesab.github.io/popbayes/reference/fit_trend.md)
(or
[`read_bugs()`](https://frbcesab.github.io/popbayes/reference/read_bugs.md)),
this function extracts estimated parameters into a `data.frame`.

The resulting `data.frame` has no particular use in `popbayes` but it
can be useful for users.

## Usage

``` r
bugs_to_df(data)
```

## Arguments

- data:

  a named `list` of BUGS outputs. The output of
  [`fit_trend()`](https://frbcesab.github.io/popbayes/reference/fit_trend.md)
  or
  [`read_bugs()`](https://frbcesab.github.io/popbayes/reference/read_bugs.md)

## Value

A `data.frame`.

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
#> module glm loaded
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

## Import BUGS outputs for one count series ----
bugs <- popbayes::read_bugs(series = "garamba__alcelaphus_buselaphus", 
                            path   = temp_path)

## Extract estimated parameters ----
popbayes::bugs_to_df(bugs)
#>                            series parameter          mean           sd
#> 1  garamba__alcelaphus_buselaphus      N[1]  1.790628e+04 1.689316e+03
#> 2  garamba__alcelaphus_buselaphus      N[2]  3.883761e+03 5.036603e+02
#> 3  garamba__alcelaphus_buselaphus      N[3]  3.330386e+03 3.923482e+02
#> 4  garamba__alcelaphus_buselaphus      N[4]  3.010467e+03 4.223448e+02
#> 5  garamba__alcelaphus_buselaphus      N[5]  2.665824e+03 3.322338e+02
#> 6  garamba__alcelaphus_buselaphus      N[6]  3.193770e+03 5.397259e+02
#> 7  garamba__alcelaphus_buselaphus      N[7]  3.713310e+03 6.753730e+02
#> 8  garamba__alcelaphus_buselaphus      N[8]  3.467460e+03 3.706429e+02
#> 9  garamba__alcelaphus_buselaphus      N[9]  2.890285e+03 2.078029e+02
#> 10 garamba__alcelaphus_buselaphus     N[10]  2.752082e+03 2.100026e+02
#> 11 garamba__alcelaphus_buselaphus     N[11]  2.761856e+03 2.631722e+02
#> 12 garamba__alcelaphus_buselaphus     N[12]  2.634945e+03 3.012109e+02
#> 13 garamba__alcelaphus_buselaphus     N[13]  1.310768e+03 7.332821e+01
#> 14 garamba__alcelaphus_buselaphus     N[14]  1.570273e+03 8.972287e+01
#> 15 garamba__alcelaphus_buselaphus     N[15]  2.392812e+03 1.519921e+02
#> 16 garamba__alcelaphus_buselaphus  deviance  2.363015e+02 5.181431e+00
#> 17 garamba__alcelaphus_buselaphus     meanr -4.902770e-02 2.819670e-03
#> 18 garamba__alcelaphus_buselaphus      r[1] -2.188956e-01 2.249388e-02
#> 19 garamba__alcelaphus_buselaphus      r[2] -1.522620e-01 6.868541e-02
#> 20 garamba__alcelaphus_buselaphus      r[3] -5.193814e-02 5.945068e-02
#> 21 garamba__alcelaphus_buselaphus      r[4] -2.391272e-02 3.275661e-02
#> 22 garamba__alcelaphus_buselaphus      r[5]  8.711951e-02 6.419089e-02
#> 23 garamba__alcelaphus_buselaphus      r[6]  7.427235e-02 5.737202e-02
#> 24 garamba__alcelaphus_buselaphus      r[7] -1.924144e-02 5.199198e-02
#> 25 garamba__alcelaphus_buselaphus      r[8] -8.945231e-02 5.030599e-02
#> 26 garamba__alcelaphus_buselaphus      r[9] -2.466457e-02 4.405818e-02
#> 27 garamba__alcelaphus_buselaphus     r[10]  1.906944e-03 6.695877e-02
#> 28 garamba__alcelaphus_buselaphus     r[11] -4.904788e-02 6.548532e-02
#> 29 garamba__alcelaphus_buselaphus     r[12] -8.665564e-02 1.547881e-02
#> 30 garamba__alcelaphus_buselaphus     r[13]  9.028380e-02 3.556965e-02
#> 31 garamba__alcelaphus_buselaphus     r[14]  1.402760e-01 2.739756e-02
#> 32 garamba__alcelaphus_buselaphus       sdr  1.120436e-01 9.450154e-03
#> 33 garamba__alcelaphus_buselaphus    vrrmax  4.614849e-01 3.892327e-02
#>             2.5%           25%           50%           75%         97.5%
#> 1   1.455630e+04  1.679222e+04  1.794258e+04  1.906091e+04  2.108549e+04
#> 2   2.935841e+03  3.535667e+03  3.872586e+03  4.213489e+03  4.910685e+03
#> 3   2.582710e+03  3.062950e+03  3.323238e+03  3.591467e+03  4.112219e+03
#> 4   2.234585e+03  2.716025e+03  2.991303e+03  3.288185e+03  3.882509e+03
#> 5   2.016537e+03  2.443837e+03  2.661995e+03  2.884955e+03  3.337119e+03
#> 6   2.235988e+03  2.808924e+03  3.168554e+03  3.544353e+03  4.326298e+03
#> 7   2.527789e+03  3.229544e+03  3.671713e+03  4.145033e+03  5.160717e+03
#> 8   2.764007e+03  3.213610e+03  3.460817e+03  3.713733e+03  4.213639e+03
#> 9   2.493077e+03  2.747081e+03  2.889112e+03  3.031210e+03  3.300216e+03
#> 10  2.350436e+03  2.608777e+03  2.749921e+03  2.894179e+03  3.168772e+03
#> 11  2.262168e+03  2.580459e+03  2.756181e+03  2.937379e+03  3.293515e+03
#> 12  2.075658e+03  2.422325e+03  2.629233e+03  2.834604e+03  3.248096e+03
#> 13  1.167060e+03  1.261315e+03  1.310487e+03  1.360599e+03  1.456304e+03
#> 14  1.395943e+03  1.510051e+03  1.569787e+03  1.630492e+03  1.747357e+03
#> 15  2.094038e+03  2.291104e+03  2.392646e+03  2.495365e+03  2.688947e+03
#> 16  2.274005e+02  2.326555e+02  2.358766e+02  2.394561e+02  2.477067e+02
#> 17 -5.427471e-02 -5.095713e-02 -4.913123e-02 -4.719951e-02 -4.332071e-02
#> 18 -2.629863e-01 -2.338235e-01 -2.188880e-01 -2.039872e-01 -1.750085e-01
#> 19 -2.851717e-01 -1.986203e-01 -1.521632e-01 -1.062113e-01 -1.635170e-02
#> 20 -1.692003e-01 -9.253924e-02 -5.134203e-02 -1.145511e-02  6.390797e-02
#> 21 -8.831008e-02 -4.594500e-02 -2.345445e-02 -1.398169e-03  3.869870e-02
#> 22 -4.173912e-02  4.479642e-02  8.691001e-02  1.300715e-01  2.155638e-01
#> 23 -3.637369e-02  3.514278e-02  7.402677e-02  1.130141e-01  1.882260e-01
#> 24 -1.189778e-01 -5.503172e-02 -1.982559e-02  1.549653e-02  8.470463e-02
#> 25 -1.867284e-01 -1.233874e-01 -9.041489e-02 -5.631910e-02  1.278464e-02
#> 26 -1.110432e-01 -5.448632e-02 -2.426217e-02  5.643973e-03  5.989999e-02
#> 27 -1.311905e-01 -4.252263e-02  2.182841e-03  4.747425e-02  1.330028e-01
#> 28 -1.784189e-01 -9.310263e-02 -4.833946e-02 -5.093828e-03  7.814803e-02
#> 29 -1.159740e-01 -9.737567e-02 -8.697705e-02 -7.610833e-02 -5.614843e-02
#> 30  2.011630e-02  6.656763e-02  9.023546e-02  1.142945e-01  1.594547e-01
#> 31  8.602679e-02  1.219673e-01  1.406142e-01  1.586507e-01  1.932516e-01
#> 32  9.430415e-02  1.055843e-01  1.117477e-01  1.181311e-01  1.316164e-01
#> 33  3.884197e-01  4.348804e-01  4.602663e-01  4.865581e-01  5.421014e-01
#>        Rhat n.eff
#> 1  1.020905   770
#> 2  1.000976 27000
#> 3  1.001021 25000
#> 4  1.001295  4500
#> 5  1.001089 12000
#> 6  1.002497   980
#> 7  1.002446  1000
#> 8  1.000969 27000
#> 9  1.001274  4800
#> 10 1.001179  6900
#> 11 1.001156  7700
#> 12 1.000983 27000
#> 13 1.001130  8900
#> 14 1.001017 27000
#> 15 1.001020 25000
#> 16 1.003186   680
#> 17 1.011542   890
#> 18 1.001527  2700
#> 19 1.001400  3400
#> 20 1.001833  1700
#> 21 1.001010 27000
#> 22 1.002627   900
#> 23 1.000971 27000
#> 24 1.002824   810
#> 25 1.001208  6100
#> 26 1.000964 27000
#> 27 1.000971 27000
#> 28 1.001108 10000
#> 29 1.001067 14000
#> 30 1.001226  5700
#> 31 1.000965 27000
#> 32 1.001592  2400
#> 33 1.001592  2400
# }
```
