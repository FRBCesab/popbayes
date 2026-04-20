# Import a list of BUGS outputs previously exported

This function imports a list of BUGS outputs previously exported by
[`fit_trend()`](https://frbcesab.github.io/popbayes/reference/fit_trend.md).
Users can import one, several, or all models.

## Usage

``` r
read_bugs(series = NULL, path = ".")
```

## Arguments

- series:

  a vector of `character` strings. One or several count series names. If
  `NULL` (default) BUGS outputs for all count series will be imported.
  Users can run
  [`list_series()`](https://frbcesab.github.io/popbayes/reference/list_series.md)
  to get the correct spelling of count series names.

- path:

  a `character` string. The directory in which BUGS outputs have been
  saved by the function
  [`fit_trend()`](https://frbcesab.github.io/popbayes/reference/fit_trend.md).

## Value

An n-element `list` (where `n` is the number of count series). See
[`fit_trend()`](https://frbcesab.github.io/popbayes/reference/fit_trend.md)
for further information.

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

## Import BUGS outputs for one count series ----
popbayes::read_bugs(series = "garamba__alcelaphus_buselaphus", 
                    path   = temp_path)
#> $garamba__alcelaphus_buselaphus
#> Inference for Bugs model at "model.bug", fit using jags,
#>  2 chains, each with 50000 iterations (first 10000 discarded), n.thin = 3
#>  n.sims = 26666 iterations saved. Running time = 5.319 secs
#>            mu.vect  sd.vect      2.5%       25%       50%       75%     97.5%
#> N[1]     17669.341 1698.331 14355.107 16554.160 17642.370 18790.382 20993.278
#> N[2]      3866.780  500.377  2906.677  3527.571  3858.572  4192.570  4879.499
#> N[3]      3323.201  390.597  2553.193  3059.935  3321.882  3588.000  4089.805
#> N[4]      3016.418  426.742  2231.369  2718.500  3001.103  3296.416  3896.537
#> N[5]      2661.412  323.286  2032.658  2444.082  2660.962  2874.564  3301.868
#> N[6]      3204.082  535.215  2247.256  2825.434  3167.097  3548.097  4350.689
#> N[7]      3740.233  670.450  2578.185  3264.135  3691.336  4159.880  5191.805
#> N[8]      3479.296  366.191  2783.991  3229.807  3473.678  3722.847  4220.713
#> N[9]      2890.884  208.245  2485.579  2750.671  2889.790  3030.026  3302.413
#> N[10]     2749.841  205.812  2355.299  2609.320  2748.213  2886.813  3162.014
#> N[11]     2754.737  256.361  2265.158  2577.409  2751.704  2928.925  3266.973
#> N[12]     2627.303  295.276  2072.194  2420.136  2622.664  2825.832  3217.897
#> N[13]     1311.344   73.735  1168.855  1261.414  1310.842  1361.003  1456.015
#> N[14]     1570.257   89.685  1396.692  1509.286  1569.734  1630.685  1747.368
#> N[15]     2392.817  153.391  2090.754  2290.062  2392.436  2495.648  2692.551
#> meanr       -0.049    0.003    -0.054    -0.051    -0.049    -0.047    -0.043
#> r[1]        -0.218    0.022    -0.263    -0.233    -0.217    -0.202    -0.175
#> r[2]        -0.150    0.068    -0.283    -0.196    -0.150    -0.103    -0.017
#> r[3]        -0.050    0.060    -0.167    -0.090    -0.049    -0.009     0.065
#> r[4]        -0.025    0.034    -0.090    -0.047    -0.025    -0.002     0.042
#> r[5]         0.090    0.063    -0.033     0.046     0.089     0.133     0.216
#> r[6]         0.076    0.057    -0.035     0.038     0.076     0.115     0.191
#> r[7]        -0.021    0.051    -0.121    -0.056    -0.021     0.014     0.081
#> r[8]        -0.091    0.050    -0.189    -0.125    -0.092    -0.058     0.009
#> r[9]        -0.025    0.044    -0.110    -0.055    -0.025     0.005     0.060
#> r[10]        0.000    0.067    -0.133    -0.044     0.001     0.046     0.130
#> r[11]       -0.049    0.066    -0.179    -0.094    -0.049    -0.005     0.079
#> r[12]       -0.086    0.015    -0.115    -0.097    -0.087    -0.076    -0.056
#> r[13]        0.090    0.036     0.020     0.066     0.090     0.114     0.160
#> r[14]        0.140    0.027     0.086     0.122     0.140     0.159     0.193
#> sdr          0.112    0.009     0.094     0.105     0.112     0.118     0.131
#> vrrmax       0.461    0.039     0.388     0.434     0.460     0.486     0.540
#> deviance   236.107    5.246   227.123   232.396   235.631   239.352   247.627
#>           Rhat n.eff
#> N[1]     1.003   650
#> N[2]     1.001 27000
#> N[3]     1.001 27000
#> N[4]     1.001  7600
#> N[5]     1.001  4400
#> N[6]     1.002  2500
#> N[7]     1.001 27000
#> N[8]     1.001  2900
#> N[9]     1.001 10000
#> N[10]    1.001  8000
#> N[11]    1.001 27000
#> N[12]    1.001  7400
#> N[13]    1.001 27000
#> N[14]    1.001  7700
#> N[15]    1.001 27000
#> meanr    1.002  1100
#> r[1]     1.002  2800
#> r[2]     1.001 27000
#> r[3]     1.001  9100
#> r[4]     1.001 27000
#> r[5]     1.001  7300
#> r[6]     1.002  1800
#> r[7]     1.001 14000
#> r[8]     1.001  6400
#> r[9]     1.001  3100
#> r[10]    1.001  4300
#> r[11]    1.001  3400
#> r[12]    1.001  7000
#> r[13]    1.001  8200
#> r[14]    1.001 27000
#> sdr      1.001  6900
#> vrrmax   1.001  6900
#> deviance 1.001  5500
#> 
#> For each parameter, n.eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
#> 
#> DIC info (using the rule: pV = var(deviance)/2)
#> pV = 13.8 and DIC = 249.9
#> DIC is an estimate of expected predictive error (lower deviance is better).
#> 

## Import BUGS outputs for all count series ----
popbayes::read_bugs(path = temp_path)
#> $garamba__alcelaphus_buselaphus
#> Inference for Bugs model at "model.bug", fit using jags,
#>  2 chains, each with 50000 iterations (first 10000 discarded), n.thin = 3
#>  n.sims = 26666 iterations saved. Running time = 5.319 secs
#>            mu.vect  sd.vect      2.5%       25%       50%       75%     97.5%
#> N[1]     17669.341 1698.331 14355.107 16554.160 17642.370 18790.382 20993.278
#> N[2]      3866.780  500.377  2906.677  3527.571  3858.572  4192.570  4879.499
#> N[3]      3323.201  390.597  2553.193  3059.935  3321.882  3588.000  4089.805
#> N[4]      3016.418  426.742  2231.369  2718.500  3001.103  3296.416  3896.537
#> N[5]      2661.412  323.286  2032.658  2444.082  2660.962  2874.564  3301.868
#> N[6]      3204.082  535.215  2247.256  2825.434  3167.097  3548.097  4350.689
#> N[7]      3740.233  670.450  2578.185  3264.135  3691.336  4159.880  5191.805
#> N[8]      3479.296  366.191  2783.991  3229.807  3473.678  3722.847  4220.713
#> N[9]      2890.884  208.245  2485.579  2750.671  2889.790  3030.026  3302.413
#> N[10]     2749.841  205.812  2355.299  2609.320  2748.213  2886.813  3162.014
#> N[11]     2754.737  256.361  2265.158  2577.409  2751.704  2928.925  3266.973
#> N[12]     2627.303  295.276  2072.194  2420.136  2622.664  2825.832  3217.897
#> N[13]     1311.344   73.735  1168.855  1261.414  1310.842  1361.003  1456.015
#> N[14]     1570.257   89.685  1396.692  1509.286  1569.734  1630.685  1747.368
#> N[15]     2392.817  153.391  2090.754  2290.062  2392.436  2495.648  2692.551
#> meanr       -0.049    0.003    -0.054    -0.051    -0.049    -0.047    -0.043
#> r[1]        -0.218    0.022    -0.263    -0.233    -0.217    -0.202    -0.175
#> r[2]        -0.150    0.068    -0.283    -0.196    -0.150    -0.103    -0.017
#> r[3]        -0.050    0.060    -0.167    -0.090    -0.049    -0.009     0.065
#> r[4]        -0.025    0.034    -0.090    -0.047    -0.025    -0.002     0.042
#> r[5]         0.090    0.063    -0.033     0.046     0.089     0.133     0.216
#> r[6]         0.076    0.057    -0.035     0.038     0.076     0.115     0.191
#> r[7]        -0.021    0.051    -0.121    -0.056    -0.021     0.014     0.081
#> r[8]        -0.091    0.050    -0.189    -0.125    -0.092    -0.058     0.009
#> r[9]        -0.025    0.044    -0.110    -0.055    -0.025     0.005     0.060
#> r[10]        0.000    0.067    -0.133    -0.044     0.001     0.046     0.130
#> r[11]       -0.049    0.066    -0.179    -0.094    -0.049    -0.005     0.079
#> r[12]       -0.086    0.015    -0.115    -0.097    -0.087    -0.076    -0.056
#> r[13]        0.090    0.036     0.020     0.066     0.090     0.114     0.160
#> r[14]        0.140    0.027     0.086     0.122     0.140     0.159     0.193
#> sdr          0.112    0.009     0.094     0.105     0.112     0.118     0.131
#> vrrmax       0.461    0.039     0.388     0.434     0.460     0.486     0.540
#> deviance   236.107    5.246   227.123   232.396   235.631   239.352   247.627
#>           Rhat n.eff
#> N[1]     1.003   650
#> N[2]     1.001 27000
#> N[3]     1.001 27000
#> N[4]     1.001  7600
#> N[5]     1.001  4400
#> N[6]     1.002  2500
#> N[7]     1.001 27000
#> N[8]     1.001  2900
#> N[9]     1.001 10000
#> N[10]    1.001  8000
#> N[11]    1.001 27000
#> N[12]    1.001  7400
#> N[13]    1.001 27000
#> N[14]    1.001  7700
#> N[15]    1.001 27000
#> meanr    1.002  1100
#> r[1]     1.002  2800
#> r[2]     1.001 27000
#> r[3]     1.001  9100
#> r[4]     1.001 27000
#> r[5]     1.001  7300
#> r[6]     1.002  1800
#> r[7]     1.001 14000
#> r[8]     1.001  6400
#> r[9]     1.001  3100
#> r[10]    1.001  4300
#> r[11]    1.001  3400
#> r[12]    1.001  7000
#> r[13]    1.001  8200
#> r[14]    1.001 27000
#> sdr      1.001  6900
#> vrrmax   1.001  6900
#> deviance 1.001  5500
#> 
#> For each parameter, n.eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
#> 
#> DIC info (using the rule: pV = var(deviance)/2)
#> pV = 13.8 and DIC = 249.9
#> DIC is an estimate of expected predictive error (lower deviance is better).
#> 
# }
```
