# Import a list of count series previously exported

This function imports a list of count series data previously exported by
[`format_data()`](https://frbcesab.github.io/popbayes/reference/format_data.md).
Users can import one, several, or all count series data.

## Usage

``` r
read_series(series = NULL, path = ".")
```

## Arguments

- series:

  a vector of `character` strings. One or several count series names to
  be imported. If `NULL` (default), all available count series will be
  imported.

- path:

  a `character` string. The directory in which count series have been
  saved by the function
  [`format_data()`](https://frbcesab.github.io/popbayes/reference/format_data.md).

## Value

An n-element `list` (where `n` is the number of count series). See
[`format_data()`](https://frbcesab.github.io/popbayes/reference/format_data.md)
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

## Import all count series ----
count_series <- popbayes::read_series(path = temp_path)

## Import one count series ----
a_bus <- popbayes::read_series(series = "garamba__alcelaphus_buselaphus",
                               path   = temp_path)
```
