# Extract the count series corresponding to a location and/or a species

This function identifies the count series relative to a species and/or a
location in a named list like the output of function
[`format_data()`](https://frbcesab.github.io/popbayes/reference/format_data.md).
If both species and location are provided, the series of counts of the
species at the specified location is extracted. Otherwise, all series
corresponding to the specified criterion (species or location) are
extracted.

## Usage

``` r
filter_series(data, species = NULL, location = NULL)
```

## Arguments

- data:

  a named `list`. The output of function
  [`format_data()`](https://frbcesab.github.io/popbayes/reference/format_data.md).

- species:

  a `character` string. A species name.

- location:

  a `character` string. A site name.

## Value

A subset of `data`, i.e. a named `list`.

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

## Number of count series ----
length(garamba_formatted)
#> [1] 10

## Retrieve count series names ----
popbayes::list_series(path = temp_path)
#>  [1] "garamba__alcelaphus_buselaphus"  "garamba__giraffa_camelopardalis"
#>  [3] "garamba__hippotragus_equinus"    "garamba__kobus_ellipsiprymnus"  
#>  [5] "garamba__kobus_kob"              "garamba__loxodonta_africana"    
#>  [7] "garamba__ourebia_ourebi"         "garamba__redunca_redunca"       
#>  [9] "garamba__syncerus_caffer"        "garamba__tragelaphus_scriptus"  

## Get data for Alcelaphus buselaphus (at all sites) ----
x <- popbayes::filter_series(garamba_formatted, 
                             species = "Alcelaphus buselaphus")
#> ✔ Found 1 series with "Alcelaphus buselaphus".

## Get data at Garamba (for all species) ----
x <- popbayes::filter_series(garamba_formatted, 
                             location = "Garamba")
#> ✔ Found 10 series with "Garamba".

## Get data for Alcelaphus buselaphus at Garamba only ----
x <- popbayes::filter_series(garamba_formatted, 
                             location = "Garamba",
                             species  = "Alcelaphus buselaphus")
#> ✔ Found 1 series with "Alcelaphus buselaphus" and "Garamba".
```
