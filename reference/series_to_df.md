# Extract original/converted count series data from a list

From the output of the function
[`format_data()`](https://frbcesab.github.io/popbayes/reference/format_data.md)
(or
[`filter_series()`](https://frbcesab.github.io/popbayes/reference/filter_series.md)),
this function extracts `data.frame` containing converted counts
(`converted = TRUE`) or original counts (`converted = FALSE`) for one,
several, or all count series.

The resulting `data.frame` has no particular use in `popbayes` but it
can be useful for users.

## Usage

``` r
series_to_df(data, converted = TRUE)
```

## Arguments

- data:

  a named `list`. The output of
  [`format_data()`](https://frbcesab.github.io/popbayes/reference/format_data.md)
  or
  [`filter_series()`](https://frbcesab.github.io/popbayes/reference/filter_series.md).

- converted:

  a `logical`. If `TRUE` (default) extracts converted counts, otherwise
  returns original counts.

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

## Extract converted count data ----
converted_data <- popbayes::series_to_df(garamba_formatted, 
                                         converted = TRUE)

## Extract original count data ----
original_data <- popbayes::series_to_df(garamba_formatted, 
                                        converted = FALSE)

dim(converted_data)
#> [1] 141  12
dim(original_data)
#> [1] 141  11
dim(garamba)
#> [1] 141  11
```
