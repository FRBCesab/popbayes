# Format count series

This function provides an easy way to get count series ready to be
analyzed by the package `popbayes`. It must be used prior to all other
functions.

This function formats the count series (passed through the argument
`data`) by selecting and renaming columns, checking columns format and
content, and removing missing data (if `na_rm = TRUE`). It converts the
original data frame into a list of count series that will be analyzed
later by the function
[`fit_trend()`](https://frbcesab.github.io/popbayes/reference/fit_trend.md)
to estimate population trends.

To be usable for the estimation of population trends, counts must be
accompanied by information on precision. The population trend model
requires a 95% confident interval (CI). If estimates are total counts or
guesstimates, this function will construct boundaries of the 95% CI by
applying the rules set out in
<https://frbcesab.github.io/popbayes/articles/popbayes.html>. If counts
were estimated by a sampling method the user needs to specify a measure
of precision. Precision is preferably provided in the form of a 95% CI
by means of two fields: `lower_ci` and `upper_ci`. It may also be given
in the form of a standard deviation (`sd`), a variance (`var`), or a
coefficient of variation (`cv`). If the fields `lower_ci` and `upper_ci`
are both absent (or `NA`), fields `sd`, `var`, and `cv` are examined in
this order. When one is found valid (no missing value), a 95% CI is
derived assuming a normal distribution. The field `stat_method` must be
present in `data` to indicate if counts are **total counts** (`'T'`),
**sampling** (`'S'`), or **guesstimate** (`'X'`).

If a series mixes aerial and ground counts, a field `field_method` must
also be present and must contain either `'A'` (aerial counts), or `'G'`
(ground counts). As all counts must eventually refer to the same field
method for a correct estimation of trend, a conversion will be performed
to homogenize counts. This conversion is based on a **preferred field
method** and a **conversion factor** both specific to a
species/category. The preferred field method specifies the conversion
direction. The conversion factor is the multiplicative factor that must
be applied to an aerial count to get an equivalent ground count (note
that if the preferred field method is `'A'`, ground counts will be
divided by the conversion factor to get the equivalent aerial count).

The argument `rmax` represents the maximum change in log population size
between two dates (i.e. the relative rate of increase). It will be used
by
[`fit_trend()`](https://frbcesab.github.io/popbayes/reference/fit_trend.md)
but must be provided in this function.

These three parameters, named `pref_field_method`, `conversion_A2G`, and
`rmax` can be present in `data` or in a second `data.frame` (passed
through the argument `info`). Alternatively, the package `popbayes`
provides their values for some African large mammals.

**Note:** If the field `field_method` is absent in `data`, counts are
assumed to be obtained with one field method.

## Usage

``` r
format_data(
  data,
  info = NULL,
  date = "date",
  count = "count",
  location = "location",
  species = "species",
  stat_method = "stat_method",
  lower_ci = "lower_ci",
  upper_ci = "upper_ci",
  sd = NULL,
  var = NULL,
  cv = NULL,
  field_method = NULL,
  pref_field_method = NULL,
  conversion_A2G = NULL,
  rmax = NULL,
  path = ".",
  na_rm = FALSE
)
```

## Arguments

- data:

  a `data.frame` with at least five columns: `location`, `species`,
  `date`, `count`, and `stat_method`.

  The `stat_method` field indicates the method used to estimate counts.
  It can contain: `T` (total counts), `X` (guesstimate), and/or `S`
  (sampling).

  If individual counts were estimated by **sampling**, additional
  column(s) providing a measure of precision is also required (e.g.
  `lower_ci` and `upper_ci`, or `sd`, `cv`, `var`). Precision metrics
  can be different between counts. For instance, some sampling counts
  can have a `sd` value and others `lower_ci` and `upper_ci`. In that
  case three columns are required (`lower_ci`, `upper_ci`, and `sd`).
  See above section **Description** for further information on the
  computation of the 95% confident interval of estimates.

  If the individuals were counted by different methods, an additional
  field `field_method` is also required. It can contain: `G` (ground
  counts) and/or `A` (aerial counts). See above section **Description**
  for further information on the counts conversion.

  Others fields can be present either in `data` or `info` (see below).

- info:

  (optional) a `data.frame` with species in rows and the following
  columns: `species` (species name), `pref_field_method`,
  `conversion_A2G`, and `rmax`. See above section **Description** for
  further information on these fields. Default is `NULL` (i.e. these
  information must be present in `data` if not available in `popbayes`).

- date:

  a `character` string. The column name in `data` of the date. This
  column `date` must be in a numerical form with possibly a decimal
  part. Default is `'date'`.

- count:

  a `character` string. The column name in `data` of the number of
  individuals. This column must be numerical. Default is `'count'`.

- location:

  a `character` string. The column name in `data` of the site. This
  field is used to distinguish count series from different sites (if
  required) and to create an unique series name. Default is
  `'location'`.

- species:

  a `character` string. The column name in `data` (and in `info` if
  provided) of the species. This field is used to distinguish count
  series for different species (if required) and to create an unique
  series name. Default is `'species'`.

- stat_method:

  a `character` string. The column name in `data` of the method used to
  estimate individuals counts. It can contain `'T'` (total counts),
  `'X'` (guesstimate), and/or `'S'` (sampling). If some counts are coded
  as `'S'`, precision column(s) must also be provided (see below).
  Default is `'stat_method'`.

- lower_ci:

  (optional) a `character` string. The column name in `data` of the
  lower boundary of the 95% CI of the estimate (i.e. `count`). If
  provided, the upper boundary of the 95% CI (argument `upper_ci`) must
  be also provided. This argument is only required if some counts have
  been estimated by a sampling method. But user may prefer use other
  precision measures, e.g. standard deviation (argument `sd`), variance
  (argument `var`), or coefficient of variation (argument `cv`). Default
  is `'lower_ci'`.

- upper_ci:

  (optional) a `character` string. The column name in `data` of the
  upper boundary of the 95% CI of the estimate (i.e. `count`). If
  provided, the lower boundary of the 95% CI (argument `lower_ci`) must
  be also provided. Default is `'upper_ci'`.

- sd:

  (optional) a `character` string. The column name in `data` of the
  standard deviation of the estimate. Default is `NULL`.

- var:

  (optional) a `character` string. The column name in `data` of the
  variance of the estimate. Default is `NULL`.

- cv:

  (optional) a `character` string. The column name in `data` of the
  coefficient of variation of the estimate. Default is `NULL`.

- field_method:

  (optional) a `character` string. The column name in `data` of the
  field method used to count individuals. Counts can be ground counts
  (coded as `'G'`) or aerial counts (coded as `'A'`). This argument is
  optional if individuals have been counted by the same method. See
  above section **Description** for further information on the count
  conversion. Default is `NULL`.

- pref_field_method:

  (optional) a `character` string. The column name in `data` of the
  preferred field method of the species. This argument is only required
  is `field_method` is not `NULL` (i.e. individuals have been counted by
  different methods). Alternatively, this value can be passed in `info`
  (or internally retrieved if the species is listed in the package). See
  above section **Description** for further information on the count
  conversion. Default is `NULL`.

- conversion_A2G:

  (optional) a `character` string. The column name in `data` of the
  count conversion factor of the species. This argument is only required
  if `field_method` is not `NULL` (i.e. individuals have been counted by
  different methods). Alternatively this value can be passed in `info`
  (or internally retrieved if the species is listed in the package). See
  above section **Description** for further information on the count
  conversion. Default is `NULL`.

- rmax:

  (optional) a `character` string. The column name in `data` of the
  species demographic potential (i.e. the relative rate of increase of
  the population). This is the change in log population size between two
  dates and will be used later by
  [`fit_trend()`](https://frbcesab.github.io/popbayes/reference/fit_trend.md).
  Default is `NULL`.

- path:

  a `character` string. The directory to save formatted data. This
  directory must exist and can be an absolute or a relative path.
  Default is the current working directory.

- na_rm:

  a `logical.` If `TRUE`, counts with `NA` values will be removed.
  Default is `FALSE` (returns an error to inform user if `NA` are
  detected).

## Value

An n-elements `list` (where `n` is the number of count series). The name
of each element of this list is a combination of location and species.
Each element of the list is a `list` with the following content:

- `location` a `character` string. The name of the series site.

- `species` a `character` string. The name of the series species.

- `date` a `numerical` vector. The sequence of dates of the series.

- `n_dates` an `integer.` The number of unique dates.

- `stat_methods` a `character` vector. The different stat methods of the
  series.

- `field_methods` (optional) a `character` vector. The different field
  methods of the series.

- `pref_field_method` (optional) a `character` string. The preferred
  field method of the species (`'A'` or `'G'`).

- `conversion_A2G` (optional) a `numeric`. The conversion factor of the
  species used to convert counts to its preferred field method.

- `rmax` a `numeric`. The maximum population growth rate of the species.

- `data_original` a `data.frame`. Original data of the series with
  renamed columns. Some rows may have been deleted (if `na_rm = TRUE`).

- `data_converted` a `data.frame`. Data containing computed boundaries
  of the 95% CI (`lower_ci_conv` and `upper_ci_conv`). If counts have
  been obtained by different field methods, contains also converted
  counts (`count_conv`) based on the preferred field method and
  conversion factor of the species. This `data.frame` will be used by
  the function
  [`fit_trend()`](https://frbcesab.github.io/popbayes/reference/fit_trend.md)
  to fit population models.

**Note:** Some original series can be discarded if one of these two
conditions is met: 1) the series contains only zero counts, and 2) the
series contains only a few dates (\< 4 dates).

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

## Print content of the first count series ----
names(garamba_formatted[[1]])
#>  [1] "location"          "species"           "dates"            
#>  [4] "n_dates"           "stat_methods"      "field_methods"    
#>  [7] "pref_field_method" "conversion_A2G"    "rmax"             
#> [10] "data_original"     "data_converted"   

## Print original data ----
garamba_formatted[[1]]$"data_original"
#>    location               species date stat_method field_method
#> 1   Garamba Alcelaphus buselaphus 1976           S            A
#> 2   Garamba Alcelaphus buselaphus 1983           S            A
#> 3   Garamba Alcelaphus buselaphus 1984           S            A
#> 4   Garamba Alcelaphus buselaphus 1986           S            A
#> 5   Garamba Alcelaphus buselaphus 1991           S            A
#> 6   Garamba Alcelaphus buselaphus 1993           S            A
#> 7   Garamba Alcelaphus buselaphus 1995           S            A
#> 8   Garamba Alcelaphus buselaphus 1998           S            A
#> 9   Garamba Alcelaphus buselaphus 2000           S            A
#> 10  Garamba Alcelaphus buselaphus 2002           S            A
#> 11  Garamba Alcelaphus buselaphus 2003           S            A
#> 12  Garamba Alcelaphus buselaphus 2004           S            A
#> 13  Garamba Alcelaphus buselaphus 2012           T            A
#> 14  Garamba Alcelaphus buselaphus 2014           T            A
#> 15  Garamba Alcelaphus buselaphus 2017           T            A
#>    pref_field_method conversion_A2G   rmax count_orig lower_ci_orig
#> 1                  G          2.302 0.2748       7750          6280
#> 2                  G          2.302 0.2748       1932          1120
#> 3                  G          2.302 0.2748       1224           782
#> 4                  G          2.302 0.2748       1705          1116
#> 5                  G          2.302 0.2748        987           663
#> 6                  G          2.302 0.2748       3444          1290
#> 7                  G          2.302 0.2748       2819          1620
#> 8                  G          2.302 0.2748       1685          1287
#> 9                  G          2.302 0.2748       1169           945
#> 10                 G          2.302 0.2748       1139           907
#> 11                 G          2.302 0.2748       1595          1142
#> 12                 G          2.302 0.2748       1204           811
#> 13                 G          2.302 0.2748        552            NA
#> 14                 G          2.302 0.2748        698            NA
#> 15                 G          2.302 0.2748       1051            NA
#>    upper_ci_orig
#> 1           9220
#> 2           2744
#> 3           1666
#> 4           2294
#> 5           1311
#> 6           5598
#> 7           4018
#> 8           2083
#> 9           1393
#> 10          1371
#> 11          2048
#> 12          1597
#> 13            NA
#> 14            NA
#> 15            NA

## Print converted data ----
garamba_formatted[[1]]$"data_converted"
#>    location               species date stat_method field_method
#> 1   Garamba Alcelaphus buselaphus 1976           S            A
#> 2   Garamba Alcelaphus buselaphus 1983           S            A
#> 3   Garamba Alcelaphus buselaphus 1984           S            A
#> 4   Garamba Alcelaphus buselaphus 1986           S            A
#> 5   Garamba Alcelaphus buselaphus 1991           S            A
#> 6   Garamba Alcelaphus buselaphus 1993           S            A
#> 7   Garamba Alcelaphus buselaphus 1995           S            A
#> 8   Garamba Alcelaphus buselaphus 1998           S            A
#> 9   Garamba Alcelaphus buselaphus 2000           S            A
#> 10  Garamba Alcelaphus buselaphus 2002           S            A
#> 11  Garamba Alcelaphus buselaphus 2003           S            A
#> 12  Garamba Alcelaphus buselaphus 2004           S            A
#> 13  Garamba Alcelaphus buselaphus 2012           T            A
#> 14  Garamba Alcelaphus buselaphus 2014           T            A
#> 15  Garamba Alcelaphus buselaphus 2017           T            A
#>    pref_field_method conversion_A2G   rmax count_conv lower_ci_conv
#> 1                  G          2.302 0.2748  17840.500     14456.560
#> 2                  G          2.302 0.2748   4447.464      2578.240
#> 3                  G          2.302 0.2748   2817.648      1800.164
#> 4                  G          2.302 0.2748   3924.910      2569.032
#> 5                  G          2.302 0.2748   2272.074      1526.226
#> 6                  G          2.302 0.2748   7928.088      2969.580
#> 7                  G          2.302 0.2748   6489.338      3729.240
#> 8                  G          2.302 0.2748   3878.870      2962.674
#> 9                  G          2.302 0.2748   2691.038      2175.390
#> 10                 G          2.302 0.2748   2621.978      2087.914
#> 11                 G          2.302 0.2748   3671.690      2628.884
#> 12                 G          2.302 0.2748   2771.608      1866.922
#> 13                 G          2.302 0.2748   1270.704      1207.169
#> 14                 G          2.302 0.2748   1606.796      1526.456
#> 15                 G          2.302 0.2748   2419.402      2298.432
#>    upper_ci_conv field_method_conv
#> 1      21224.440                 G
#> 2       6316.688                 G
#> 3       3835.132                 G
#> 4       5280.788                 G
#> 5       3017.922                 G
#> 6      12886.596                 G
#> 7       9249.436                 G
#> 8       4795.066                 G
#> 9       3206.686                 G
#> 10      3156.042                 G
#> 11      4714.496                 G
#> 12      3676.294                 G
#> 13      1524.845                 G
#> 14      1928.155                 G
#> 15      2903.282                 G
```
