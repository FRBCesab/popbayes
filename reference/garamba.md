# African mammals survey in the Garamba National Park

This dataset contains individual counts of 10 African mammal species in
the Garamba National Park (Democratic Republic of the Congo) from 1976
to 2017.

## Usage

``` r
garamba
```

## Format

A `data.frame` with 141 rows (counts) and the following 8 variables:

- location:

  the location of the survey (Garamba)

- species:

  the binomial name of the species

- date:

  the date of the survey

- stat_method:

  the method used to estimate individuals counts. One of `T` (total
  counts), `G` (guesstimate), and `S` (sampling counts)

- field_method:

  the field method used to collect data. One of `A` (aerial counts), and
  `G` (ground counts)

- count:

  number of individuals

- lower_ci:

  lower boundary of the 95% confidence interval of the counts (only for
  sampling counts)

- upper_ci:

  upper boundary of the 95% confidence interval of the counts (only for
  sampling counts)

- pref_field_method:

  the preferred field method of the species. One of `A` for Aerial
  counts, and `G` for Ground counts

- conversion_A2G:

  the conversion multiplicative factor (corresponding to the
  detectability category) used to convert aerial to ground counts

- rmax:

  the maximum population growth rate

## References

Hillman Smith K & Kalpers J (2015) *Garamba Conservation in Peace &
War*. Hillman Smith publisher, 448pp. ISBN: 9789966185105.  
Monico M (2014) *Aerial Survey Report March 2014 - Garamba National
Park, DRC*. African Parks Network/ICCN/Pan-African Elephant Aerial
Survey, 42pp.  
Spies K *et al.* (2017) *Aerial Survey Report April 2017 - Garamba
National Park, DRC*. African Parks Network/EU/ICCN. 38pp.

## Examples

``` r
data("garamba")
head(garamba, 30)
#>    location                species date stat_method field_method count lower_ci
#> 1   Garamba  Alcelaphus buselaphus 1976           S            A  7750     6280
#> 2   Garamba  Alcelaphus buselaphus 1983           S            A  1932     1120
#> 3   Garamba  Alcelaphus buselaphus 1984           S            A  1224      782
#> 4   Garamba  Alcelaphus buselaphus 1986           S            A  1705     1116
#> 5   Garamba  Alcelaphus buselaphus 1991           S            A   987      663
#> 6   Garamba  Alcelaphus buselaphus 1993           S            A  3444     1290
#> 7   Garamba  Alcelaphus buselaphus 1995           S            A  2819     1620
#> 8   Garamba  Alcelaphus buselaphus 1998           S            A  1685     1287
#> 9   Garamba  Alcelaphus buselaphus 2000           S            A  1169      945
#> 10  Garamba  Alcelaphus buselaphus 2002           S            A  1139      907
#> 11  Garamba  Alcelaphus buselaphus 2003           S            A  1595     1142
#> 12  Garamba  Alcelaphus buselaphus 2004           S            A  1204      811
#> 13  Garamba  Alcelaphus buselaphus 2012           T            A   552       NA
#> 14  Garamba  Alcelaphus buselaphus 2014           T            A   698       NA
#> 15  Garamba  Alcelaphus buselaphus 2017           T            A  1051       NA
#> 16  Garamba Giraffa camelopardalis 1976           S            A   350      100
#> 17  Garamba Giraffa camelopardalis 1983           S            A   175       12
#> 18  Garamba Giraffa camelopardalis 1984           S            A   237       93
#> 19  Garamba Giraffa camelopardalis 1986           S            A   153       13
#> 20  Garamba Giraffa camelopardalis 1991           S            A   346      143
#> 21  Garamba Giraffa camelopardalis 1993           S            A   347      347
#> 22  Garamba Giraffa camelopardalis 1995           S            A   178      178
#> 23  Garamba Giraffa camelopardalis 1998           S            A   144       71
#> 24  Garamba Giraffa camelopardalis 2000           S            A   118       54
#> 25  Garamba Giraffa camelopardalis 2002           S            A    62       49
#> 26  Garamba Giraffa camelopardalis 2003           S            A    62       62
#> 27  Garamba Giraffa camelopardalis 2004           S            A   185       33
#> 28  Garamba Giraffa camelopardalis 2012           T            A    11       NA
#> 29  Garamba Giraffa camelopardalis 2014           T            A    42       NA
#> 30  Garamba Giraffa camelopardalis 2017           T            A    42       NA
#>    upper_ci pref_field_method conversion_A2G   rmax
#> 1      9220                 G          2.302 0.2748
#> 2      2744                 G          2.302 0.2748
#> 3      1666                 G          2.302 0.2748
#> 4      2294                 G          2.302 0.2748
#> 5      1311                 G          2.302 0.2748
#> 6      5598                 G          2.302 0.2748
#> 7      4018                 G          2.302 0.2748
#> 8      2083                 G          2.302 0.2748
#> 9      1393                 G          2.302 0.2748
#> 10     1371                 G          2.302 0.2748
#> 11     2048                 G          2.302 0.2748
#> 12     1597                 G          2.302 0.2748
#> 13       NA                 G          2.302 0.2748
#> 14       NA                 G          2.302 0.2748
#> 15       NA                 G          2.302 0.2748
#> 16      600                 A          3.011 0.1750
#> 17      338                 A          3.011 0.1750
#> 18      381                 A          3.011 0.1750
#> 19      293                 A          3.011 0.1750
#> 20      549                 A          3.011 0.1750
#> 21      766                 A          3.011 0.1750
#> 22      388                 A          3.011 0.1750
#> 23      217                 A          3.011 0.1750
#> 24      182                 A          3.011 0.1750
#> 25       75                 A          3.011 0.1750
#> 26      137                 A          3.011 0.1750
#> 27      337                 A          3.011 0.1750
#> 28       NA                 A          3.011 0.1750
#> 29       NA                 A          3.011 0.1750
#> 30       NA                 A          3.011 0.1750
```
