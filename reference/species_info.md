# Species information dataset

This dataset contains information about 15 African mammal species. It
can be used in the function
[`format_data()`](https://frbcesab.github.io/popbayes/reference/format_data.md)
to convert individual counts estimated from a field method to a
preferred field method. The field method can be `A` (aerial counts) or
`G` (ground counts). See
[`format_data()`](https://frbcesab.github.io/popbayes/reference/format_data.md)
for further information. It also contains the maximum population growth
rate (i.e. the maximum change in log population size).

User can take this dataset as a template to add information for missing
species. Note that only `species`, `pref_field_method`,
`conversion_A2G`, and `rmax` are required.

## Usage

``` r
species_info
```

## Format

A `data.frame` with 15 rows (African mammals species) and the following
variables:

- order:

  the order of the species

- family:

  the family of the species

- species:

  the species binomial name

- english:

  the species English name

- french:

  the species French name

- category:

  the detectability category of the species. One of `MLB` for
  Medium-sized Light and Brown species (20-150kg), `LLB` for Large Light
  and Brown species (\>150kg), `LD` for Large Dark (\>150kg),
  `Elephant`, and `Giraffe`

- pref_field_method:

  the preferred field method of the species. One of `A` for Aerial
  counts, and `G` for Ground counts

- conversion_A2G:

  the conversion multiplicative factor (corresponding to the
  detectability category) used to convert aerial to ground counts

- rmax:

  the maximum population growth rate

## Examples

``` r
data("species_info")
species_info
#>           order       family                species             english
#> 1  Artiodactyla      Bovidae     Aepyceros melampus              Impala
#> 2  Artiodactyla      Bovidae  Alcelaphus buselaphus          Hartebeest
#> 3  Artiodactyla      Bovidae  Connochaetes taurinus     Blue wildebeest
#> 4  Artiodactyla      Bovidae     Damaliscus lunatus               Tiang
#> 5  Artiodactyla      Bovidae     Eudorcas rufifrons Red-fronted gazelle
#> 6  Artiodactyla   Giraffidae Giraffa camelopardalis             Giraffe
#> 7  Artiodactyla      Bovidae    Hippotragus equinus       Roan Antelope
#> 8  Artiodactyla      Bovidae   Kobus ellipsiprymnus           Waterbuck
#> 9  Artiodactyla      Bovidae              Kobus kob                 Kob
#> 10  Proboscidea Elephantidae     Loxodonta africana    Savanna Elephant
#> 11 Artiodactyla      Bovidae         Ourebia ourebi               Oribi
#> 12 Artiodactyla      Bovidae        Redunca redunca      Bohor reedbuck
#> 13 Artiodactyla      Bovidae        Syncerus caffer     African Buffalo
#> 14 Artiodactyla      Bovidae  Tragelaphus derbianus         Giant eland
#> 15 Artiodactyla      Bovidae   Tragelaphus scriptus            Bushbuck
#>              french category pref_field_method conversion_A2G   rmax
#> 1            Impala      MLB                 G          6.747 0.4010
#> 2            Bubale      LLB                 G          2.302 0.2748
#> 3         Gnou bleu      LLB                 G          2.302 0.2679
#> 4              Topi      MLB                 G          6.747 0.2990
#> 5  Gazelle Rufifron      MLB                 G          6.747 0.5270
#> 6            Girafe  Giraffe                 A          3.011 0.1750
#> 7       Hippotrague      LLB                 G          2.302 0.2420
#> 8      Cobe Defassa      MLB                 G          6.747 0.2702
#> 9    Cobe de Buffon      MLB                 G          6.747 0.3802
#> 10         Elephant Elephant                 A          0.659 0.1120
#> 11           Ourebi      MLB                 G          6.747 0.5988
#> 12          Redunca      MLB                 G          6.747 0.4010
#> 13           Buffle       LD                 A          0.561 0.2080
#> 14   Eland de derby      LLB                 G          2.302 0.1500
#> 15    Guib Harnache      MLB                 G          6.747 0.4487
```
