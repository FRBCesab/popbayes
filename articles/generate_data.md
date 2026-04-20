# Get species information

## Calculating rmax

The demographic potential of a species is limited. The intrinsic rate of
increase $rmax$ is the maximum increase in log population size that a
species can attain in a year. According to Sinclair (2003), it is
related to the body mass of adult females $W$ by the formula:

$$rmax = 1.375 \times W^{- 0.315}$$

Body masses are found in the literature in publications such as Kingdon
& Hoffman (2013), Cornelis *et al.* (2014), Illius & Gordon (1992),
Sinclair (1996), Suraud *et al.* (2012), or Foley & Faust (2010).
Alternatively, $rmax$ can be obtained from specific demographic
analyses. In the following table, we have listed the $rmax$ values
obtained with one or the other method giving precedence to specific
analyses when available.

``` r
species_name <- c("Impala", "Tiang", "Blue wildebeest", "Roan", "Buffalo", 
                  "Eland", "Giraffe", "Elephant")

adult_female_body_mass <- c(55, 127, 230, 250, 400, 450, 702, 2873)

species <- data.frame(adult_female_body_mass, row.names = species_name)

species$rmax <- 1.375 * adult_female_body_mass ^ (-0.315)
species["Eland", "rmax"]    <- 0.150                 # Sinclair (1996)
species["Elephant", "rmax"] <- 0.07                  # Foley & Faust (2010)
species["Giraffe", "rmax"]  <- 0.125                 # Suraud et al. (2012)

species
#>                 adult_female_body_mass      rmax
#> Impala                              55 0.3891244
#> Tiang                              127 0.2989541
#> Blue wildebeest                    230 0.2479468
#> Roan                               250 0.2415192
#> Buffalo                            400 0.2082830
#> Eland                              450 0.1500000
#> Giraffe                            702 0.1250000
#> Elephant                          2873 0.0700000
```

## Building conversion table

Aerial and ground counts are not directly comparable because some
species are better detected from the ground and others from the air. It
is generally considered that small light species are better detected
from the ground while large dark species are better detected from the
air. We took advantage of series of counts carried out in parallel using
the two field methods to calculate conversion factors to be applied to
aerial counts to obtain ground count equivalents. This permits
reconciling the two types of counts within a mixed series. We also
specify the preferred field method for each category of species.

``` r
categories <- c("Medium light and brown species (20-150kg)",
                "Large light and brown species (>150kg)",
                "Large dark (>150kg)", "Giraffe", "Elephant")

short_names        <- c("MLB", "LLB", "LD", "Giraffe", "Elephant")
conversion_facts   <- c(6.747, 2.302, 0.561, 3.011, 0.659)
pref_field_methods <- c("G", "G", "A", "A", "A")

category_info <- data.frame("category"          = categories,
                            "acronym"           = short_names,
                            "conversion_fact"   = conversion_facts, 
                            "pref_field_method" = pref_field_methods)
category_info
#>                                    category  acronym conversion_fact
#> 1 Medium light and brown species (20-150kg)      MLB           6.747
#> 2    Large light and brown species (>150kg)      LLB           2.302
#> 3                       Large dark (>150kg)       LD           0.561
#> 4                                   Giraffe  Giraffe           3.011
#> 5                                  Elephant Elephant           0.659
#>   pref_field_method
#> 1                 G
#> 2                 G
#> 3                 A
#> 4                 A
#> 5                 A
```

## Categorizing species

Here we relate the species to a color/size category.

``` r
species$"category" <- c("MLB", "MLB", "LLB", "LLB", "LD", "LLB", "Giraffe", 
                        "Elephant")
species
#>                 adult_female_body_mass      rmax category
#> Impala                              55 0.3891244      MLB
#> Tiang                              127 0.2989541      MLB
#> Blue wildebeest                    230 0.2479468      LLB
#> Roan                               250 0.2415192      LLB
#> Buffalo                            400 0.2082830       LD
#> Eland                              450 0.1500000      LLB
#> Giraffe                            702 0.1250000  Giraffe
#> Elephant                          2873 0.0700000 Elephant
```

  

## References

- Cornelis D *et al.* (2014) Species account: African buffalo (*Syncerus
  caffer*). In: *Ecology, Evolution and Behaviour of Wild Cattle:
  Implications for Conservation* (Eds M Melletti & J Burton). Cambridge
  University Press, Cambridge. DOI:
  [10.1017/CBO9781139568098](https://doi.org/10.1017/CBO9781139568098).

- Foley CAH & Faust LJ (2010) Rapid population growth in an elephant
  *Loxodonta africana* population recovering from poaching in Tarangire
  National Park, Tanzania. *Oryx*, **44**, 205-212. DOI:
  [10.1017/S0030605309990706](https://doi.org/10.1017/S0030605309990706).

- Illius AW & Gordon IJ (1992) Modelling the nutritional ecology of
  ungulate herbivores: evolution of body size and competitive
  interactions. *Oecologia*, **89**, 428-434. DOI:
  [10.1017/S0030605309990706](https://doi.org/10.1007/BF00317422).

- Kingdon J & Hoffman M (2013) *Mammals of Africa. Volume VI: Pigs,
  Hippopotamuses, Chevrotain, Giraffes, Deer and Bovids*. Bloomsbury
  Publishing, London, United Kingdom, 680 pp.

- Sinclair ARE (1996) Mammal populations: fluctuation, regulation, life
  history theory, and their implications for conservation. In:
  *Frontiers of population ecology* (Eds RB Floyd & AW Sheppard),
  pp. 127–154. CSIRO: Melbourne, Australia.

- Sinclair (2003) Mammal population regulation, keystone processes and
  ecosystem dynamics. *Philosophical Transactions: Biological Sciences*,
  **358**, 1729-1740. DOI:
  [10.1098/rstb.2003.1359](https://doi.org/10.1098/rstb.2003.1359).

- Suraud JP *et al.* (2012) Higher than expected growth rate of the
  endangered West African giraffe *Giraffa camelopardalis peralta*: a
  successful human–wildlife cohabitation. *Oryx*, **46**, 577-583. DOI:
  [10.1017/S0030605311000639](https://doi.org/10.1017/S0030605311000639).
