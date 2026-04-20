# Compute rmax from adult female body mass

The demographic potential of a species is limited. The intrinsic rate of
increase `rmax` is the maximum increase in log population size that a
species can attain in a year. According to Sinclair (2003), it is
related to the body mass of adult females by: \\1.375 \times
W^{-0.315}\\

## Usage

``` r
w_to_rmax(w)
```

## Arguments

- w:

  a numerical vector. Adult female body mass (in kg).

## Value

A numerical vector of `rmax` values.

## References

Sinclair (2013) Mammal population regulation, keystone processes and
ecosystem dynamics. *Philosophical Transactions: Biological Sciences*,
**358**, 1729-1740.

## Examples

``` r
## Set adult female body mass ----
body_masses <- c(55, 127)

## Add species names ----
names(body_masses) <- c("Impala", "Tiang")

## Compute species rmax ----
w_to_rmax(body_masses)
#>    Impala     Tiang 
#> 0.3891244 0.2989541 
```
