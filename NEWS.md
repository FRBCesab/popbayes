# popbayes (development version)

* popbayes now uses cli instead of usethis for error messages (@olivroy, #33).

# popbayes 1.2

* Fix issue in `format_data()` when user imports a `tibble`
* Fix some typos in documentation

# popbayes 1.1

* Change `stat_method` category `G` by `X` (eXpert knowledge) to avoid confusion
with the category `G` (Ground counts) in the `field_method` variable
* Set arguments `pref_field_method`, `field_method`, `conversion_A2G`, and `rmax` 
to `NULL` by default in `format_data()`
* Allow `NA` values in column `field_method` if `stat_method = "X"`
* Function `format_data()` now works at the count series level (not the whole 
data set). This allow users to define different values for these arguments for
different counts series (with the same species)
* Function `format_data()` returns an error if some confident interval bounds 
are strictly equal

# popbayes 1.0

* First release of the package
* Submission to CRAN
