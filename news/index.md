# Changelog

## popbayes 1.3

- `popbayes` now uses `cli` instead of `usethis` for error messages
  ([@olivroy](https://github.com/olivroy),
  [\#33](https://github.com/FRBCesab/popbayes/pull/33)).
- Improve website ([@olivroy](https://github.com/olivroy),
  [\#34](https://github.com/FRBCesab/popbayes/pull/34))
- Fix error messages in
  [`fit_trend()`](https://frbcesab.github.io/popbayes/reference/fit_trend.md)
  ([@fjoyce](https://github.com/fjoyce),
  [\#35](https://github.com/FRBCesab/popbayes/issues/35))

## popbayes 1.2

- Fix issue in
  [`format_data()`](https://frbcesab.github.io/popbayes/reference/format_data.md)
  when user imports a `tibble`
- Fix some typos in documentation

## popbayes 1.1

CRAN release: 2022-03-04

- Change `stat_method` category `G` by `X` (eXpert knowledge) to avoid
  confusion with the category `G` (Ground counts) in the `field_method`
  variable
- Set arguments `pref_field_method`, `field_method`, `conversion_A2G`,
  and `rmax` to `NULL` by default in
  [`format_data()`](https://frbcesab.github.io/popbayes/reference/format_data.md)
- Allow `NA` values in column `field_method` if `stat_method = "X"`
- Function
  [`format_data()`](https://frbcesab.github.io/popbayes/reference/format_data.md)
  now works at the count series level (not the whole data set). This
  allow users to define different values for these arguments for
  different counts series (with the same species)
- Function
  [`format_data()`](https://frbcesab.github.io/popbayes/reference/format_data.md)
  returns an error if some confident interval bounds are strictly equal

## popbayes 1.0

CRAN release: 2021-11-05

- First release of the package
- Submission to CRAN
