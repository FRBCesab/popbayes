## Test environments

* Local
  * macOS 11.6 install, R 4.1.2
* GitHub Actions
  * Mac OS X 10.15, r-release (R 4.1.2), r-oldrel
  * Windows Server 2019, r-devel, r-release (R 4.1.2), r-oldrel
  * Ubuntu 20.04.3 LTS, r-devel, r-release (R 4.1.2), r-oldrel
* WinBuilder
  * r-devel - 1 note (new submission)
  * r-release - 1 note (new submission)
  * r-oldrel - 1 note (new submission)
* R-hub
  * Windows Server 2008 R2 SP1, R-devel, 32/64 bit - 1 note (new submission)
  * Ubuntu Linux 20.04.1 LTS, R-release, GCC - 1 note (new submission)
  * Fedora Linux, R-devel, clang, gfortran - ERROR: the system requirement 'JAGS' is not available on this platform.


## R CMD check results

0 error | 0 warning | 1 note

* checking CRAN incoming feasibility ... NOTE
  
  Maintainer: 'Nicolas Casajus <nicolas.casajus@fondationbiodiversite.fr>' - New submission

* No mis-spelled words in 'DESCRIPTION'


## Downstream dependencies

There are currently no downstream dependencies for this package.


## Resubmit comments

* Size of tarball: 21557165 bytes - Fixed: size is now 2.1 MB
* Add `\value{No return value}` in `diagnostic()`, `plot_series()` and `plot_trend()` to document function output
* Replace all `\dontrun{}` by `\donttest{}` in JAGS examples (model fitting takes > 5s)
* Add `on.exit(par(opar))` in `plot_series()` and `plot_trend()` to restore user graphical parameters
* Thanks!
