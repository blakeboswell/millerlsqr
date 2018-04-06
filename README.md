# `millerlsqr`

## About

The script [R/miller_lsqr.R](https://github.com/blakeboswell/millerlsqr/blob/master/R/miller_lsqr.R) is a translation of Alan Miller's [AS274 Fortran implementation](https://github.com/cran/biglm/blob/master/src/boundedQRf.f) employed in the R package `biglm`.

The other scripts / functions are copied and minimally modified from `biglm` for ease of testing purposes.  The scripts are also structured as a package for convenience .. I don't recommend installing and using this as a package.

The purpose of this script is to expose AS274 in base R.  I have to admit though, my opinion of Fortran has improved as a result of this effort.
