## About

[miller_lsqr.R](https://github.com/blakeboswell/millerlsqr/blob/master/R/miller_lsqr.R) is a translation of the [AS274 Fortran implementation](https://github.com/cran/biglm/blob/master/src/boundedQRf.f) in the R package `biglm`.

The other scripts / functions are copied and minimally modified from `biglm` for ease of testing purposes.  The scripts are also structured as a package for convenience .. I don't recommend installing and using this as a package.

## Why?

I wanted to understand AS274.  Translating from Fortran to R improved that understanding.