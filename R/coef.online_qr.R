
coef.online_qr <- function(qr, nvar = NULL, ...) {
    
  p <- length(qr$D)
    
  if (is.null(nvar)) {
    nvar <- p
  }

  if (nvar < 1 | nvar > p){
    stop("Invalid value of `nvar`")
  }

  if (!qr$tol_set){
    qr$singchk()
  }

  tmp <- qr$regcf()
  
  if (tmp$ier != 0) {
    stop("Error in REGCF: can't happen")
  }
  
  tmp$beta
  
}
