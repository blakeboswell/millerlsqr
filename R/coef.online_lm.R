coef.online_lm <- function(object, ...) {
  if (!object$qr$tol_set) {
    object$qr$singchk()
  }
  rval <- coef(object$qr)
  rval[object$qr$D == 0] <- NA
  names(rval) <- object$names
  rval
}
