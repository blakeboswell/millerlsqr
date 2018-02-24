
update_miller <- function(X, y, weight, qr, intercept = FALSE) {
  ## the intercept logic is not used .. X will already have an
  ## intercept column if the intercept is to be included

  p    <- length(qr$D) + as.integer(intercept)
  k    <- 1 + as.integer(intercept)
  xrow <- rep(0.0, p)
  
  for(i in 1:nrow(X)) {
    if(intercept){
      xrow[1] <- 1.0
    }
    xrow[k:p] <- X[i, ]
    qr$includ(xrow, y[[i]], weight[[i]])
  }
  
}


miller_lsqr <- function(np) {
  
  small <- 1.e-69
  zero  <- 0.0

  nrbar  <- np * (np - 1) / 2
  D      <- rep(zero, np)
  rbar   <- rep(zero, nrbar)
  thetab <- rep(zero, np)
  sserr  <- zero
  
  tol     <- rep(zero, np)
  tol_set <- FALSE
  nobs    <- 0L
  
  
  some_checks <- function(np, nrbar) {
    ier <- 0L
    if(np < 1){ ier <- 1L }
    if(nrbar < np * (np - 1) / 2) {
      ier <- ier + 2
    }
    ier
  }
  
  
  includ <- function(xrow, yelem, weight = 1) {
    
    ##  ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2
    ##  Modified from algorithm AS 75.1
    ##
    ##  Calling this routine updates `D`, `rbar`, `thetab`` and `sserr`` by the
    ##  inclusion of `xrow`, `yelem` with the specified `weight`.   The number
    ##  of columns (variables) may exceed the number of rows (cases).
    
    ##  Some checks.
    
    if(length(xrow) < np){
      warning("invalid xrow encountered")
      return(invisible())
    }
    
    ier <- some_checks(np, nrbar)
    if(ier != 0) {
      return(invisible())
    }
    
    nobs  <<- nobs + 1
    nextr <- 1L
    
    for(i in 1:np) {
      
    ##  Skip unnecessary transformations.   Test on exact zeroes must be
    ##  used or stability can be destroyed.
      
      if(abs(weight) < small) {
        return(invisible())
      }
      
      xi <- xrow[i]
      if(abs(xi) < small) {
        nextr <- nextr + np - i
        next
      }
      
      di   <- D[i]
      wxi  <- weight * xi
      dpi  <- di + wxi * xi
      cbar <- di / dpi
      sbar <- wxi / dpi
      
      weight <- cbar * weight
      
      D[i] <<- dpi
      
      k <- i + 1
      while (k <= np) {
        xk          <- xrow[k]
        xrow[k]     <- xk - xi * rbar[nextr]
        rbar[nextr] <<- cbar * rbar[nextr] + sbar * xk
        nextr       <- nextr + 1
        k <- k + 1
      }
      
      xk    <- yelem
      yelem <- xk - xi * thetab[i]
      thetab[i] <<- cbar * thetab[i] + sbar * xk
    }
    
    ## `yelem` * `sqrt(weight)` is now equal to Brown & Durbin's recursive residual.
    
    sserr <<- sserr + weight * yelem ^ 2
  }
  
  
  tolset <- function() {
    
    ##  ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2
    ##
    ##  Sets up array TOL for testing for zeroes in an orthogonal
    ##  reduction formed using AS75.1.
    
    ## Some checks.
    
    ier <- some_checks(np, nrbar)
    if(ier != 0) {
      return(invisible())
    }
    
    ##  EPS is a machine-dependent constant.   For compilers which use
    ##  the IEEE format for floating-point numbers, recommended values
    ##  are 1.E-06 for single precision and 1.D-12 for double precision.
    
    eps   <- 1.0e-12
    
    ##  Set `tol[i]` = sum of absolute values in column `i` of `rbar` after
    ##  scaling each element by the square root of its row multiplier.
    
    work  <- sqrt(D)

    for(col in 1:np) {
      pos   <- col - 1
      total <- work[col]
      row   <- 1
      
      while(row < col) {
        total <- total + abs(rbar[pos]) * work[row]
        pos   <- pos + np - row - 1
        row   <- row + 1
      }
      
      tol[col] <<- eps * total
    }
    
    tol_set <<- TRUE
    
  }

  
    
  singchk <- function() {
    
    ##  ALGORITHM AS274  APPL. STATIST. (1992) VOL.41, NO. 2
    ##
    ##  Checks for singularities, reports, and adjusts orthogonal
    ##  reductions produced by AS75.1.
    
    ## Some checks.
    
    ier <- some_checks(np, nrbar)
    if(ier != 0) {
      return(invisible())
    }
    
    linedep <- rep(FALSE, np)
    work    <- sqrt(D)
    
    for(col in 1:np) {
      
      ##  Set elements within `rbar` to `zero`` if they are less than `tol[col]` in
      ##  absolute value after being scaled by the square root of their row
      ##  multiplier.
      
      temp <- tol[col]
      pos  <- col - 1
      row  <- 1
      
      while (row < col) {
        if(abs(rbar[pos]) * work[row] < temp) {
          rbar[pos] <- zero
        }
        pos <- pos + np - row - 1
        row <- row + 1
      }
      
      ##  If diagonal element is near zero, set it to zero, set appropriate
      ##  element of `lindep`, and use `includ` to augment the projections in
      ##  the lower rows of the orthogonalization.
      
      if(work[col] <= temp) {
        
        linedep[col] <- TRUE
        ier       <<- ier - 1
        
        if (col < np) {
          pos2 <- pos + np - col + 1
          x    <- rep(zero, np)
          
          for(k in (col + 1):np){
            x[k] <- rbar[pos+k-col]
          }
          
          y      <- thetab[col]
          weight <- D[col]
          
          for(k in (pos+1):(pos2-1)){
            rbar[k] <<- zero
          }
          
          D[col]      <<- zero
          thetab[col] <<- zero
          includ(x, y, weight)
          nobs <<- nobs - 1
          
        } else {
          sserr <<- sserr + D[col] * thetab[col] ^ 2
        }
        
      }
    }
  }
  
  
  regcf <- function() {
    
    ##  ALGORITHM AS274  APPL. STATIST. (1992) VOL 41, NO. x
    ##
    ##  Modified version of AS75.4 to calculate regression coefficients
    ##  for the first NREQ variables, given an orthogonal reduction from
    ##  AS75.1.
    
    ## Some checks.
    
    ier <- some_checks(np, nrbar)
    if(ier != 0) {
      return(list(beta = NULL, ier = ier))
    }
    
    beta <- rep(zero, np)
    
    for(i in np:1) {
      if(D[i] < small) {
        beta[i] <- zero
        D[i]    <<- zero
      } else {
        
        beta[i] <- thetab[i]
        nextr   <- (i - 1) * (np + np - i) / 2 + 1
        
        j <- i + 1
        while(j <= np) {
          beta[i] <- beta[i] - rbar[nextr] * beta[j]
          nextr   <- nextr + 1
          j <- j + 1
        }
      }
    }
    
    list(
      beta = beta,
      ier  = ier
    )
  }
  
  
  e <- environment()
  class(e) <- "online_qr"
  e
  
}

