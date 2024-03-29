### Numerical gradient====
grad <- function(mod, par, ..., Interval=1e-6) {
  mat  <- matrix(par, nrow=length(par), ncol=length(par))
  for (i in 1:ncol(mat)) {
    mat[i,i] <- par[i] + Interval
  }
  df <- vector("numeric", length=length(par))
  f_x <- mod(par, ...)
  for (i in 1:ncol(mat)) {
    df[i] <- (mod(mat[,i], ...) - f_x) / Interval
  }
  df[which(!is.finite(df))] <- 0
  return(df)
  return(df)
}
steep <- function(mod, par, ..., est, df) {
  mod(est - {par[1]*df*{(abs(df) >  sd(df)) * 1}} -
            {par[2]*df*{(abs(df) <= sd(df)) * 1}}, ...)
}
steepA <- function(mod, par, ..., est, df) {
  mod(est - {par*df}, ...)
}
Gammas <- function(mod, par, ..., est, df, vi, si, iter, epsilon) {
  v    <- {{par[1] * vi} + {{1 - par[1]} * df     }}/{1 - {par[1] ^ iter}}
  s    <- {{par[2] * si} + {{1 - par[2]} * {df*df}}}/{1 - {par[2] ^ iter}}
  newV <- est - {{par[3]*v}/{epsilon+sqrt(s)}}
  return(mod(newV, ...))
}
sHAR <- function(par, mod, ...) {
  theta <- rnorm(length(par))
  d     <- theta / sqrt(sum(theta*theta))
  u     <- suppressWarnings(optim(.5, steepA, mod=mod, df=d, ..., est=par,
                                  method="Nelder-Mead")$par)
  prop  <- par - {u * d}
  return(list("prop"=prop, "u"=u, "d"=d))
}
reltol <- function(x, tol) tol * (abs(x) + tol)

### Adam====
ADAM <- function(fn, startvalue, gr, ..., Interval=1e-6,
                 maxit=100, tol=1e-8, verbose=T) {
  # Opening message
  if(verbose==T) {
    cat("Steepest Adam will run for ", maxit, " iterations at most.\n\n", sep="")
    startTime = proc.time()
  }

  # Initial settings
  epsilon <- Interval
  if(verbose==T) {
    pb <- txtProgressBar(min=0, max=maxit, style=3)
  }
  par     <- s <- v <- G <- matrix(NA, nrow=maxit, ncol=length(startvalue))
  GG      <- matrix(NA, nrow=maxit, ncol=3)
  f       <- vector("numeric", length=maxit)
  convergence = F

  # First step
  G[1,] <- gr(fn, startvalue, ..., Interval=Interval)
  alpha <- suppressWarnings(unlist(optim(0, fn=steepA, mod=fn, df=G[1,],
                                         ..., est=startvalue,
                                         method="Nelder-Mead")$par))
  GG[1,]   <- c(.5, .5, alpha)
  v[1,]    <- G[1,]*GG[1,1]
  s[1,]    <- G[1,]*GG[1,2]
  par[1,]  <- startvalue - {alpha*G[1,]}
  f[1]     <- fn(par[1,], ...)
  if(verbose==T) {
    setTxtProgressBar(pb, 1)
  }

  # Run ADAM algorithm
  for(i in 2:maxit) {
    G[i,]  <- gr(fn, par[i-1,], ..., Interval=Interval)
    GG[i,] <- suppressWarnings(unlist(optim(par=c(0,0,0), fn=Gammas, mod=fn,
                                             df=G[i,], ..., est=par[i-1,],
                                             vi=v[i-1,], si=s[i-1,], iter=i,
                                             epsilon=epsilon)$par))
    if (sum(abs(GG[i,])) == 0) {
      GG[i,] <- suppressWarnings(unlist(ucminf(par=c(0,0,0), fn=Gammas, mod=fn,
                                               df=G[i,], ..., est=par[i-1,],
                                               vi=v[i-1,], si=s[i-1,], iter=i,
                                               epsilon=epsilon)$par))
    }
    gammav  <- GG[i,1]
    gammas  <- GG[i,2]
    alpha   <- GG[i,3]
    v[i,]   <- {{gammav * v[i-1,]} + {{1 - gammav} *  G[i,]       }}/{1 - {gammav^i}}
    s[i,]   <- {{gammas * s[i-1,]} + {{1 - gammas} * {G[i,]*G[i,]}}}/{1 - {gammas^i}}
    par[i,] <- par[i-1,] - {{alpha*v[i,]}/{epsilon+sqrt(s[i,])}}
    f[i]    <- fn(par[i,], ...)
    # Check convergence
    if (reltol(f[i], tol) > abs(f[i] - f[i-1])) {
      convergence = T
      if (verbose==T) {
        setTxtProgressBar(pb, maxit)
        cat("\nConvergence achieved!")
      }
      break
    } else {
        if (verbose==T) {
          setTxtProgressBar(pb, i)
        }
      }
  }
  if ({convergence==F} & {verbose==T}) {
    cat("\nConvergence may not have been achieved!")
  }
  if (i < maxit) {
    f   <- f[-c({i+1}:maxit)]
    par <- as.matrix(par[-c({i+1}:maxit),])
    G   <- G[-c({i+1}:maxit),]
    s   <- s[-c({i+1}:maxit),]
    v   <- v[-c({i+1}:maxit),]
    GG  <- GG[-c({i+1}:maxit),]
  }
  if (verbose==T) {
    close(pb)
    cat("\n")
    # Final messages
    stopTime = proc.time()
    elapsedTime = stopTime - startTime
    cat("It took ",round(elapsedTime[3],2),
        " secs for the run to finish. \n", sep="")
  }

  # Return results
  Results <- list("Cost"=f, "Estimates"=par, "Gradients"=G,
                  "DSG"=s, "Momentum"=v, "par"=par[which.min(f),],
                  "value"=f[which.min(f)], "steps"=GG)
  return(Results)
}

### Steepest 2-group Gradient Descent====
STDM <- function(fn, startvalue, gr, ..., Interval=1e-6, maxit=100, tol=1e-8,
                 verbose=T) {
  # Opening message
  if(verbose==T) {
    cat("Steepest 2-group Gradient Descent will run for ",
        maxit, " iterations at most.\n\n", sep="")
    startTime = proc.time()
  }

  # Initial settings
  if(verbose==T) {
    pb   <- txtProgressBar(min=0, max=maxit, style=3)
  }
  #par  <- df <- array(dim = c(maxit,length(startvalue)))
  par  <- df <- matrix(NA, nrow=maxit, ncol=length(startvalue))
  f    <- vector("numeric", length=maxit)
  step <- array(dim = c(maxit,2))
  convergence = F

  # First step
  df[1,]   <- gr(fn, startvalue, ..., Interval=Interval)
  step[1,] <- suppressWarnings(unlist(optim(c(0,0), fn=steep, mod=fn,
                                            ..., df=df[1,], est=startvalue,
                                            method="Nelder-Mead")$par))
  par[1,]   <- startvalue -
               {step[1,1]*df[1,]*{(abs(df[1,]) >  sd(df[1,])) * 1}} -
               {step[1,2]*df[1,]*{(abs(df[1,]) <= sd(df[1,])) * 1}}
  f[1]      <- fn(par[1,], ...)
  if(verbose==T) {
    setTxtProgressBar(pb, 1)
  }

  # Start estimation
  for(run in 2:maxit) {
    # Calculate gradient and estimate step parameter
    df[run,] <- gr(fn, par[run-1,], ..., Interval=Interval)
    step[run,] <- suppressWarnings(unlist(optim(step[run-1,], fn=steep, mod=fn,
                                                ..., df=df[run,], est=par[run-1,],
                                                method="Nelder-Mead")$par))
    if (sum(abs(step[run,])) == 0) {
      step[run,] <- suppressWarnings(unlist(ucminf(par=c(0,0), fn=steep, mod=fn,
                                                   ..., df=df[run,],
                                                   est=par[run-1,])$par))
    }
    par[run,] <- par[run-1,] -
                 {step[run,1]*df[run,]*{(abs(df[run,]) >  sd(df[run,])) * 1}} -
                 {step[run,2]*df[run,]*{(abs(df[run,]) <= sd(df[run,])) * 1}}
    f[run]    <- fn(par[run,], ...)
    # Check convergence
    if (reltol(f[run], tol) > abs(f[run] - f[run-1])) {
      convergence = T
      if(verbose==T) {
        setTxtProgressBar(pb, maxit)
        cat("\nConvergence achieved!")
      }
      break
    } else {
        if(verbose==T) {
          setTxtProgressBar(pb, run)
        }
      }
  }
  if ({convergence==F} && {verbose==T}) {
    cat("\nConvergence may not have been achieved!")
  }
  if (run < maxit) {
    f    <- f[-c({run+1}:maxit)]
    par  <- par[-c({run+1}:maxit),]
    df   <- df[-c({run+1}:maxit),]
    step <- step[-c({run+1}:maxit),]
  }
  if(verbose==T) {
    close(pb)
    cat("\n")
    # Final messages
    stopTime = proc.time()
    elapsedTime = stopTime - startTime
    cat("It took ",round(elapsedTime[3],2)," secs for the run to finish. \n",
        sep="")
  }

  # Return results
  Results <- list("Cost"=f, "Estimates"=par, "Gradients"=df,
                  "par"=par[which.min(f),], "value"=f[which.min(f)],
                  "steps"=step)
  return(Results)
}

### Steepest Hit-and-Run Gradient Descent====
HRGD <- function(fn, startvalue, gr, ..., Interval=1e-6, maxit=100) {
  # Opening message
  cat("Hit-and-Run Gradient Descent will run for ",
      maxit, " iterations at most.\n\n", sep="")
  startTime = proc.time()

  # Initial settings
  pb   <- txtProgressBar(min=0, max=maxit, style=3)
  par  <- df <- array(dim = c(maxit,length(startvalue)))
  f    <- vector("numeric", length=maxit)
  step <- vector("numeric", length=maxit)
  convergence = F

  # First step
  df[1,]   <- gr(fn, startvalue, ..., Interval=Interval)
  SS       <- suppressWarnings(unlist(optim(c(0,0), fn=steep, mod=fn,
                                            ..., df=df[1,], est=startvalue,
                                            method="Nelder-Mead")$par))
  step[1]   <- mean(SS)
  par[1,]   <- startvalue -
               {SS[1]*df[1,]*{(abs(df[1,]) >  sd(df[1,])) * 1}} -
               {SS[2]*df[1,]*{(abs(df[1,]) <= sd(df[1,])) * 1}}
  f[1]      <- fn(par[1,], ...)
  setTxtProgressBar(pb, 1)

  # Start estimation
  for(run in 2:maxit) {
    # Calculate gradient and estimate step parameter
    temp <- sHAR(par=par[run-1,], mod=fn, ...)
    df[run,] <- temp$d
    step[run] <- temp$u
    par[run,] <- temp$prop
    f[run]    <- fn(par[run,], ...)
    # Check convergence
    if ({step[run] - step[run-1]} == 0) {
      convergence = T
      setTxtProgressBar(pb, maxit)
      cat("\nConvergence achieved!")
      break
    } else { setTxtProgressBar(pb, run) }
  }
  if (convergence == F) {
    cat("\nConvergence may not have been achieved!")
  }
  if (run < maxit) {
    f    <- f[-c({run+1}:maxit)]
    par  <- par[-c({run+1}:maxit),]
    df   <- df[-c({run+1}:maxit),]
    step <- step[-c({run+1}:maxit)]
  }
  close(pb)
  cat("\n")

  # Final messages
  stopTime = proc.time()
  elapsedTime = stopTime - startTime
  cat("It took ",round(elapsedTime[3],2)," secs for the run to finish. \n", sep="")

  # Return results
  Results <- list("Cost"=f, "Estimates"=par, "Gradients"=df,
                  "par"=par[which.min(f),], "value"=f[which.min(f)], "steps"=step)
  return(Results)
}

