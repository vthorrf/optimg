# Wrapper
gradient <- function(par, fn, ..., Interval=1e-6, order=1) {
  FN <- function(par) fn(par, ...)
  return(gradN(f=FN, par=par, h=Interval, order=order))
}
