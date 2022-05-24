# Wrapper
optimg <- function(par, fn, gr=NULL, ..., method=c("SGD","STGD","LMM","ADAM"),
                   Interval=1e-6, maxit=100, tol=1e-8, full=F, verbose=F) {
  ### Initial checks
  ## Method
  if(is.vector(method) | is.null(method)) {
    method = method[1]
  }
  if({length(par)==1} & {method=="STGD"}) {
    method="SGD"
  }
  ## Cost function
  FN <- function(par) fn(par, ...)
  ## Gradient function
  if (!is.null(gr)) {
    GR <- function(par) gr(par, ...)
  } else if(is.null(gr)) {
    GR <- function(par) grad(FN, par, Interval)
  } else stop("Specification of gradient function is incorrect.")
  
  ### Select method
  if(method=="SGD"){
    Result <- SGD(FN, par, GR, Interval, maxit, tol, verbose)
  } else if(method=="STGD"){
    Result <- STGD(FN, par, GR, Interval, maxit, tol, verbose)
  } else if(method=="LMM"){
    Result <- LMM(FN, par, GR, Interval, maxit, tol, verbose)
  } else if(method=="ADAM"){
    Result <- ADAM(FN, par, GR, Interval, maxit, tol, verbose)
  } else stop("Unkown method. Please check optimg's documentation.")
  
  ### Return results
  par <- as.matrix(Result$Estimates)[Result$MaxInt,]
  value <- Result$Cost[Result$MaxInt]
  counts <- Result$MaxInt
  convergence <- Result$convergence
  if(full==T) {
    return(list("par"=par, "value"=value, "counts"=counts,
                "convergence"=convergence, "Full"=Result))
  } else if(full==F) {
    Results <- list("par"=par, "value"=value, "counts"=counts, "convergence"=convergence)
    return(Results)
  } else stop("Incorrect value for argument 'full'.")
}
