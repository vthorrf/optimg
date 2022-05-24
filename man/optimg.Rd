\name{optimg}
\alias{optimg}
\title{General-Purpose Gradient-Based Optimization}
\description{
General-purpose optimization based on gradient descent algorithms. It is greatly inspired on the stats::optim function, aiming at increased usability and adaptability for optim's users.
}
\usage{
optimg(par, fn, gr=NULL, ..., method=c("SGD","STGD","LMM","ADAM"),
       Interval=1e-6, maxit=100, tol=1e-8, full=F, verbose=F)
}
\arguments{
   \item{par}{Initial values for the parameters to be optimized over.}
   \item{fn}{A function to be minimized. It should return a scalar result.}
   \item{gr}{A function to return the gradient. Should include "Interval" as an argument, even if not used. If it is NULL, a finite-difference approximation will be used.}
   \item{...}{Further arguments to be passed to fn and gr.}
   \item{method}{The method to be used. See ‘Details’.}
   \item{Interval}{The epsilon difference to be used by gr.}
   \item{maxit}{The maximum number of iterations. Defaults to 100.}
   \item{tol}{Relative convergence tolerance. The algorithm stops if it is unable to reduce the value by a factor of reltol * (abs(val) + reltol) at a step. Defaults to 1e-8.}
   \item{full}{Boolean argument for returning all results, or only the basics (see Value). Defaults to FALSE.}
   \item{verbose}{Boolean argument for printing update status. Defaults to FALSE.}
}
\details{
As with the optim function, ... must be matched exactly. Also as with optim, optimg performs minimization. All the methods implemented are the "steepest" variation of the original methods. This means that tuning parameters are optimized at each step automatically. This makes the algorithms slower, but also more effective, especially for complex models.

The default method is the traditional Steepest Gradient Descent ("SGD") algorithm.

Method "STGD" is a variation of the Steepest Gradient Descent method which optimizes different step sizes for two groups of gradients: those within a standard deviation (below or above the mean), and those beyond a standard deviation (below or above the mean). If the specified function has only one parameter, SGD is used instead.

Method "LMM" is the Levenberg-Marquardt method. This method computes an approximate regularized truncated Hessian at each step. The LMM is a robust version of the Gaussian-Newton method, which means that it can better handle far off initial values.

Method "ADAM" is the Adaptive Moment Estimation method. This method computes adaptive learning rates for each parameter. Adam stores both exponentially decaying average of past squared gradients, as well as measures of momentum.
}
\value{
If full = FALSE, a list with components:
   \item{par}{The best set of parameters found.}
   \item{value}{The value of fn corresponding to par.}
   \item{counts}{A scalar giving the number of iterations before convergence is reached.}
   \item{convergence}{An integer code. 0 indicates successful completion. 1 indicates that the iteration limit maxit had been reached.}
If full = TRUE, apart from the above, components regarding the changes of par, gradient, value, and auxiliary parameters will also be returned.
}
\examples{
### Fit a simple linear regression with RMSE as cost function
require(optimg)

# Predictor
x <- seq(-3,3,len=100)
# Criterion
y <- rnorm(100, 2 + {1.2*x}, 1)

# RMSE cost function
fn <- function(par, Y, X) {
  mu <- par[1] + {par[2] * X}
  rmse <- sqrt(mean({Y-mu}^2))
  return(rmse)
}

# Compare optimization methods
optim(c(0,0), fn, Y=y, X=x, method="Nelder-Mead")
optim(c(0,0), fn, Y=y, X=x, method="BFGS")
optim(c(0,0), fn, Y=y, X=x, method="CG")
optim(c(0,0), fn, Y=y, X=x, method="L-BFGS-B")
optimg(c(0,0), fn, Y=y, X=x, method="SGD")
optimg(c(0,0), fn, Y=y, X=x, method="STGD")
optimg(c(0,0), fn, Y=y, X=x, method="LMM")
optimg(c(0,0), fn, Y=y, X=x, method="ADAM")
}
