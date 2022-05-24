\name{gradient}
\alias{gradient}
\title{Numerical Gradient of a Function}
\description{
Calculate the gradient of a function by finite difference coefficient.
}
\usage{
gradient(par, fn, ..., Interval=1e-6, order=1)
}
\arguments{
   \item{par}{Initial values for the parameters to be optimized over.}
   \item{fn}{A function to be minimized. It should return a scalar result.}
   \item{...}{Further arguments to be passed to fn.}
   \item{Interval}{The epsilon difference to be used by gr.}
   \item{order}{An integer expressing the order of the derivative.}
}
\details{
This function calculates a numerical approximation of the kth-order derivative of the supplied function at the points par. The arguments in ... are passed to fn, and it is assumed that fn returns a single value.
}
\value{
A real scalar or vector of the approximated gradient(s).
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

# Return gradients
gradient(c(0,0), fn, Y=y, X=x, order=1)
gradient(c(0,0), fn, Y=y, X=x, order=2)
gradient(c(0,0), fn, Y=y, X=x, order=3)
}
