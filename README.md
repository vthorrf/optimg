optimg: General-purpose Gradient-based Optimization (version 0.1.2)
=============
<p align="center">
    <a href="https://www.repostatus.org/#active"><img src="https://www.repostatus.org/badges/latest/active.svg" alt="Repository status"/></a>
    <a href="https://github.com/vthorrf/optimg/releases"><img src="https://img.shields.io/github/v/release/vthorrf/optimg?display_name=tag&sort=semver"/></a>
    <a href="https://www.r-pkg.org/pkg/optimg"><img src="https://www.r-pkg.org/badges/version/optimg" alt="CRAN version"/></a>
    <a href="https://www.r-pkg.org/pkg/optimg"><img src="https://cranlogs.r-pkg.org/badges/grand-total/optimg" alt="CRAN RStudio mirror downloads"/></a>
    <a href="https://cranchecks.info/pkgs/optimg"><img src="https://cranchecks.info/badges/worst/optimg" alt="CRAN checks"/></a>
</p>

This package is a general purpose tool for helping users to implement gradient descent methods for function optimization. Currently, the Steepest 2-Groups Gradient Descent and the Adaptive Moment Estimation (Adam) are the methods implemented. Other methods will be implemented in the future.

This package should be considered experimental in this point of development.

# Installation #

Using the 'remotes' package:

    install.packages("remotes")
    remotes::install_github("vthorrf/optimg")
