
<!-- README.md is generated from README.Rmd. Please edit that file -->

# testOTM

<!-- badges: start -->

<!-- badges: end -->

This repository contains a R package which implements the multivariate
ranks and quantiles, together with their application in statistical
testing, that based on the semi-discrete optimal transportation. The
underlying geometric computation is performed by the
[Geogram](http://alice.loria.fr/index.php/software/4-library/75-geogram.html)
library.

## Installation

<!--
You can install the released version of testOTM from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("testOTM")
```
-->

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Francis-Hsu/testOTM")
```

## Example

This is a basic example which shows you how to use `testOTM`:

``` r
library(testOTM)

# generate some data
n = 100
X = cbind(rnorm(n, 0, 1), rnorm(n, 0, 0.5))

# scale the data to [0, 1] range
X = scaling.min.max(X)

# compute the optimal transport map from U[0, 1]^2 to the data
X.OTM = otm.fit(X)

# plot the restricted Voronoi diagram and the restricted Delaunay triangulation
par(mfrow = c(1, 2))
plot(X.OTM, which = "Both", draw.center = F, draw.map = T)
```

![](man/figures/README-example-1.png)<!-- -->

## Reference

<div id="refs" class="references">

<div id="ref-GS2019">

Ghosal, Promit, and Bodhisattva Sen. 2019. “Multivariate Ranks and
Quantiles Using Optimal Transportation and Applications to
Goodness-of-Fit Testing.” <http://arxiv.org/abs/1905.05340>.

</div>

<div id="ref-LS2018">

Lévy, Bruno, and Erica L. Schwindt. 2018. “Notions of Optimal Transport
Theory and How to Implement Them on a Computer.” *Computers & Graphics*
72: 135–48. <https://doi.org/10.1016/j.cag.2018.01.009>.

</div>

</div>
