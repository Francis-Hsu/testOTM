# testOTM
This repository contains a package which implements the multivariate ranks and quantiles, together with their application in statistical testing, that based on the semi-discrete optimal transportation. The underlying geometric computation is performed by the [Geogram](http://alice.loria.fr/index.php/software/4-library/75-geogram.html) library.

## Installation
In R, run `devtools::install_github("Francis-Hsu/testOTM")`.

## Example
```R
library(testOTM)

# generate some data
n = 100
X = cbind(rnorm(n, 0, 1), rnorm(n, 0, 0.5))

# scale the data to [0, 1] range
X = scaling.min.max(X)

# compute the optimal transport map from U[0, 1]^2 to the data
X.OTM = otm.fit(X)

# plot the restricted Voronoi diagram and the restricted Delaunay triangulation 
plot(X.OTM, which = "Both", draw.map = T)
```

## Reference
Promit Ghosal, Bodhisattva Sen (2019). [*Multivariate Ranks and Quantiles using Optimal Transportation and Applications to Goodness-of-fit Testing*](https://arxiv.org/abs/1905.05340).

Bruno Lévy, Erica L. Schwindt (2018). [*Notions of optimal transport theory and how to implement them on a computer*](https://doi.org/10.1016/j.cag.2018.01.009).
