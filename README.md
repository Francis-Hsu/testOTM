# testOTM
This repository contains a package which implements multivariate ranks and quantiles, together with their application in statistical testing, that based upon semi-discrete optimal transport theory. The underlying geometric computation is done through [Geogram](http://alice.loria.fr/index.php/software/4-library/75-geogram.html).

## Installation
In R console, run `devtools::install_github("Francis-Hsu/testOTM")`.

## Example
```R
library(testOTM)

# generate some data points
X = cbind(c(0.8, 0.7, 0.6, 0.1, 0.2), 
          c(0.7, 0.2, 0.8, 0.7, 0.1))
          
# compute the optimal transport map from a 2D uniform measure to the data
X.OTM = OTM_2D(X) 

# plot the restricted Voronoi diagram and the restricted Delaunay triangulation 
plot(X.OTM, which = "Both", draw.map = T)
```

## Reference
Promit Ghosal, Bodhisattva Sen (2019). [*Multivariate Ranks and Quantiles using Optimal Transportation and Applications to Goodness-of-fit Testing*](https://arxiv.org/abs/1905.05340).

Bruno Lévy, Erica L. Schwindt (2018). [*Notions of optimal transport theory and how to implement them on a computer*](https://doi.org/10.1016/j.cag.2018.01.009).
