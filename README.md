
<!-- README.md is generated from README.Rmd. Please edit that file -->

# exdiva

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/almeidasilvaf/exdiva)](https://github.com/almeidasilvaf/exdiva/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/almeidasilvaf/exdiva)](https://github.com/almeidasilvaf/exdiva/pulls)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Codecov test
coverage](https://codecov.io/gh/almeidasilvaf/exdiva/branch/devel/graph/badge.svg)](https://app.codecov.io/gh/almeidasilvaf/exdiva?branch=devel)
<!-- badges: end -->

The goal of `exdiva` is to …

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `exdiva` from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("exdiva")
```

And the development version from
[GitHub](https://github.com/almeidasilvaf/exdiva) with:

``` r
BiocManager::install("almeidasilvaf/exdiva")
```

## Citation

Below is the citation output from using `citation('exdiva')` in R.
Please run this yourself to check for any updates on how to cite
**exdiva**.

``` r
print(citation('exdiva'), bibtex = TRUE)
```

Please note that the `exdiva` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `exdiva` project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

- Continuous code testing is possible thanks to [GitHub
  actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
  through *[usethis](https://CRAN.R-project.org/package=usethis)*,
  *[remotes](https://CRAN.R-project.org/package=remotes)*, and
  *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized
  to use [Bioconductor’s docker
  containers](https://www.bioconductor.org/help/docker/) and
  *[BiocCheck](https://bioconductor.org/packages/3.18/BiocCheck)*.
- Code coverage assessment is possible thanks to
  [codecov](https://codecov.io/gh) and
  *[covr](https://CRAN.R-project.org/package=covr)*.
- The [documentation website](http://almeidasilvaf.github.io/exdiva) is
  automatically updated thanks to
  *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
- The documentation is formatted thanks to
  *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.18/biocthis)*.
