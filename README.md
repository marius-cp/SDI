
# SDI (Score Decompostion Inference)

An R package to quantify the sampling uncertainty of Bregman type score
decomposition.

## Installation

Development version: the current version is available from GitHub.

``` r
# install.packages("devtools")
devtools::install_github("marius-cp/SDI")
```

## Example

- The following is based on the inflation application in the paper
  **Robust Mean Forecast Evaluation with Score Decompositions**; see
  [GitHub
  repository](https://github.com/marius-cp/InferenceScoreDecompositions).
  Link to the paper (to be included).

``` r
library(SDI)
pathdata <- "/Users/mp/Library/CloudStorage/Dropbox/study/PhD/Projects/RM_CORP/github/InferenceScoreDecompositions/02_application/inflation/murphy_replication/inflation_mean.rds"

dat <- readRDS(pathdata) %>% 
  dplyr::filter(stemp_rlz<="2020-12-01")
```

- Apply the SDI function:

``` r
sdi <- 
  SDI(
  X_1 = dat$michigan, 
  X_2 = dat$spf, 
  Y = dat$rlz, 
  S = function(x,y) (x-y)^2, 
  V = function(x,y) (x-y)*2, 
  Spp = function(x,y) 2
)
```

- Extract the score decomposition for the Michigan and SPF forecasts:

``` r
sdi$asyvardm$dec_1 # michigan
```

    ## # A tibble: 1 × 4
    ##       s   mcb    dsc   unc
    ##   <dbl> <dbl>  <dbl> <dbl>
    ## 1  2.00 0.455 0.0690  1.62

``` r
sdi$asyvardm$dec_2 #spf
```

    ## # A tibble: 1 × 4
    ##       s   mcb   dsc   unc
    ##   <dbl> <dbl> <dbl> <dbl>
    ## 1  1.46 0.361 0.516  1.62

- Extract the p-values of test for equal scores, equal miscalibration
  and equal discrimination:

``` r
sdi$asyvardm$pvals
```

    ##      s_pval  mcb_pval   dsc_pval
    ## 1 0.2835178 0.8105969 0.03703164
