#' @importFrom magrittr "%>%"
NULL

#' Score Decomposition Inference (SDI)
#'
#' This function provides a wrapper for calculating the variance of the decomposition components for two (rival) forecasts
#' based on Mincer-Zarnowitz regressions, along with performing null tests for miscalibration and discrimination.
#'
#' @param X_1 Numeric vector for the first forecast.
#' @param X_2 Numeric vector for the second forecast.
#' @param Y Numeric vector for the outcome variable.
#' @param S Function, the scoring function used to assess the forecasts.
#' @param V Function, the identification function, where V is the derivative of S.
#' @param Spp Function, second derivative of the scoring function, default is NA.
#' @param selectionvector A (optional) 1x4 matrix object specifying the selection vector. Default is a 1x4 matrix of ones.
#' @param vcov_estimator Function, variance-covariance estimator from the sandwich package.
#'
#' @return A list of class "SDI" containing several components:
#'   - `asyvardm`: a list with the variance decomposition results including `asy_vars`, `pvals`, `Lambda`, `Omega`, `eta`, `Gamma`, `dec_1`, and `dec_2`.
#'   - `mcbnulltest_X1`: results of the miscalibration null test for the first forecast.
#'   - `mcbnulltest_X2`: results of the miscalibration null test for the second forecast.
#'   - `mcbIUUI`: combined p-values for the miscalibration tests.
#'   - `dscnulltestX1`: placeholder for discrimination null test for the first forecast.
#'   - `dscnulltestX2`: placeholder for discrimination null test for the second forecast.
#'   - `dscIUUI`: placeholder for combined p-values for the discrimination test.
#'   - `mySDIvar`: the selected variance of the SDI.
#'   - `dm`: a data frame with the Diebold-Mariano test results including `dm_asy_var`, `dm_stat`, and `dm_pval`.
#'
#' @examples
#' # Define a simple scoring function
#' S <- function(y, x) { (y - x)^2 }
#' # Define its derivative(s)
#' V <- function(x, y) { 2 * (y - x) }
#' Spp <- function(x, y) { 2 }
#' # Generate sample data
#' Y <- rnorm(100)
#' X_1 <- Y + rnorm(100, sd = 0.5)
#' X_2 <- Y + rnorm(100, sd = 0.5)
#' # Calculate score decomposition inference
#' results <- SDI(X_1, X_2, Y, S, V, Spp)
#' print(results)
#'
#' @export

SDI <- function(
    X_1,
    X_2,
    Y,
    S, # scoring function
    V ,# identification function, V:= S'
    Spp = NA,  # S''
    selectionvector = matrix(c(1,1,1,1), nrow = 1, ncol = 4), # selection vector
    vcov_estimator = sandwich::vcovHAC
    ){

  tt <- nrow(Y)

  asyvardm <- asy_var_dm(
    X_1=X_1,
    X_2=X_2,
    Y=Y,
    S=S,
    V=V,
    Spp = Spp,
    vcov_estimator = vcov_estimator
  )

  mcbnulltest_X1 <- mcb_null_test(
    X = X_1,
    Y = Y,
    S=S,
    V=V,
    Spp=Spp,
    vcov_estimator = vcov_estimator
  )

  mcbnulltest_X2 <- mcb_null_test(
    X = X_2,
    Y = Y,
    S=S,
    V=V,
    Spp=Spp,
    vcov_estimator = vcov_estimator
  )

  mcbIUUI <- pmax(
    asyvardm$pvals$mcb_pval,
    2*pmin(
      mcbnulltest_X1$mcb_pval,
      mcbnulltest_X2$mcb_pval
    )
  )

  dscnulltest_X1 <- dsc_null_test(
    X = X_1,
    Y = Y,
    S=S,
    V=V,
    Spp=Spp,
    vcov_estimator = vcov_estimator
  )

  dscnulltest_X2 <- dsc_null_test(
    X = X_2,
    Y = Y,
    S=S,
    V=V,
    Spp=Spp,
    vcov_estimator = vcov_estimator
  )

  dscIUUI <- pmax(
    asyvardm$pvals$dsc_pval,
    2*pmin(
      dscnulltest_X1$dsc_pval,
      dscnulltest_X2$dsc_pval
    )
  )



  if(is.null(selectionvector)){
    mySDIvar <- NULL
  } else if (
    !is.matrix(selectionvector) | !((dim(selectionvector)[1]==1 ) & dim(selectionvector)[2] %in% c(1,4))
  ) {
    # i
    warning("The selection vector must be a 1x4 matrix object in R. We use NULL instead.")
    mySDIvar <- NULL
  } else if (dim(selectionvector)[1]==1 & dim(selectionvector)[2]==1) {
    mySDIvar <- as.numeric(selectionvector)*with(asyvardm,Gamma%*%Omega%*%t(Gamma))*as.numeric(selectionvector)
  } else {
    mySDIvar <- selectionvector%*%with(asyvardm,Gamma%*%Omega%*%t(Gamma))%*%t(selectionvector)
  }

  #original DM
  # sometimes original DM fails in my simulations.
  # If so, we are super conservative and set the test statistic to 0 (ie pval to 1).
  dm_asy_var <- tryCatch(
    {
      vcov_estimator(lm(as.numeric(S(x = X_1, y = Y) - S(x = X_2, y = Y)) ~ 1)) %>%
        as.numeric()
    },
    error = function(e) {999}
  )
  dm_stat <- mean(S(x=X_1,y=Y) - S(x=X_2,y=Y)) / sqrt(dm_asy_var) %>% as.numeric()
  dm_pval <- 2 * stats::pt(-abs(dm_stat), df = tt - 1) %>% as.numeric()
  dm <- tibble::tibble(
    dm_asy_var=dm_asy_var,
    dm_stat=dm_stat,
    dm_pval=dm_pval)


  out <-
    list(
      asyvardm=asyvardm,
      mcbnulltest_X1=mcbnulltest_X1,
      mcbnulltest_X2=mcbnulltest_X2,
      mcbIUUI=mcbIUUI,
      dscnulltest_X1=dscnulltest_X1,
      dscnulltest_X2=dscnulltest_X2,
      dscIUUI=dscIUUI,
      mySDIvar=mySDIvar,
      dm=dm
    )

  structure(out, class = "SDI")
}
