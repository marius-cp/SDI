#' @importFrom magrittr "%>%"
#' @importFrom tibble tibble
#' @importFrom sandwich NeweyWest vcovHAC
NULL

#' @title Mincer–Zarnowitz Quantile Regression Test
#'
#' @description
#' Performs Mincer–Zarnowitz tests based on quantile regression \eqn{y_t = \beta_0 + \beta_1 x_t + u_t}.
#' Computes Wald statistics to test for zero MCB and/or zero DSC components:
#' \itemize{
#'   \item Null hypothesis \strong{MCB = 0:} joint test \eqn{(\beta_0, \beta_1) = (0, 1)} for mean (or median/quantile) calibration.
#'   \item Null hypothesis \strong{DSC = 0:} test \eqn{\beta_1 = 0} for discrimination skill.
#' }
#' The test for \strong{MCB = 0} is VQR backtest by Gaglianone et al (2011).
#'
#' @param x Forecast values.
#' @param y Realized outcomes.
#' @param alpha Quantile level for the regression (e.g., 0.25).
#'
#' @return
#' A \code{list} with:
#' \itemize{
#'   \item \code{estimate} – summary of the quantile regression fit.
#'   \item \code{WaldStat.MCB}, \code{Wald.pval.MCB} – test for \eqn{(\beta_0, \beta_1) = (0, 1)}.
#'   \item \code{WaldStat.DSC}, \code{Wald.pval.DSC} – test for \eqn{\beta_1 = 0}.
#'   \item \code{fittedvalues} – tibble with \code{x} and fitted \code{quantile regression} values.
#' }
#'
#' @examples
#' set.seed(123)
#' y <- rnorm(100)
#' x <- y + rnorm(100)
#' MZ.test.quantile(x, y, alpha = 0.5)
#'
#' @seealso \code{\link[quantreg]{rq}} for quantile regression estimation.
#'
#' @references
#' Gaglianone, W. P., Lima, L. R., Linton, O., & Smith, D. R. (2011). Evaluating Value-at-Risk models via quantile regression. \emph{Journal of Business and Economic Statistics}, 29(1), 150–160.
#' \doi{10.1198/jbes.2010.08288}
#'
#' @export
MZ.test.quantile <- function(x, y, alpha){

  # Quantile regression
  MZ.QR.fit <- quantreg::rq(y ~ x, tau = alpha)
  Summary.MZ.QR.fit <- summary(MZ.QR.fit, covariance = TRUE, se = "nid")


  beta_hat <- coef(MZ.QR.fit)
  V_beta   <- Summary.MZ.QR.fit$cov

  ## --- Test for MCB = 0: H0: (beta0, beta1) = (0, 1) ---
  R_MCB <- diag(2)              # selects (beta0, beta1)
  r_MCB <- c(0, 1)
  diff_MCB <- R_MCB %*% beta_hat - r_MCB

  Wald_MCB <- as.numeric(t(diff_MCB) %*% solve(R_MCB %*% V_beta %*% R_MCB) %*% diff_MCB)
  pval_MCB <- 1 - pchisq(Wald_MCB, df = 2)


  ## --- Test for DSC = 0: H0: beta1 = 0 ---
  R_DSC <- c(0, 1)              # selects slope only
  r_DSC <- 0
  diff_DSC <- R_DSC %*% beta_hat - r_DSC

  Wald_DSC <- as.numeric(t(diff_DSC) %*% solve(R_DSC %*% V_beta %*% R_DSC) %*% diff_DSC)
  pval_DSC <- 1 - pchisq(Wald_DSC, df = 1)

  return(
    list(
      estimate        = Summary.MZ.QR.fit,
      WaldStat.MCB    = Wald_MCB,
      Wald.pval.MCB   = pval_MCB,
      WaldStat.DSC    = Wald_DSC,
      Wald.pval.DSC   = pval_DSC,
      fittedvalues    = tibble("x" = x, "fit" = MZ.QR.fit$fitted.values)
    )
  )
}


#' Decomposition of Quantile Scores
#'
#' This function performs a decomposition of the average score given forecasts,
#' recalibrated forecasts, and actual outcomes. It allows for flexible use of any quantile scoring function `S`
#' provided by the user, to calculate the mean score, marginal contribution bias (MCB),
#' and the decomposition of the scoring contribution (DSC).
#'
#' @param x Numeric vector of forecast values.
#' @param x_rc Numeric vector of recalibrated forecast values.
#' @param y Numeric vector of actual outcome values.
#' @param S Function, a user-defined generalized piecewise linear (GPL) scoring function that takes thre arguments (forecast, outcome, alpha)
#'   and returns a numeric value indicating the score.
#' @param alpha Quantile level.
#'
#' @return A tibble with columns:
#'   - `s`: the average score using the original forecast.
#'   - `mcb`: miscalibration component
#'   - `dsc`: discrimination component
#'   - `unc`: uncertainty component
#'
#' @examples
#' # tbc
#'
#' @export

qdecomposition <- function(
    x,
    x_rc,
    y,
    S,
    alpha
) {
  # inputs:
  #x forecast
  #x_rc recalibrated forecast
  #y realization
  # output: score and score decomposition

  s <- mean(S(x, y, alpha))
  #res <- y-x
  #c_rc_ucond <- optim(par = 0, fn = function(c) score(x + c, y),method = "Brent", lower = min(res), upper = max(res))$par
  #s_rc_ucond <- score(x + c_rc_ucond, y);s_rc_ucond
  s_rc <- mean(S(x_rc, y, alpha))
  s_rc
  s_mg <- mean(S(stats::quantile(y, probs=alpha), y, alpha))
  s_mg
  mcb <- s - s_rc
  #umcb <- s - s_rc_ucond
  #cmcb <- s_rc_ucond - s_rc
  dsc <- s_mg - s_rc
  dsc
  unc <- s_mg
  #(mcb-dsc+unc == s)
  return(tibble::tibble(s, mcb, dsc, unc))
}


#' Flexible MZ-style Quantile Regression Using Different GPL Scoring Functions
#'
#' This function fits a quantile MZ-style regression using different generalized piecewise linear (GPL) scoring functions of order alpha S(x, y, alpha), e.g., for the asymmetric piecewise linear scoring function \eqn{S(x,y) = (\mathbb{1}\{y < x\} - \alpha)(x-y)}.,
#' where the forecast is x and the realization is y. All scoring functions must be Bregman scores.
#' It returns the decomposition of the average score and the coefficients of the MZ regression under the given scoring function.
#'
#' @param x Numeric vector: Forecasts to be evaluated.
#' @param y Numeric vector: Realized values (must have the same length as x).
#' @param S function: Scoring function, which must be a function S(x, y, alpha).
#' @param alpha Quantile level, \eqn{\alpha \in (0,1)}.
#'
#' @return A list containing the following components:
#'   \item{x_rc}{Numeric vector: MZ recalibrated mean forecasts under the given scoring function.}
#'   \item{coef}{Numeric vector: Coefficients of the MZ regression under the given scoring function.}
#'   \item{dec}{tibble: Decomposition of the average score.}
#'
#' @details
#' The function uses optimization techniques to estimate the parameters of the mean forecast model based on the specified loss function.
#' Users should ensure that the scoring function adheres to the properties of Bregman scores.
#'
#' @examples
#' # tbc
#'
#' @seealso
#' \code{\link{optim}}, \code{\link{lm}}
#'
#' @keywords mean forecast loss function optimization
#' @export
#' @importFrom  magrittr `%>%`
#'
qMZdec <- function(x, y, alpha, S) {
  # loss for optim function
  loss <- function(params, x, y, S) {
    beta0 <- params[1]
    beta1 <- params[2]
    haty <- beta0 + beta1 * x
    error <- S(x=haty,y=y, alpha=alpha)
    return(sum(error))
  }

  # Initial guess for parameters
  initial_params <- c(0, 1)
  # Define lower bounds for parameters (beta0 >= 0)
  lower_bounds <- c(-Inf, -Inf)

  if((identical(body(S), body(function(x,y, alpha) (1 * (y < x) - alpha) * (x - y))))){
    mzreg <- quantreg::rq(y~x, tau=alpha)
    x_rc <- stats::fitted(mzreg)
    coef <- as.vector(stats::coefficients(mzreg))
  } else {
    result <- stats::optim(par = initial_params, fn = loss, x = x, y = y, alpha=alpha, S=S)
    # Optimal parameters
    optimal_beta0 <- result$par[1]
    optimal_beta1 <- result$par[2]

    # fitted values
    x_rc <- optimal_beta0 + optimal_beta1*x

    # coefficients
    coef <- c(optimal_beta0, optimal_beta1)
  }

  dec <- qdecomposition(x=x, y=y, x_rc=x_rc, S=S, alpha=alpha)

  return(list(x_rc=x_rc,coef=coef, dec = dec))

}



#' Calculate variance of the decomposition components for two (rival) quantile forecasts based on Mincer-Zarnowitz regressions.
#'
#' This function computes variance of the score decomposition for two rival quantile forecasts X_1 and X_2
#' based on the provided data. It calculates components related to the score (s), miscalibration (mcb), and discrimination (dsc).
#'
#' @param X_1 Numeric vector  for the first forecast.
#' @param X_2 Numeric vector for the second forecast.
#' @param Y Numeric vector for the outcome variable.
#' @param S Function, the scoring function used to assess the mean forecasts.
#' @param V Function, the identification function, where V is the derivative of S.
#' @param alpha Quantile level, \eqn{\alpha \in (0,1)}.
#' @param vcov_estimator Function, variance-covariance estimator from the sandwich package.
#' @param tol It can happen the determinat is not meaningfully calculated. We check for it with `det(Omega) <= 10^(-tol)`.
#' In this case the `Omega` output is a 5x5 matrix of `NaN` values.
#'
#' @return A list containing several components:
#'   - `asy_vars`: a data frame with the variances of s, mcb, and dsc for both forecasts.
#'   - `pvals`: p-values for testing the significance of differences.
#'   - `Lambda`: matrix used in the variance calculations.
#'   - `temp_Omega`: raw matrix of a_{1,t}, ..., e_t that is used to estimat cov matrix
#'   - `Omega`: covariance matrix based on the scoring function.
#'   - `eta`: transformation matrix used in variance calculations.
#'   - `Gamma`: product of eta and Lambda matrices.
#'   - `dec_1`: decomposition results for the first forecast.
#'   - `dec_2`: decomposition results for the second forecast.
#'
#' @details It can happen that covariance matrix estimation fails.
#' In this cases the `Omega` output is a 5x5 matrix of `NaN` values.
#' The p-values are set to 1 in this case to facilitate conservative test decisions.
#'
#' @importFrom magrittr `%>%`
#'
#' @examples
#' # tbc
#'
#' @export
q_asy_var_dm <- function(
    X_1,
    X_2,
    Y,
    S,
    V ,
    alpha,
    vcov_estimator = sandwich::vcovHAC,
    tol = 50
){

  tt <- length(Y)

  # fcast 1
  MZdec1 <- qMZdec(
    x=X_1,
    y=Y,
    S=S,
    alpha=alpha
  )
  mz_1 <- MZdec1$x_rc
  dec_1 <- MZdec1$dec
  s_1 <- dec_1$s
  mcb_1 <- dec_1$mcb
  dsc_1 <- dec_1$dsc
  unc_1 <- dec_1$unc

  # fcast 2
  MZdec2 <- qMZdec(
    x=X_2,
    y=Y,
    S=S,
    alpha = alpha
  )
  mz_2 <- MZdec2$x_rc
  dec_2 <- MZdec2$dec
  s_2 <- dec_2$s
  mcb_2 <- dec_2$mcb
  dsc_2 <- dec_2$dsc
  unc_2 <- dec_2$unc

  s_dif <- s_1-s_2
  mcb_dif <- mcb_1-mcb_2
  dsc_dif <- dsc_1-dsc_2


  u=g1_1=g1_2=g2_1=g2_2=0

  Lambda <- matrix(
    c(
      1, 0, 0, 0,         0, 0, 0, 0,           0, 0,
      0, 1, 0, 0,         0, 0, 0, 0,           0, 0,
      0, 0, g1_1, g1_2,   0, 0, 0, 0,           0, 0,
      0, 0, 0, 0,         1, 0, 0, 0,           0, 0,
      0, 0, 0, 0,         0, 1, 0, 0,           0, 0,
      0, 0, 0, 0,         0, 0, g2_1, g2_2,     0, 0,
      0, 0, 0, 0,         0, 0, 0, 0,           1, 0,
      0, 0, 0, 0,         0, 0, 0, 0,           0, u
    ),
    byrow = T,
    nrow = 8,
    ncol = 2*(2+2)+2
  )
  eta <- matrix(
    c(
      1,-1,1,    0,0,0,   0,0,
      0,-1,1,    0,0,0,   -1,1,
      0,0,0,     1,-1,1,  0,0,
      0,0,0,     0,-1,1,  -1,1
    ),
    byrow = T,
    nrow = 4,
    ncol = 8
  )

  Gamma <- eta%*%Lambda[,c(-3,-4,-7,-8,-10)]

  marg <- function(z) stats::quantile(z, probs=alpha)

  temp_Omega <-
    as.matrix(data.frame(
      a_1 = S(y = Y, x = X_1, alpha = alpha) ,
      c_1 = S(y = Y, x = (mz_1), alpha = alpha) ,
      #v1_1 = V(x = (mz_1), y = Y) ,
      #v2_1 = V(x = (mz_1), y = Y) * X_1,
      a_2 = S(y = Y, x = X_2, alpha = alpha) ,
      c_2 = S(y = Y, x = (mz_2), alpha = alpha) ,
      #v1_2 = V(x = (mz_2), y = Y) ,
      #v2_2 = V(x = (mz_2), y = Y) * X_2,
      e = S(x = marg(Y), y = Y, alpha = alpha)
      #u = V(x = marg(Y), y = Y)
    ))

  Omega <- tryCatch(vcov_estimator(lm(temp_Omega ~ 1)), error = function(e){FALSE})

  Omega <- if(identical(Omega, FALSE)) {
    matrix(NaN, nrow = 5, ncol = 5)
  } else if ( is.matrix(Omega) && det(Omega) <= 10^(-tol)){
    matrix(NaN, nrow = 5, ncol = 5)
  } else {
    Omega
  }



  omega1 <- t(c(1,-1,-1,1) )# for classical DM test
  omega2 <- t(c(1,0,-1,0)) # for MCB test
  omega3 <- t(c(0,1,0,-1)) # for DSC test



  # asym var score diff (equivalent to variance in the DM test)
  s_var <- omega1%*%Gamma%*%Omega%*%t(Gamma)%*%t(omega1)
  s_var

  # asym var mcb diff
  mcb_var <- omega2%*%Gamma%*%Omega%*%t(Gamma)%*%t(omega2)
  mcb_var

  # asym var dsc diff
  dsc_var <- omega3%*%Gamma%*%Omega%*%t(Gamma)%*%t(omega3)
  dsc_var

  tab <-
    data.frame(
      s_1, s_2,  s_var,
      mcb_1, mcb_2, mcb_var,
      dsc_1, dsc_2, dsc_var,
      unc_1, unc_2
    )

  pvals <- get_pval(dt=tab,tt=tt)
  is_nan <- sapply(pvals, is.nan)
  pvals[is_nan] <- 1

  out <-
    list(
      asy_vars=tab,
      pvals = pvals,
      Lambda = Lambda,
      eta = eta,
      Gamma=Gamma,
      temp_Omega = temp_Omega,
      Omega = Omega,
      dec_1=dec_1,
      dec_2=dec_2
    )

  return(out)


}



