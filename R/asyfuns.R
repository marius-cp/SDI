#' @importFrom magrittr "%>%"
#' @importFrom tibble tibble
#' @importFrom sandwich NeweyWest vcovHAC
NULL

#' Calculate variance of the decomposition components for two (rival) mean forecasts based on Mincer-Zarnowitz regressions.
#'
#' This function computes the variance of the score decomposition for two rival mean forecasts X_1 and X_2 based on the provided data.
#' It calculates the score (s), as well as the components of the score decomposition, i.e., miscalibration (mcb), and discrimination (dsc), and uncertainty (unc).
#'
#' @param X_1 Numeric vector  for the first forecast.
#' @param X_2 Numeric vector for the second forecast.
#' @param Y Numeric vector for the outcome variable.
#' @param S Function for the scoring function used to assess the mean forecasts.
#' @param V Function for the identification function, where V is the derivative of S.
#' @param Spp Function for second derivative of the scoring function, default is NA.
#' @param vcov_estimator Function for variance-covariance estimator from the sandwich package.
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
#' # Define a simple scoring function
#' S <- function(y, x) { (y - x)^2 }
#' # Define its derivative
#' V <- function(x, y) { 2 * (y - x) }
#' Spp = function(x, y) 2
#' # Generate sample data
#' Y <- rnorm(100)
#' X_1 <- Y + rnorm(100, sd = 0.5)
#' X_2 <- Y + rnorm(100, sd = 0.5)
#' # Calculate asymmetric variance
#' results <- asy_var_dm(X_1, X_2, Y, S, V, Spp)
#' print(results)
#'
#' @export


asy_var_dm <- function(
    X_1,
    X_2,
    Y,
    S,
    V ,
    Spp,
    vcov_estimator = sandwich::vcovHAC,
    tol = 50
    ){

  tt <- length(Y)

  # fcast 1
  MZdec1 <- MZdec(
    x=X_1,
    y=Y,
    S=S
  )
  mz_1 <- MZdec1$x_rc
  dec_1 <- MZdec1$dec
  s_1 <- dec_1$s
  mcb_1 <- dec_1$mcb
  dsc_1 <- dec_1$dsc
  unc_1 <- dec_1$unc

  # fcast 2
  MZdec2 <- MZdec(
    x=X_2,
    y=Y,
    S=S
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

  # g(\hat\theta_{1,T}) first and second component (g1_1,g2_1)
  #g1_1 <- as.numeric(LLN_B(X=X_1,Y=Y)[1]) # --> to be made flexible with SPP!!!
  #g2_1 <- as.numeric(LLN_B(X=X_1,Y=Y)[2]) # --> to be made flexible with SPP!!!

  # g(\hat\theta_{2,T}) first and second component (g1_2,g2_2)
  #g1_2 <- as.numeric(LLN_B(X=X_2,Y=Y)[1]) # --> to be made flexible with SPP!!!
  #g2_2 <- as.numeric(LLN_B(X=X_2,Y=Y)[2]) # --> to be made flexible with SPP!!!

  # g(r) which is independent of forecast
  # u <-  as.numeric(LLN_E(Y=Y)) # --> to be made flexible with SPP!!!

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

  marg <- function(z) mean(z)

  temp_Omega <-
    as.matrix(data.frame(
      a_1 = S(y = Y, x = X_1) ,
      c_1 = S(y = Y, x = (mz_1)) ,
      #v1_1 = V(x = (mz_1), y = Y) ,
      #v2_1 = V(x = (mz_1), y = Y) * X_1,
      a_2 = S(y = Y, x = X_2) ,
      c_2 = S(y = Y, x = (mz_2)) ,
      #v1_2 = V(x = (mz_2), y = Y) ,
      #v2_2 = V(x = (mz_2), y = Y) * X_2,
      e = S(x = marg(Y), y = Y)
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



#' Test for zero Miscalibration (MCB=0)
#'
#' This function implements tests for the null hypothesis of zero miscalibration (MCB=0)
#' using forecast (X) and realization (Y) data. It employs two methods for testing:
#' a generalized Chi-square test of the loss function decomposition as proposed in the accombined paper
#' and the classical Mincer-Zarnowitz test.
#'
#' @param X A numeric vector representing the forecast values.
#' @param Y A numeric vector representing the actual realization values.
#' @param S scoring function
#' @param V identification function, first derivative of score V:= S'
#' @param Spp  second derviative of score S''
#' @param vcov_estimator Function, variance-covariance estimator from the sandwich package.
#' @return A tibble containing the results of the mis-calibration tests, including
#'         the MCB value, MCB variance, MCB p-value, Mincer-Zarnowitz p-value,
#'         and F-test p-value.
#'
#' @details The function first regresses Y on X to extract fitted values and residuals,
#'          then calculates the variance-covariance matrix using a HAC estimator.
#'          It proceeds to compute the MCB value, its variance, and p-value using the
#'          generalized Chi-square asymptotics. The Mincer-Zarnowitz test is conducted
#'          as a comparative methodology. The F-test p-value is also computed for
#'          assessing the joint hypothesis that the intercept is zero and the slope is one
#'          in the linear regression of Y on X.
#'
#'          Note: The function requires the `sandwich`, `MASS`, `car`, and `CompQuadForm`
#'          packages for various calculations.
#'
#' @examples
#' # Example usage:
#' S <- function(x, y) (x - y)^2
#' V <- function(x, y) 2*(x - y)
#' Spp <- function(x, y) 2
#' set.seed(1234)
#' X <- rnorm(100)
#' Y <- X + rnorm(100)
#' test_results <- mcb_null_test(X, Y, S, V, Spp)
#' print(test_results)
#'
#' @export

mcb_null_test <-
  function(
    X, #forecast
    Y, # realization,
    S,
    V ,# identification function, V:= S'
    Spp, # S''
    vcov_estimator = sandwich::vcovHAC
  ){

    if(!is.function(S)) {stop("Please provide a function for scoreing function 'S'.")}
    if(!is.function(V)) {stop("Please provide a function for identification function 'V'.")}
    if(!is.function(Spp)) {stop("Please provide a function for S''.")}


    tt <- length(X)

    # todo: flexible MZ reg functions!

    MZ <- MZdec(
      x = X,
      y = Y,
      S = S
    )

    dec <- decomposition(
      x=X,
      x_rc=MZ$x_rc,
      y=Y,
      S=S
    )



    MZreg <- lm(Y~X)
    W <- t(t(stats::model.matrix(MZreg))) # use lm to quickly generate W=(1,X)'

    Sigma <-tt*vcov_estimator(
      lm(as.vector(V(x=MZ$x_rc,y=Y))*W ~ 1)
    )
    Sigma

    Omega<- matrix(0, nrow=2,ncol=2)
    for (i in 1:nrow(W)) {
      temp <- Spp(X[i],Y[i])*(W[i,])%*%t(W[i,])
      Omega <- Omega+temp
    }
    Omega <- Omega/tt

    N <- t(t(MASS::mvrnorm(n=1,mu=c(0,0),Sigma=Sigma)))

    # F test
    freg <- lm((Y-X) ~ X)
    fsum <- summary(freg)
    F_pval <-
      car::linearHypothesis(
        freg,
        c("(Intercept)=0","X=0")
        #.vcov = sandwich::NeweyWest(reg, lag = 3)
      ) $`Pr(>F)`[2]

    sol <-
      tibble::tibble(
        mcb = dec$mcb,
        mcb_var =  as.numeric(t(N)%*%solve(Omega)%*%N),
        mcb_pval = CompQuadForm::imhof(
          2*tt*dec$mcb,
          lambda = eigen(solve(Omega)%*%Sigma)$values
        )$Qq %>% round(4),
        mz_pval = MZ.test.mean(haty = X, y=Y)$Wald.pval,
        F_pval = F_pval#pf(fsum$fstatistic[1], 1, tt-2, lower.tail = FALSE)
      )
    return(sol)
  }


#' Test for zero discrimination (DSC=0)
#'
#' This function implements tests for the null hypothesis of DSC=0, using
#' forecast (X) and realization (Y) data. It employs the generalized Chi-square
#' test of the score decomposition. Also the Wald and F test pendants are
#' reported.
#' @param X A numeric vector representing the forecast values.
#' @param Y A numeric vector representing the actual realization values.
#' @param S Scoring function.
#' @param V Identification function, first derivative of the score (V := S').
#' @param Spp Second derivative of the score (S'').
#' @param vcov_estimator Function, variance-covariance estimator from the sandwich package.
#' @return A tibble containing the results of the dynamic stochastic calibration tests, including
#'         the DSC value, DSC variance, DSC p-value, Mincer-Zarnowitz p-value, and Wald test p-value.
#'
#' @details The function first regresses Y on X to extract fitted values and residuals,
#'          then calculates the variance-covariance matrix using a HAC estimator.
#'          It proceeds to compute the DSC value, its variance, and p-value using the
#'          generalized Chi-square asymptotics.
#'
#'          Note: The function requires the `sandwich`, `MASS`, `car`, and `CompQuadForm`
#'          packages for various calculations.
#'
#' @examples
#' # Example usage:
#' S <- function(x, y) (x - y)^2
#' V <- function(x, y) 2*(x - y)
#' Spp <- function(x, y) 2
#' set.seed(1234)
#' X <- rnorm(100)
#' Y <- X + rnorm(100)
#' test_results <- dsc_null_test(X, Y, S, V, Spp)
#' print(test_results)
#'
#' @export
dsc_null_test <-
  function(
    X, #forecast
    Y, # realization,
    S,
    V ,# identification function, V:= S'
    Spp, # S''
    vcov_estimator = sandwich::vcovHAC
  ){

    # Check if S, V, and Spp are functions
    if(!is.function(S)) {stop("Please provide a function for scoring function 'S'.")}
    if(!is.function(V)) {stop("Please provide a function for identification function 'V'.")}
    if(!is.function(Spp)) {stop("Please provide a function for S''.")}

    # Number of observations
    tt <- length(X)

    # Fit a linear model for Y on X
    MZ <- MZdec(
      x = X,
      y = Y,
      S = S
    )

    # Perform the score decomposition
    dec <- decomposition(
      x = X,
      x_rc = MZ$x_rc,
      y = Y,
      S = S
    )

    # Get the model matrix and transpose it
    MZreg <- lm(Y ~ X)
    W <- t(stats::model.matrix(MZreg)) %>% t()

    # Compute the Wald p-value
    wald_pval <-
      car::linearHypothesis(
        MZreg,
        c("X=0"),
        vcov. = vcov_estimator(MZreg)
      ) $`Pr(>F)`[2]

    # Compute the F-test p-value
    model_summary <- summary(MZreg)
    f_statistic <- model_summary$fstatistic[1]
    f_p_value <- stats::pf(f_statistic, model_summary$fstatistic[2], model_summary$fstatistic[3], lower.tail = FALSE)

    # Initialize H_i_T matrix
    H_i_T <- matrix(0, nrow = 2, ncol = 2)
    for (i in 1:nrow(W)) {
      temp <- Spp(X[i], Y[i]) * (W[i, ]) %*% t(W[i, ])
      H_i_T <- H_i_T + temp
    }
    H_i_T <- H_i_T / tt

    # Initialize H_T matrix
    H_T <- matrix(0, nrow = 1, ncol = 1)
    for (i in 1:nrow(W)) {
      temp <- Spp(X[i], Y[i]) # Use reference forecast instead? should not matter under H_0.
      H_T <- H_T + temp
    }
    H_T <- H_T / tt

    # Initialize score contributions matrix
    score_contributions <- matrix(0, nrow = tt, ncol = 2 + 1)
    for (t in 1:tt) {
      r_t <- mean(Y)
      theta_i_t <- W[t, ]
      score_contributions[t, 1] <- V(r_t, Y[t])  # reference forecast
      score_contributions[t, 2:(2 + 1)] <- W[t, ] * ( V(x= MZ$x_rc[t], y = Y[t]) )  # recalibrated forecast
    }

    # Compute the variance-covariance matrix using the HAC estimator
    Pi_T <- tt * vcov_estimator(lm(score_contributions[, 2:3] ~ 1)) # only use relevant covariances (omit ref forecast as it drops out)

    # Generate multivariate normal random variables
    N <- MASS::mvrnorm(n = 1, mu = c(0, 0), Sigma = Pi_T) %>% t()

    # Compute the p-value using the Imhof method for quadratic forms
    dsc_pval <- CompQuadForm::imhof(
      2 * (tt) * dec$dsc,
      lambda = eigen(((solve(H_i_T) - (matrix(c(solve(H_T), 0, 0, 0), nrow = 2))) %*% Pi_T))$values
    )$Qq

    # Compute the DSC variance
    dsc_var = (N) %*% (solve(H_i_T) - (matrix(c(solve(H_T), 0, 0, 0), nrow = 2))) %*% t(N) %>% as.numeric()

    # Create a tibble with the results
    sol <-
      tibble(
        dsc = dec$dsc,
        dsc_var = dsc_var,
        dsc_pval = dsc_pval,
        f_pval = f_p_value,
        wald_pval = wald_pval
      )

    # Return the results
    return(sol)
  }




#' Flexible MZ Mean Regression Using Different Bregman Scoring Functions
#'
#' This function fits a mean MZ regression using different Bregman scoring functions S(x, y), e.g., for the squared error \eqn{S(x, y) = (x - y)^2},
#' where the forecast is x and the realization is y. All scoring functions must be Bregman scores.
#' It returns the decomposition of the average score and the coefficients of the MZ regression under the given scoring function.
#'
#' @param x Numeric vector: Forecasts to be evaluated.
#' @param y Numeric vector: Realized values (must have the same length as x).
#' @param S function: Scoring function, which must be a function S(x, y).
#'
#' @return A list containing the following components:
#'   \item{x_rc}{Numeric vector: MZ recalibrated mean forecasts under the given scoring function.}
#'   \item{coef}{Numeric vector: Coefficients of the MZ regression under the given scoring function.}
#'   \item{dec}{tibble: Decomposition of the average score.}
#'
#' @details
#' Users should ensure that the scoring function adheres to the properties of Bregman scores.
#'
#' @examples
#' x <- c(1, 2, 3, 4, 5)
#' y <- c(1.1, 2.2, 3.2, 4.1, 5.3)
#' result <- MZdec(
#'   x = x,
#'   y = y,
#'   S = function(x, y) (x - y)^2
#' )
#'
#' @seealso
#' \code{\link{optim}}, \code{\link{lm}}
#'
#' @export
#' @importFrom  magrittr `%>%`
#'
MZdec <- function(x, y, S) {
  # loss for optim function
  loss <- function(params, x, y, S) {
    beta0 <- params[1]
    beta1 <- params[2]
    haty <- beta0 + beta1 * x
    error <- S(x=haty,y=y)
    return(sum(error))
  }

  # Initial guess for parameters
  initial_params <- c(0, 1)
  # Define lower bounds for parameters (beta0 >= 0)
  lower_bounds <- c(-Inf, -Inf)

  if((identical(body(S), body(function(x, y) (x - y)^2)))){
    mzreg <- stats::lm(y ~ x)
    x_rc <- stats::fitted(mzreg)
    coef <- as.vector(stats::coefficients(mzreg))

    # Minimize squared error using optim function
    result <- stats::optim(par = initial_params, fn = loss, x = x, y = y, S=S)

    # Optimal parameters
    mse_optimal_beta0 <- result$par[1]
    mse_optimal_beta1 <- result$par[2]

    # print warning if closed form (OLS) and optin deviate
    if(
      !isTRUE(
        all.equal(
          round(c(mse_optimal_beta0,mse_optimal_beta1),3),
          round(as.vector(stats::coefficients(mzreg)),3)
        )
      )
    ){
      print(
        paste(
          "optim() and lm() return different coefficients:",
          round(mse_optimal_beta0, 3), ",", round(mse_optimal_beta1, 3),
          "vs",
          round(stats::coefficients(mzreg)[1], 3),",",round(stats::coefficients(mzreg)[2], 3)
        )
      )
      #print(c("optim() and lm() return different coefficent"), round(c(mse_optimal_beta0,mse_optimal_beta1),3), "and", round(as.vector(coefficients(mzreg)),3))
    }

  } else {
    result <- stats::optim(par = initial_params, fn = loss, x = x, y = y, S=S)
    # Optimal parameters
    optimal_beta0 <- result$par[1]
    optimal_beta1 <- result$par[2]

    # fitted values
    x_rc <- optimal_beta0 + optimal_beta1*x

    # coefficients
    coef <- c(optimal_beta0, optimal_beta1)
  }

  # only important for volatility !
  # qlike will not work in this case !
   if(any(x_rc<0)){print("!!! Reklaibration returns negative value !!! ")}

  dec <- decomposition(x=x, y=y, x_rc=x_rc, S=S)

  return(list(x_rc=x_rc,coef=coef, dec = dec))

}



#' Decomposition of the Bregman Score
#'
#' This function performs a decomposition of any Bregman scoring function given forecasts,
#' recalibrated forecasts, and actual outcomes. It allows for flexible use of any Bregman scoring function `S`
#' provided by the user, to calculate the average score, the miscalibration component (MCB), the
#' discrimination component (DSC), and the uncertainty component.
#'
#' @param x Numeric vector of forecast values.
#' @param x_rc Numeric vector of recalibrated forecast values.
#' @param y Numeric vector of actual outcome values.
#' @param S Function, a user-defined scoring function that takes two arguments (forecast, outcome)
#'   and returns a numeric value indicating the score.
#'
#' @return A tibble with columns:
#'   - `s`: the average score using the original forecast.
#'   - `mcb`: miscalibration component
#'   - `dsc`: discrimination component
#'   - `unc`: uncertainty component
#'
#' @examples
#' # Define a simple scoring function
#' S <- function(x, y) { (y - x)^2 }
#' # Generate sample data
#' y <- rnorm(100)
#' x <- y + rnorm(100, sd = 0.5)
#' x_rc <- y + rnorm(100, sd = 0.2)
#' # Decompose the scoring
#' decomposition_results <- decomposition(x, x_rc, y, S)
#' print(decomposition_results)
#'
#' @export

decomposition <- function(
    x,
    x_rc,
    y,
    S
) {
  # inputs:
  #x forecast
  #x_rc recalibrated forecast
  #y realization
  # output: score and score decomposition

  s <- mean(S(x, y))
  #res <- y-x
  #c_rc_ucond <- optim(par = 0, fn = function(c) score(x + c, y),method = "Brent", lower = min(res), upper = max(res))$par
  #s_rc_ucond <- score(x + c_rc_ucond, y);s_rc_ucond
  s_rc <- mean(S(x_rc, y))
  s_rc
  s_mg <- mean(S(mean(y), y))
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


#' Calculate p-values based on specified theory.
#'
#' This function calculates p-values based on the provided theory for different components.
#' It takes an asy_vars object containing relevant variables and the number of observations.
#'
#' @param dt An asy_vars object from function .
#' @param tt Number of observations.
#'
#' @return A data frame containing p-values for different components: s_pval, s_pvalDM, mcb_pval, and dsc_pval.
#'
#' @details
#' This function computes p-values for s, mcb, and dsc components based on the specified theory.
#' It also calculates the original DM p-value for the s component.
#'
#' @export

get_pval <- function(
    dt,   # must be asy_vars object
    tt # number ov obs
){
  # pvals based on our theory
  s_stat <- (dt$s_1 - dt$s_2) / sqrt(dt$s_var)
  s_pval <- 2 * stats::pt(-abs(s_stat), df = tt - 1)

  mcb_stat <- (dt$mcb_1 - dt$mcb_2) / sqrt(dt$mcb_var)
  mcb_pval <- 2 * stats::pt(-abs(mcb_stat), df = tt - 1)

  dsc_stat <- (dt$dsc_1 - dt$dsc_2) / sqrt(dt$dsc_var)
  dsc_pval <- 2 * stats::pt(-abs(dsc_stat), df = tt - 1)

  # original DM p-value
  #s_dm_orig_stat <- (dt$s_1 - dt$s_2) / sqrt(dt$s_varDM)
  #s_pvalDM <- 2 * pt(-abs(s_dm_orig_stat), df = tt - 1)

  out <-
    data.frame(
      s_pval, #s_pvalDM,
      mcb_pval, dsc_pval
    )

  return(out)
}


#' Mincer-Zarnowitz Mean Test
#'
#' This function performs the classical Mincer-Zarnowitz test to evaluate the accuracy of economic forecasts.
#'
#' @param haty Numeric vector of forecasted values.
#' @param y Numeric vector of actual values.
#' @param vcov_estimator Function, variance-covariance estimator from the sandwich package.
#'
#' @return A list containing several components:
#'   - `estimate`: Estimated coefficients from the Mincer-Zarnowitz regression.
#'   - `CI`: Confidence intervals for the estimated coefficients.
#'   - `WaldStat`: Wald statistic for testing the joint hypothesis that the intercept is zero and the slope is one.
#'   - `Wald.pval`: p-value associated with the Wald statistic.
#'   - `fittedvalues`: A tibble with the forecasted values and the corresponding fitted values from the regression.
#'
#' @references Mincer, J. and Zarnowitz, V. (1969). The evaluation of economic forecasts. In Economic Forecasts and Expectations: Analysis of Forecasting Behavior and Performance, pages 3–46. National Bureau of Economic Research, Inc.
#'
#' @importFrom magrittr `%>%`
#' @importFrom tibble tibble
#' @importFrom lmtest coefci
#' @importFrom sandwich NeweyWest vcovHAC
#' @importFrom stats lm coef pchisq
#'
#' @examples
#' # Generate sample data
#' y <- rnorm(100)
#' haty <- y + rnorm(100, sd = 0.5)
#' # Perform Mincer-Zarnowitz test
#' results <- MZ.test.mean(haty, y)
#' print(results)
#'
#' @export
MZ.test.mean <- function(
    haty,
    y,
    vcov_estimator = sandwich::vcovHAC
    ){
  MZ.fit <- stats::lm(y~haty)
  coefci <- lmtest::coefci(MZ.fit, vcov. = vcov_estimator) %>% round(3)
  # below solve() returns the inverse of the covar matrix
  W <- t(stats::coef(MZ.fit) - c(0,1)) %*% solve(vcov_estimator(MZ.fit)) %*% (stats::coef(MZ.fit) - c(0,1))
  pval <- 1 - stats::pchisq(W,2) %>% round(3)
  return(list(estimate=stats::coef(MZ.fit),
              CI=coefci,
              WaldStat=as.numeric(W),
              Wald.pval=as.numeric(pval),
              fittedvalues = tibble::tibble("x"=haty ,"fit" = MZ.fit$fitted.values)))
}
