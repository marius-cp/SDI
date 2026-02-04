#' @importFrom magrittr "%>%"
#' @importFrom tibble tibble
#' @importFrom mvtnorm pmvnorm
NULL


#' Galvão–Yoon (2024) HAC helper functions for time-series quantile regression
#'
#' @description
#' Internal helper functions used to compute HAC covariance matrices for quantile regression
#' following Galvao and Yoon (2024).  The functions below are adapted from the supplementary material / code of
#'  Galvao, A. F., & Yoon, J. (2024), JASA 119(547), 2305–2316.
#'  They are included here for internal use only.
#'
#' @details
#' Source: Galvao, A. F., & Yoon, J. (2024). HAC Covariance Matrix Estimation in Quantile Regression.
#' \emph{Journal of the American Statistical Association}, 119(547), 2305–2316.
#' \doi{10.1080/01621459.2023.2257365}
#'
#' @keywords internal
#' @name hac_qr_galvao_yoon_internal

NULL

ar1 <- function(yy,nt,ng){
  s1 <- 0; s2 <- 0
  for(i in 1:ng){
    sel = ((i-1)*nt+1):(i*nt)
    yh  = yy[sel]
    yh1 = yh[1:(nt-1)] - mean(yh[1:(nt-1)])
    s1 = s1 + sum(yh1*(yh[2:nt]-mean(yh[2:nt])))
    s2 = s2 + sum(yh1*yh1)
  }
  delta = (s1/s2)*nt/(nt-1) + 1/(nt-1)
  er = NULL
  for(i in 1:ng){
    sel = ((i-1)*nt+1):(i*nt)
    yh = yy[sel]
    eh = yh[2:nt] - mean(yh[2:nt]) - delta*yh[1:(nt-1)]
    er = c(er,eh)
  }
  sig = stats::sd(er)
  return(list(delta = delta, sigma = sig))
}

band <- function(yy,nt,ng,dy,ker){
  # bandwidth estimation. use sample correlation coefficients
  if(dy==1){
    delta = ar1(yy,nt,ng)$delta
    if(ker %in% c(1,2)){
      ah = ((2*delta)/(1-(delta^2)))^2
      band = round((1.1447*((ah*nt*ng)^{1/3})),0)
    }
    if(ker %in% c(3)){
      ah = ((2*delta)/((1-delta)^{2}))^2
      band = round((1.3221*((ah*nt*ng)^{1/5})),0)
    }
  }
  if(dy>1){
    delta = array(0,c(dy,1))
    sig = array(0,c(dy,1))
    for(j in 1:dy){
      aa = ar1(yy[,j],nt,ng)
      delta[j] = aa$delta
      sig[j] = aa$sigma
    }
    if(ker %in% c(1,2)){
      e1 = 4*(delta^2)*(sig^4)/(((1-delta)^6)*((1+delta)^2))
      e2 = (sig^4)/((1-delta)^4)
      ah = sum(e1)/sum(e2)
      band = round((1.1447*((ah*nt*ng)^{1/3})),0)
    }
    if(ker %in% 3){
      e1 = 4*(delta^2)*(sig^4)/((1-delta)^8)
      e2 = (sig^4)/((1-delta)^4)
      ah = sum(e1)/sum(e2)
      band = round((1.3221*((ah*nt*ng)^{1/5})),0)
    }
  }
  band = pmax(band,1)
  return(list(delta = delta, a.hat=ah, band=band))
}

band.new <- function(yy,hh,xx,tau,nt,ng,dy,ker,reg=0){
  # alternative bandwidth estimation method
  del = ar1(yy,nt,ng)$delta
  c1 = mvtnorm::pmvnorm(lower=-Inf,upper=c(0,0), mean=rep(-stats::qnorm(tau),2), sigma=matrix(c(1,del,del,1),nrow=2)) - tau^2
  if(abs(c1[1]) < 1e-8){c11 = 0}
  if(abs(c1[1]) >= 1e-8){c11 = c1[1]}
  del.h = c11/(tau*(1-tau))

  if(dy==1){
    if(reg==0){del.x = 1}
    if(reg==1){del.x = ar1(xx,nt,ng)$delta}
    delta = del.h*del.x
    if(ker %in% c(1,2)){
      ah = ((2*delta)/(1-(delta^2)))^2
      band = round((1.1447*((ah*nt*ng)^{1/3})),0)
    }
    if(ker %in% c(3)){
      ah = ((2*delta)/((1-delta)^{2}))^2
      band = round((1.3221*((ah*nt*ng)^{1/5})),0)
    }
  }
  if(dy>1){
    delta = array(0,c(dy,1))
    sig = array(0,c(dy,1))
    for(j in 1:dy){
      xh = xx[,j]
      if(prod(xh[1]==xh)==1){del.x = 1}
      if(prod(xh[1]==xh)!=1){del.x = ar1(xh,nt,ng)$delta}
      delta[j] = del.h*del.x
      sig[j] = ar1((hh*xh),nt,ng)$sigma
    }
    if(ker %in% c(1,2)){
      e1 = 4*(delta^2)*(sig^4)/(((1-delta)^6)*((1+delta)^2))
      e2 = (sig^4)/((1-delta)^4)
      ah = sum(e1)/sum(e2)
      band = round((1.1447*((ah*nt*ng)^{1/3})),0)
    }
    if(ker %in% 3){
      e1 = 4*(delta^2)*(sig^4)/((1-delta)^8)
      e2 = (sig^4)/((1-delta)^4)
      ah = sum(e1)/sum(e2)
      band = round((1.3221*((ah*nt*ng)^{1/5})),0)
    }
  }
  band = pmax(band,1)
  return(list(delta = delta, a.hat=ah, band=band))
}

step.long <- function(nt,ng,hx,mt,ker){
  # long-run variance estimation
  # ker = 1,2,3 means truncated, Bartlett, and QS kernels
  Jm = hx[,,1]
  if(ker==1){ # truncated kernel
    for(s in 1:mt){
      Jm = Jm + 2*((nt-s)/nt)*hx[,,(s+1)]
    }
  }
  if(ker==2){ # Bartlett kernel
    for(s in 1:mt){
      Jm = Jm + 2*(1-abs(s/mt))*((nt-s)/nt)*hx[,,(s+1)]
    }
  }
  if(ker==3){ # QS kernel
    for(s in 1:(nt-2)){
      z = s/mt
      Jm = Jm + 2*(25/(12*(pi^2)*(z^2)))*(sin(6*pi*z/5)/(6*pi*z/5)-cos(6*pi*z/5))*((nt-s)/nt)*hx[,,(s+1)]
    }
  }
  return(list(J.n = hx[,,1], J.r = Jm))
}

hac.se = function(x,e0,e1,uhat,tau,nt,mt,ker){
  # calculate HAC estimator
  fxxinv = e1*nt
  q = ncol(x)
  xx1 = array(0,c(q,q,(nt-1)))
  hx1 = array(0,c(q,q,(nt-1)))
  for(s in 0:(nt-2)){
    xh1 = x[1:(nt-s),]
    xh2 = x[(1+s):nt,]
    uh  = tau - (uhat < 0)
    xu1 = xh1*uh[1:(nt-s)]
    xu2 = xh2*uh[(1+s):nt]
    xx1[,,(s+1)] = crossprod(xh1,xh2)/(nt-s)
    hx1[,,(s+1)] = crossprod(xu1,xu2)/(nt-s)
  }
  res <- step.long(nt,ng=1,hx1,mt,ker)
  D0 <- res$J.r
  Ds0 <- res$J.n	# no serial correction
  cov = (fxxinv %*% D0 %*% fxxinv)/nt
  covs = (fxxinv %*% Ds0 %*% fxxinv)/nt # no serial correction
  return(list(cov = cov, covs = covs, D0 = D0, D1.inv = fxxinv, uhat = uhat))
}

band.true <- function(rho,rho.x,tau,nt,dy,ker){
  # calculate bandwidth values using true autocorrelation coefficients
  ctrue <- mvtnorm::pmvnorm(lower=-Inf,upper=c(0,0), mean=rep(-stats::qnorm(tau),2), sigma=matrix(c(1,rho,rho,1),nrow=2)) - tau^2
  del.h <- ctrue[1]/(tau*(1-tau))
  if(dy==1){
    delta <- del.h*rho.x	# true value of delta
    if(ker %in% c(1,2)){
      ah = ((2*delta)/(1-(delta^2)))^2
      band = round((1.1447*((ah*nt)^{1/3})),0)
    }
    if(ker %in% c(3)){
      ah = ((2*delta)/((1-delta)^{2}))^2
      band = round((1.3221*((ah*nt)^{1/5})),0)
    }
  }
  if(dy==2){
    delta = c(del.h,(del.h*rho.x))
    sig = c(0.502,0.601)
    if(ker %in% c(1,2)){
      e1 = 4*(delta^2)*(sig^4)/(((1-delta)^6)*((1+delta)^2))
      e2 = (sig^4)/((1-delta)^4)
      ah = sum(e1)/sum(e2)
      band = round((1.1447*((ah*nt)^{1/3})),0)
    }
    if(ker %in% 3){
      e1 = 4*(delta^2)*(sig^4)/((1-delta)^8)
      e2 = (sig^4)/((1-delta)^4)
      ah = sum(e1)/sum(e2)
      band = round((1.3221*((ah*nt)^{1/5})),0)
    }
  }
  band = pmax(band,1)
  return(list(delta = delta, a.hat=ah, band=band))
}

