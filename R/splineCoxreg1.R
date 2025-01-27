#' @title Fitting the five-parameter spline Cox model giving a specified shape
#' @description
#' \code{splineCox.reg1} estimates the parameters of a five-parameter spline Cox model based on a specified shape for the baseline hazard function.
#' The function calculates the estimates for the model parameters (beta) and the baseline hazard scale parameter (gamma), using non-linear optimization.
#' If a numeric vector is provided for the \code{model} parameter, it will be normalized to have an L1 norm of 1.
#' Additionally, if \code{plot = TRUE}, the function generates a plot of the estimated baseline hazard function with 95% confidence intervals.
#' The x-axis represents time, and the y-axis represents the estimated hazard.
#' The solid line indicates the estimated hazard function, while the dashed red lines represent the confidence intervals.
#' @export
#' @importFrom stats nlm
#' @importFrom graphics lines par
#' @import joint.Cox
#' @param t.event a vector for time-to-event
#' @param event a vector for event indicator (=1 event; =0 censoring)
#' @param Z a matrix for covariates; nrow(Z)=sample size, ncol(Z)=the number of covariates
#' @param xi1 lower bound for the hazard function; the default is \code{min(t.event)}
#' @param xi3 upper bound for the hazard function; the default is \code{max(t.event)}
#' @param model A character string specifying the shape of the baseline hazard function or a numeric vector of length 5 representing custom weights.
#'              If a numeric vector is provided, it will be normalized to have an L1 norm of 1.
#'              Available options include:
#'              "increase", "constant", "decrease", "unimodal1", "unimodal2", "unimodal3", "bathtub1", "bathtub2", "bathtub3".
#'              Default is "constant"
#' @param plot A logical value indicating whether to plot the estimated baseline hazard function.
#'             If \code{TRUE}, a plot is generated displaying the estimated baseline hazard function along with its 95% confidence intervals.
#'             The x-axis represents time, and the y-axis represents the estimated hazard.
#'             The solid line indicates the estimated hazard function, while the dashed red lines represent the confidence intervals.
#'             Default is \code{TRUE}.
#' @param p0 Initial values to maximize the likelihood (1 + p parameters; baseline hazard scale parameter and p regression coefficients)
#' @return A list containing the following components:
#'   \item{model}{A shape of the baseline hazard function or the normalized custom numeric vector used.}
#'   \item{parameter}{A numeric vector of the parameters defining the baseline hazard shape.}
#'   \item{beta}{A named vector with the estimates, standard errors, and 95% confidence intervals for the regression coefficients}
#'   \item{gamma}{A named vector with the estimate, standard error, and 95% confidence interval for the baseline hazard parameter}
#'   \item{loglik}{A named vector containing the log-likelihood (\code{LogLikelihood}),
#'                 Akaike Information Criterion (\code{AIC}), and Bayesian Information Criterion (\code{BIC})}
#'   \item{plot}{A baseline hazard function plot (if \code{plot = TRUE}).}
#' @examples
#' # Example data
#' library(joint.Cox)
#' data(dataOvarian)
#' t.event = dataOvarian$t.event
#' event = dataOvarian$event
#' Z = dataOvarian$CXCL12
#'
#' reg1 <- splineCox.reg1(t.event, event, Z, model = "constant")
#' print(reg1)
#'
splineCox.reg1 <- function (t.event, event, Z, xi1 = min(t.event), xi3 = max(t.event),
                            model = "constant", p0 = rep(0, 1 + ncol(as.matrix(Z))),
                            plot = TRUE)
{
  shape.list = list(  increase  = c(0.05,0.1,0.15,0.3,0.4),
                      constant  = c(0.125,0.25,0.25,0.25,0.125),
                      decrease  = c(0.4,0.3,0.15,0.1,0.05),
                      unimodal1 = c(0.001,0.014,0.97,0.014,0.001),
                      unimodal2 = c(0.001,0.8,0.124,0.074,0.001),
                      unimodal3 = c(0.001,0.074,0.124,0.8,0.001),
                      unimodal4 = c(0.099, 0.18, 0.7, 0.02, 0.001),
                      unimodal5 = c(0.001, 0.02, 0.7, 0.18, 0.099),
                      bathtub1 = c(0.3, 0.1995, 0.001, 0.1995, 0.3),
                      bathtub2 = c(0.3, 0.001, 0.1, 0.299, 0.3),
                      bathtub3 = c(0.3, 0.299, 0.1, 0.001, 0.3),
                      bathtub4 = c(0.2, 0.001, 0.1, 0.299, 0.4),
                      bathtub5 = c(0.4, 0.299, 0.1, 0.001, 0.2)
  )
  norm_L1 = function(vec) {
    if(any(vec < 0)) {
      stop("All elements of the input vector must be non-negative.")
    }
    L1 = sum(abs(vec))
    if (L1 != 1) {
      vec = vec / L1
    }
    return(vec)
  }

  if (is.character(model) && model %in% names(shape.list)) {
    para = shape.list[[model]]
  } else if (is.numeric(model) && length(model) == 5) {
    para = norm_L1(model)
  } else {
    stop("Invalid model specification. Provide a valid model name or a numeric vector of length 5.")
  }

  d = event
  Z = as.matrix(Z)
  p = ncol(Z)
  l.func = function(phi) {
    b = phi[2:(1 + p)]
    g = exp(pmin(phi[1], 500))
    r = as.vector(M.spline(t.event, xi1 = min(t.event), xi3 = max(t.event)) %*%
                    (g * para))
    R = as.vector(I.spline(t.event, xi1 = min(t.event), xi3 = max(t.event)) %*%
                    (g * para))
    bZ = as.numeric(Z %*% b)
    l = sum(d * (log(r) + bZ))
    l = l - sum(pmin(exp(bZ) * R, exp(500)))
    -l
  }
  res = nlm(l.func, p = p0, hessian = TRUE)
  LogLik = -res$minimum
  k = p + 1
  n = length(d)
  AIC = -2 * LogLik + 2 * k
  BIC = -2 * LogLik + log(n) * k
  loglik.res = c(LogLikelihood = LogLik, AIC = AIC, BIC = BIC)
  beta.est = res$est[2:(1 + p)]
  gam.est = exp(res$est[1])
  H = -res$hessian
  V = solve(-H, tol = 10^(-100))
  beta.se = sqrt(diag(V)[2:(1 + p)])
  gam.se = sqrt(gam.est %*% V[1, 1] %*% gam.est)
  b.lower = beta.est - 1.96 * beta.se
  b.upper = beta.est + 1.96 * beta.se
  g.lower = gam.est - 1.96 * gam.se
  g.upper = gam.est + 1.96 * gam.se
  beta.res = c(estimate = beta.est, SE = beta.se, Lower = b.lower,
               Upper = b.upper)
  gam.res = c(estimate = gam.est, SE = gam.se, Lower = g.lower,
              Upper = g.upper)
  if (plot) {
    t.seq = seq(xi1, xi3, length.out = 1000)
    t.splineplot=function(time){
      as.numeric(M.spline(time,xi1,xi3)%*%(gam.est * para))
    }
    splLow=function(time){
      haz.est = as.numeric(M.spline(time,xi1,xi3) %*% (gam.est * para))
      haz.se = as.numeric(gam.se * sum(para * M.spline(time, xi1, xi3)))
      return(haz.est - 1.96 * haz.se)
    }
    splUp=function(time){
      haz.est = as.numeric(M.spline(time,xi1,xi3) %*% (gam.est * para))
      haz.se = as.numeric(gam.se * sum(para * M.spline(time, xi1, xi3)))
      return(haz.est + 1.96 * haz.se)
    }
    ylower = min(sapply(t.seq,splLow))
    yupper = max(sapply(t.seq,splUp))
    margin = (yupper - ylower) * 0.1
    ylim.range = c(ylower - margin, yupper + margin)
    plot(t.seq, sapply(t.seq,t.splineplot), type = "l", lwd = 3, ylim = ylim.range,
         xlab = "time", ylab = "baseline hazard")
    lines(t.seq,sapply(t.seq,splLow),col="red",lty = "dashed")
    lines(t.seq,sapply(t.seq,splUp),col="red", lty = "dashed")
  }
  list(model = model, parameter = para, beta = beta.res, gamma = gam.res,
       loglik = loglik.res)
}


