#' @title Fitting the five-parameter spline Cox model with a specified shape, selecting the best fit
#' @description
#' \code{splineCox.reg2} estimates the parameters of a five-parameter spline Cox model for multiple specified shapes
#' and selects the best-fitting model based on the maximization of the log-likelihood function.
#' This function supports predefined model shapes and custom numeric vectors of length 5.
#' If numeric vectors are provided, they will be normalized to have an L1 norm of 1.
#' Additionally, if \code{plot = TRUE}, the function generates a plot of the estimated baseline hazard function for the best-fitting model,
#' along with its 95% confidence intervals.
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
#' @param model A list of character strings and/or numeric vectors of length 5 specifying the shapes of the baseline hazard function to evaluate.
#'              Character options include:
#'              "increase", "constant", "decrease", "unimodal1", "unimodal2", "unimodal3", "bathtub1", "bathtub2", "bathtub3".
#'              Numeric vectors must be of length 5 and will be normalized to have an L1 norm of 1.
#'              Default is \code{names(shape.list)}, which includes all predefined models.
#' @param plot A logical value indicating whether to plot the estimated baseline hazard function.
#'             If \code{TRUE}, a plot is generated displaying the estimated baseline hazard function along with its 95% confidence intervals.
#'             The x-axis represents time, and the y-axis represents the estimated hazard.
#'             The solid line indicates the estimated hazard function, while the dashed red lines represent the confidence intervals.
#'             Default is \code{TRUE}.
#' @param p0 Initial values to maximize the likelihood (1 + p parameters; baseline hazard scale parameter and p regression coefficients)
#' @return A list containing the following components:
#'   \item{model}{A character string indicating the shape of the baseline hazard function used.}
#'   \item{parameter}{A numeric vector of the parameters defining the baseline hazard shape.}
#'   \item{beta}{A named vector with the estimates, standard errors, and 95% confidence intervals for the regression coefficients}
#'   \item{gamma}{A named vector with the estimate, standard error, and 95% confidence interval for the baseline hazard parameter}
#'   \item{loglik}{A named vector containing the log-likelihood (\code{LogLikelihood}),
#'                 Akaike Information Criterion (\code{AIC}), and Bayesian Information
#'                 Criterion (\code{BIC}) for the best-fitting model}
#'   \item{other_models}{A data frame containing the log-likelihood (\code{LogLikelihood}) for all other evaluated models,
#'                            with model names as row names.}
#'   \item{plot}{A baseline hazard function plot for the best-fitting model (if \code{plot = TRUE}).}
#' @examples
#' # Example data
#' library(joint.Cox)
#' data(dataOvarian)
#' t.event = dataOvarian$t.event
#' event = dataOvarian$event
#' Z = dataOvarian$CXCL12
#'
#' M = c("constant", "increase", "decrease")
#' reg2 <- splineCox.reg2(t.event, event, Z, model = M)
#' print(reg2)
#'
splineCox.reg2 <- function (t.event, event, Z, xi1 = min(t.event), xi3 = max(t.event),
                            model = names(shape.list), p0 = rep(0, 1 + ncol(as.matrix(Z))),
                            plot = TRUE)
{
  shape.list = list(
    increase  = c(0.05, 0.1, 0.15, 0.3, 0.4),
    constant  = c(0.125, 0.25, 0.25, 0.25, 0.125),
    decrease  = c(0.4, 0.3, 0.15, 0.1, 0.05),
    unimodal1 = c(0.001, 0.014, 0.97, 0.014, 0.001),
    unimodal2 = c(0.001, 0.8, 0.124, 0.074, 0.001),
    unimodal3 = c(0.001, 0.074, 0.124, 0.8, 0.001),
    unimodal4 = c(0.099, 0.18, 0.7, 0.02, 0.001),
    unimodal5 = c(0.001, 0.02, 0.7, 0.18, 0.099),
    bathtub1 = c(0.3, 0.1995, 0.001, 0.1995, 0.3),
    bathtub2 = c(0.3, 0.001, 0.1, 0.299, 0.3),
    bathtub3 = c(0.3, 0.299, 0.1, 0.001, 0.3),
    bathtub4 = c(0.2, 0.001, 0.1, 0.299, 0.4),
    bathtub5 = c(0.4, 0.299, 0.1, 0.001, 0.2)
  )

  results <- lapply(model, function(model) {
    if (is.character(model) && model %in% names(shape.list)) {
      splineCox.reg1(t.event, event, Z, xi1, xi3, model, p0, plot = FALSE)
    } else if (is.numeric(model) && length(model) == 5) {
      splineCox.reg1(t.event, event, Z, xi1, xi3, model, p0, plot = FALSE)
    } else {
      stop("Invalid model specification. Provide a valid model name or a numeric vector of length 5.")
    }
  })

  loglik.values = sapply(results, function(res) res$loglik["LogLikelihood"])
  best.index = which.max(loglik.values)
  best.result = results[[best.index]]
  other.loglik = loglik.values[-best.index]
  other.loglik = data.frame(Loglikelihodd = other.loglik)
  rownames(other.loglik) = model[-best.index]
  if (plot) {
    t.seq = seq(xi1, xi3, length.out = 1000)
    t.splineplot=function(time){
      as.numeric(M.spline(time,xi1,xi3)%*%(best.result$gamma["estimate"] * best.result$parameter))
    }
    splLow=function(time){
      haz.est = as.numeric(M.spline(time,xi1,xi3) %*% (best.result$gamma["estimate"] * best.result$parameter))
      haz.se = as.numeric(best.result$gamma["estimate"] * sum(best.result$parameter * M.spline(time, xi1, xi3)))
      return(haz.est - 1.96 * haz.se)
    }
    splUp=function(time){
      haz.est = as.numeric(M.spline(time,xi1,xi3) %*% (best.result$gamma["estimate"] * best.result$parameter))
      haz.se = as.numeric(best.result$gamma["estimate"] * sum(best.result$parameter * M.spline(time, xi1, xi3)))
      return(haz.est + 1.96 * haz.se)
    }
    ylower = min(sapply(t.seq,splLow))
    yupper = max(sapply(t.seq,splUp))
    margin = (yupper - ylower) * 0.1
    ylim.range = c(ylower - margin, yupper + margin)
    par(mfrow=c(1,1))
    plot(t.seq, sapply(t.seq,t.splineplot), type = "l", lwd = 3, ylim = ylim.range,
         xlab = "time", ylab = "baseline hazard")
    lines(t.seq,sapply(t.seq,splLow),col="red",lty = "dashed")
    lines(t.seq,sapply(t.seq,splUp),col="red", lty = "dashed")
  }
  c(best.result, list(other_loglik = other.loglik))
}
