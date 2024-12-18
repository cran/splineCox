#' @title Fitting the five-parameter spline Cox model with a specified shape, selecting the best fit
#' @description
#' \code{splineCox.reg2} estimates the parameters of a five-parameter spline Cox model for multiple specified shapes
#' and selects the best fitting model based on the minimization of the log-likelihood function.
#' The function calculates the estimates for the model parameters (beta) and the baseline hazard scale parameter (gamma), using non-linear optimization.
#' @export
#' @importFrom stats nlm
#' @import joint.Cox
#' @param t.event a vector for time-to-event
#' @param event a vector for event indicator (=1 event; =0 censoring)
#' @param Z a matrix for covariates; nrow(Z)=sample size, ncol(Z)=the number of covariates
#' @param xi1 lower bound for the hazard function; the default is min(t.event)
#' @param xi3 upper bound for the hazard function; the default is max(t.event)
#' @param model A character vector specifying which model shapes to consider for the baseline hazard.
#'              Available options are:
#'              "increase", "constant", "decrease", "unimodal1", "unimodal2", "unimodal3", "bathtub1", "bathtub2", "bathtub3".
#'              Default is \code{names(shape.list)} which includes all available models.
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
splineCox.reg2 <- function(t.event, event, Z, xi1 = min(t.event), xi3 = max(t.event),
                           model = names(shape.list), p0 = rep(0, 1 + ncol(as.matrix(Z))))
{
  shape.list = list(
    increase  = c(0.05, 0.1, 0.15, 0.3, 0.4),
    constant  = c(0.125, 0.25, 0.25, 0.25, 0.125),
    decrease  = c(0.4, 0.3, 0.15, 0.1, 0.05),
    unimodal1 = c(0.001, 0.014, 0.97, 0.014, 0.001),
    unimodal2 = c(0.001, 0.8, 0.124, 0.074, 0.001),
    unimodal3 = c(0.001, 0.074, 0.124, 0.8, 0.001),
    bathtub1  = c(0.3, 0.1995, 0.001, 0.1995, 0.3),
    bathtub2  = c(0.3, 0.001, 0.1009, 0.299, 0.3),
    bathtub3  = c(0.3, 0.299, 0.1009, 0.001, 0.3)
  )

  if (any(!model %in% names(shape.list))) {
    stop("Invalid model names. Choose from: ", paste(names(shape.list), collapse = ", "))
  }

  results = lapply(model, function(model) {
    splineCox.reg1(t.event, event, Z, xi1, xi3, model, p0)
  })

  loglik.values = sapply(results, function(res) res$loglik["LogLikelihood"])
  best.index = which.max(loglik.values)
  best.result = results[[best.index]]

  other.loglik = loglik.values[-best.index]
  other.loglik = data.frame(Loglikelihodd = other.loglik)
  rownames(other.loglik) = model[-best.index]
  c(best.result, list(other_loglik = other.loglik))
}
