#' @title Random generation from the B-spline copula using five M-spline basis functions
#' @description
#' \code{spline.copula.simu} generates random samples \eqn{(U, V)} from the B-spline copula
#' defined by a 5×5 coefficient matrix \code{R}. The simulation uses the inverse transform
#' method based on the conditional distribution \eqn{V \mid U=u}.
#'
#' @details
#' Given \eqn{U \sim \mathrm{Uniform}(0,1)}, the conditional distribution function of
#' \eqn{V \mid U = u} is
#' \deqn{F_{V\mid U=u}(v) = M(u)^\top R I(v),}
#' where \eqn{M(\cdot)} and \eqn{I(\cdot)} are the five M-spline and I-spline basis vectors.
#' For each \eqn{u}, a draw \eqn{V} is obtained by solving
#' \eqn{F_{V\mid U=u}(v)=W} for \eqn{v}, where \eqn{W\sim U(0,1)}.
#'
#' If \code{report_tau = TRUE}, the function also returns:
#' \itemize{
#'   \item \code{tau_emp}: empirical Kendall's tau of the simulated values;
#'   \item \code{tau_theory}: theoretical Kendall's tau computed using
#'         \code{\link{spline.copula}} with \code{Kendall = TRUE}.
#' }
#'
#' If \code{report_rho = TRUE}, the function also returns:
#' \itemize{
#'   \item \code{rho_emp}: empirical Spearman's rho of the simulated values;
#'   \item \code{rho_theory}: theoretical Spearman's rho computed using
#'         \code{\link{spline.copula}} with \code{Spearman = TRUE}.
#' }
#'
#' @param n Integer. Number of samples to generate.
#' @param R A 5×5 non-negative coefficient matrix or a preset name:
#'          \code{"PE1"}, \code{"PE2"}, \code{"PE3"},
#'          \code{"PN1"}, \code{"PN2"}, \code{"PN3"},
#'          \code{"I"},
#'          \code{"NE1"}, \code{"NE2"}, \code{"NE3"},
#'          \code{"NN1"}, \code{"NN2"}, \code{"NN3"}.
#' @param seed Optional integer for reproducibility.
#' @param report_tau Logical. If TRUE, returns empirical and theoretical Kendall's tau.
#' @param report_rho Logical. If TRUE, returns empirical and theoretical Spearman's rho.
#' @param verbose Logical. If TRUE, prints correlation results.
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item \code{U}: simulated U-values;
#'   \item \code{V}: simulated V-values;
#'   \item \code{R}: user-specified R (preset name or matrix).
#' }
#' If \code{report_tau = TRUE}, the list also returns:
#' \itemize{
#'   \item \code{tau_emp}: empirical Kendall's tau;
#'   \item \code{tau_theory}: theoretical Kendall's tau.
#' }
#' If \code{report_rho = TRUE}, the list also returns:
#' \itemize{
#'   \item \code{rho_emp}: empirical Spearman's rho;
#'   \item \code{rho_theory}: theoretical Spearman's rho.
#' }
#'
#' @importFrom stats runif cor uniroot
#' @importFrom joint.Cox M.spline I.spline
#'
#' @examples
#' set.seed(123)
#' out <- spline.copula.simu(2000, R = "PE1",
#'                           report_tau = TRUE,
#'                           report_rho = TRUE)
#' out$tau_emp
#' out$tau_theory
#' out$rho_emp
#' out$rho_theory
#'
#' @export
spline.copula.simu <- function(n, R = "PE1", seed = NULL,
                               report_tau = TRUE,
                               report_rho = TRUE,
                               verbose = FALSE) {
  
  .presets <- local({
    list(
      PE1 = matrix(c(1,0,0,0,0,
                     0,2,0,0,0,
                     0,0,2,0,0,
                     0,0,0,2,0,
                     0,0,0,0,1)/8, 5, 5, byrow = TRUE),
      PE2 = matrix(c(1,2,2,0,0,
                     2,4,4,0,0,
                     2,4,4,0,0,
                     0,0,0,10,0,
                     0,0,0,0,5)/40, 5, 5, byrow = TRUE),
      PE3 = matrix(c(5,0,0,0,0,
                     0,10,0,0,0,
                     0,0,4,4,2,
                     0,0,4,4,2,
                     0,0,2,2,1)/40, 5, 5, byrow = TRUE),
      PN1 = matrix(c(1,0,0,0,0,
                     0,0,0,2,0,
                     0,2,0,0,0,
                     0,0,2,0,0,
                     0,0,0,0,1)/8, 5, 5, byrow = TRUE),
      PN2 = matrix(c(1,1,0,0,0,
                     1,0,2,1,0,
                     0,0,2,2,0,
                     0,3,0,0,1,
                     0,0,0,1,1)/16, 5, 5, byrow = TRUE),
      PN3 = matrix(c(4,2,2,0,0,
                     3,5,3,5,0,
                     1,6,6,0,3,
                     0,3,4,8,1,
                     0,0,1,3,4)/64, 5, 5, byrow = TRUE),
      I = matrix(c(1,2,2,2,1,
                   2,4,4,4,2,
                   2,4,4,4,2,
                   2,4,4,4,2,
                   1,2,2,2,1)/64, 5, 5, byrow = TRUE),
      NE1 = matrix(c(0,0,0,0,1,
                     0,0,0,2,0,
                     0,0,2,0,0,
                     0,2,0,0,0,
                     1,0,0,0,0)/8, 5, 5, byrow = TRUE),
      NE2 = matrix(c(0,0,2,2,1,
                     0,0,4,4,2,
                     0,0,4,4,2,
                     0,10,0,0,0,
                     5,0,0,0,0)/40, 5, 5, byrow = TRUE),
      NE3 = matrix(c(0,0,0,0,5,
                     0,0,0,10,0,
                     2,4,4,0,0,
                     2,4,4,0,0,
                     1,2,2,0,0)/40, 5, 5, byrow = TRUE),
      NN1 = matrix(c(0,0,0,0,1,
                     0,0,2,0,0,
                     0,0,0,2,0,
                     0,2,0,0,0,
                     1,0,0,0,0)/8, 5, 5, byrow = TRUE),
      NN2 = matrix(c(0,0,0,1,1,
                     0,3,0,0,1,
                     0,0,2,2,0,
                     1,0,2,1,0,
                     1,1,0,0,0)/16, 5, 5, byrow = TRUE),
      NN3 = matrix(c(0,0,1,3,4,
                     0,3,4,8,1,
                     1,6,6,0,3,
                     3,5,3,5,0,
                     4,2,2,0,0)/64, 5, 5, byrow = TRUE)
    )
  })
  
  .validate_R <- function(R){
    if (!is.matrix(R) || !identical(dim(R), c(5L,5L))) stop("R must be a 5x5 matrix")
    if (any(R < 0)) stop("Entries of R must be non-negative")
    if (abs(sum(R) - 1) > 1e-8) stop("Sum of all entries of R must be 1")
    expected <- c(1/8, 1/4, 1/4, 1/4, 1/8)
    if (any(abs(rowSums(R) - expected) > 1e-8)) stop("Incorrect row sums")
    if (any(abs(colSums(R) - expected) > 1e-8)) stop("Incorrect column sums")
    R
  }
  
  .get_R <- function(R) {
    if (is.character(R)) {
      if (!R %in% names(.presets)) stop("Unknown preset name: ", R)
      R <- .presets[[R]]
    }
    .validate_R(R)
  }
  
  R_input <- R
  
  if (!is.null(seed)) set.seed(seed)
  R <- .get_R(R_input)
  
  U <- stats::runif(n)
  V <- numeric(n)
  
  for (i in seq_len(n)) {
    u <- U[i]
    w <- stats::runif(1)
    Mu <- matrix(as.numeric(M.spline(u, xi1=0, xi3=1)), nrow = 1)
    
    f <- function(v){
      Iv <- matrix(as.numeric(I.spline(v, xi1=0, xi3=1)), ncol = 1)
      drop(Mu %*% R %*% Iv) - w
    }
    
    f0 <- -w
    Iv1 <- matrix(as.numeric(I.spline(1, xi1=0, xi3=1)), ncol = 1)
    f1 <- drop(Mu %*% R %*% Iv1) - w
    if (f0 > 0 || f1 < 0) w <- min(max(w, 1e-12), 1 - 1e-12)
    
    v_sol <- tryCatch(
      stats::uniroot(f, lower=0, upper=1, tol=1e-8, maxiter=100)$root,
      error = function(e) NA_real_
    )
    V[i] <- if (is.na(v_sol)) stats::runif(1) else v_sol
  }
  
  res <- list(U = U, V = V, R = R_input)
  
  if (report_tau || report_rho) {
    if (report_tau) {
      tau_emp <- suppressWarnings(stats::cor(U, V, method = "kendall"))
      res$tau_emp <- tau_emp
    }
    if (report_rho) {
      rho_emp <- suppressWarnings(stats::cor(U, V, method = "spearman"))
      res$rho_emp <- rho_emp
    }
    
    tmp <- spline.copula(
      u        = 0.5,
      v        = 0.5,
      R        = R_input,
      Kendall  = report_tau,
      Spearman = report_rho
    )
    
    if (report_tau) {
      tau_theory <- tmp$Kendall_tau
      res$tau_theory <- tau_theory
    }
    if (report_rho) {
      rho_theory <- tmp$Spearman_rho
      res$rho_theory <- rho_theory
    }
    
    if (verbose) {
      msg <- character()
      if (report_tau) {
        msg <- c(msg,
                 sprintf("Empirical Kendall's tau = %.4f; Theoretical tau = %.4f",
                         tau_emp, tau_theory))
      }
      if (report_rho) {
        msg <- c(msg,
                 sprintf("Empirical Spearman's rho = %.4f; Theoretical rho = %.4f",
                         rho_emp, rho_theory))
      }
      message(paste(msg, collapse = "\n"))
    }
  }
  
  res
}
