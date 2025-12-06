#' @title B-spline copula using the five M-spline basis functions
#' @description
#' \code{spline.copula} computes the B-spline copula (or its density) based on the
#' five-parameter M-spline basis functions. This copula is a specific instance of
#' the B-spline copula family and can be implemented using matrix-based operations
#' with \code{M.spline} and \code{I.spline} from the \pkg{joint.Cox} package.
#'
#' @details
#' If \code{density = TRUE}, the function computes the copula density \eqn{c(u, v)};
#' otherwise, it returns the copula distribution function \eqn{C(u, v)}.
#' If \code{density = FALSE}, it returns the copula function. The implementation
#' uses five M-spline or I-spline basis functions defined on \eqn{[0,1]}. The
#' coefficient matrix is fixed internally but can be modified in advanced use.
#'
#' @importFrom stats integrate nlm
#' @importFrom graphics lines par
#' @import joint.Cox
#' @export
#'
#' @param u A numeric vector of values in \eqn{[0, 1]} (uniform marginals for the
#'   first variable).
#' @param v A numeric vector of values in \eqn{[0, 1]} (uniform marginals for the
#'   second variable).
#' @param R A 5×5 non-negative coefficient matrix defining the copula structure.
#'   The matrix must satisfy the following conditions:
#'   \itemize{
#'     \item All entries must be non-negative (\eqn{R_{kl} \ge 0}).
#'     \item The sum of all entries must be exactly 1.
#'     \item The row sums and column sums must equal
#'           \code{c(1/8, 1/4, 1/4, 1/4, 1/8)} (in order).
#'   }
#'   These conditions ensure that the resulting function is a valid copula density.
#'   You may also specify one of the built-in presets:
#'   \code{"PE1"}, \code{"PE2"}, \code{"PE3"},
#'   \code{"PN1"}, \code{"PN2"}, \code{"PN3"},
#'   \code{"I"},
#'   \code{"NE1"}, \code{"NE2"}, \code{"NE3"},
#'   \code{"NN1"}, \code{"NN2"}, \code{"NN3"}.
#'   Default is \code{"PE1"}.
#'
#' @param mat Logical; if \code{TRUE}, returns the full matrix (outer product) of
#'   copula evaluations; otherwise returns only the diagonal values, i.e.,
#'   \eqn{C(u_i, v_i)} or \eqn{c(u_i, v_i)} for \eqn{i = 1, \ldots, n}.
#'   Default is \code{FALSE}.
#' @param density Logical; if \code{TRUE}, evaluates the copula density; if
#'   \code{FALSE}, evaluates the copula function. Default is \code{FALSE}.
#' @param Kendall Logical; if \code{TRUE}, returns Kendall’s tau in addition to
#'   copula values. Default is \code{FALSE}.
#' @param Spearman Logical; if \code{TRUE}, returns Spearman’s rho in addition to
#'   copula values. Default is \code{FALSE}.
#'
#' @return
#' If both \code{Kendall = FALSE} and \code{Spearman = FALSE} (default), returns:
#' \itemize{
#'   \item A numeric vector of length \code{length(u)} if \code{mat = FALSE}.
#'   \item A numeric matrix of size \code{length(u)} \eqn{\times} \code{length(v)}
#'         if \code{mat = TRUE}.
#' }
#'
#' If \code{Kendall = TRUE} or \code{Spearman = TRUE}, returns a list containing:
#' \itemize{
#'   \item \code{value}: A numeric vector or matrix representing the evaluated
#'         copula function or copula density.
#'   \item \code{Kendall_tau}: (Optional) Kendall’s tau, included only if
#'         \code{Kendall = TRUE}.
#'   \item \code{Spearman_rho}: (Optional) Spearman’s rho, included only if
#'         \code{Spearman = TRUE}.
#' }
#'
#' @seealso \code{\link[joint.Cox]{M.spline}}, \code{\link[joint.Cox]{I.spline}}
#'
#' @examples
#' N <- 50
#' u <- v <- seq(from = 0, to = 1, length.out = N)
#' U <- rep(u, N)
#' V <- rep(v, each = N)
#'
#' c.data <- data.frame(
#'   U = U,
#'   V = V,
#'   C = spline.copula(U, V, R = "PE1", density = TRUE)
#' )
#'
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   ggplot2::ggplot(c.data, ggplot2::aes(x = U, y = V)) +
#'     ggplot2::geom_contour(
#'       ggplot2::aes(z = C, colour = ggplot2::after_stat(level)),
#'       bins = 25
#'     ) +
#'     ggplot2::xlab("u") + ggplot2::ylab("v")
#' }



spline.copula = function(u, v, R = "PE1", mat = FALSE, density = FALSE, Kendall = FALSE, Spearman = FALSE){
  
  if(mat == FALSE && length(u) != length(v)){
    warning("the lengths of u and v differ.")
  }
  
  preset <- list(
    ## Positive Exchangeable
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
    
    ## Positive Non-exchangeable
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
    
    ## Independence
    I = matrix(c(1,2,2,2,1,
                 2,4,4,4,2,
                 2,4,4,4,2,
                 2,4,4,4,2,
                 1,2,2,2,1)/64, 5, 5, byrow = TRUE),
    
    ## Negative Exchangeable
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
    
    ## Negative Non-exchangeable 
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
  
  if (is.character(R)) {
    if (!R %in% names(preset))
      stop("Unknown preset name: ", R)
    R <- preset[[R]]
  }
  
  row_sums <- rowSums(R)
  col_sums <- colSums(R)
  
  expected_rows <- c(1/8, 1/4, 1/4, 1/4, 1/8)
  expected_cols <- c(1/8, 1/4, 1/4, 1/4, 1/8)
  
  if (any(abs(row_sums - expected_rows) > 1e-8)) {
    stop("Row sums of R must be: 1/8, 1/4, 1/4, 1/4, 1/8 (in order).")
  }
  if (any(abs(col_sums - expected_cols) > 1e-8)) {
    stop("Column sums of R must be: 1/8, 1/4, 1/4, 1/4, 1/8 (in order).")
  }
  if (!is.matrix(R) || !identical(dim(R), c(5L,5L))) {
    stop("R must be a 5x5 matrix")
  }
  if (any(R < 0)) {
    stop("ALL entries in R must be non-negative")
  }
  if (abs(sum(R) - 1) > 1e-8) {
    stop("Sum of all entries in R must be 1")
  }
  
  if (density == TRUE) {
    A = M.spline(u, xi1 = 0, xi3 = 1) %*% R %*% t(M.spline(v, xi1 = 0, xi3 = 1))
  }else{
    A = I.spline(u, xi1 = 0, xi3 = 1) %*% R %*% t(I.spline(v, xi1 = 0, xi3 = 1))
  }
  
  copula_values <- if (mat) A else diag(A)
  
  if (Kendall || Spearman) {
    result <- list()
    result$value <- copula_values
    
    if (Kendall) {
      E <- matrix(0, 5, 5)
      for (k in 1:5) {
        for (i in 1:5) {
          E[k, i] <- integrate(function(x) {
            I.spline(x, 0, 1)[, k] * M.spline(x, 0, 1)[, i]
          }, lower = 0, upper = 1)$value
        }
      }
      tau <- 0
      for (k in 1:5) {
        for (l in 1:5) {
          for (i in 1:5) {
            for (j in 1:5) {
              tau <- tau + R[k, l] * R[i, j] * E[k, i] * E[l, j]
            }
          }
        }
      }
      
      result$Kendall_tau <- 4 * tau - 1
    }
    
    if (Spearman) {
      F <- numeric(5)
      for (k in 1:5) {
        F[k] <- integrate(function(u) I.spline(u, 0, 1)[, k], lower = 0, upper = 1)$value
      }
      rho <- 0
      for (k in 1:5) {
        for (l in 1:5) {
          rho <- rho + R[k, l] * F[k] * F[l]
        }
      }
      result$Spearman_rho <- 12 * rho - 3
    }
    
    return(result)
  }
  return(copula_values)
}
