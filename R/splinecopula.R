#' @title B-spline copula using the five M-spline basis functions
#' @description
#' \code{spline.copula} computes the B-spline copula (or its density) based on the five-parameter M-spline basis functions.
#' This copula is a specific instance of the B-spline copula family,
#' and can be implemented using matrix-based operations with \code{M.spline} and \code{I.spline} from the \pkg{joint.Cox} package.
#' @details
#' If \code{density = TRUE}, the function computes the copula **density** \(c(u, v)\); otherwise, it returns the **copula distribution function** \(C(u, v)\).
#' If \code{density = FALSE}, it returns the copula function. The implementation uses five M-spline or I-spline
#' basis functions defined on [0,1]. The coefficient matrix is fixed internally but can be modified in advanced use.
#'
#' @export
#' @importFrom stats nlm
#' @importFrom graphics lines par
#' @import joint.Cox
#' @param u A numeric vector of values in [0, 1] (uniform marginals for the first variable).
#' @param v A numeric vector of values in [0, 1] (uniform marginals for the second variable).
#' @param R A 5×5 non-negative coefficient matrix defining the copula structure. 
#'          The matrix must satisfy the following conditions:
#'          \itemize{
#'            \item All entries must be non-negative (\eqn{R_{kl} \ge 0}).
#'            \item The sum of all entries must be exactly 1.
#'            \item The row sums and column sums must equal:
#'                  \code{c(1/8, 1/4, 1/4, 1/4, 1/8)} (in order).
#'          }
#'          These conditions ensure that the resulting function is a valid copula density.
#'          You may also specify one of the built-in presets:
#'          "PE1", "PE2",–"PE3", "PN1", "PN2", "PN3", "I", NE1", "NE2", "NE3", "NN1", "NN2", "NN3".
#'          Default `"PE1"`.
#' @param mat Logical; if \code{TRUE}, returns the full matrix (outer product) of copula evaluations;
#' otherwise returns only the diagonal values, i.e., C(u_i, v_i) or c(u_i, v_i) for i = 1,...,n. Default is \code{FALSE}.
#' @param density Logical; if \code{TRUE}, evaluates the copula density; if \code{FALSE}, evaluates the copula function. Default is \code{FALSE}.
#' 
#' @return A numeric vector or matrix. If \code{mat=FALSE}, returns a vector of length \code{length(u)};
#' if \code{mat=TRUE}, returns a matrix of size \code{length(u)} x \code{length(v)}.
#' 
#' @seealso \code{\link[joint.Cox]{M.spline}}, \code{\link[joint.Cox]{I.spline}}
#' 
#' @examples
#' # Example data
#' library(joint.Cox)
#' library(ggplot2)
#' 
#' N = 50
#' u = v = seq(from = 0, to = 1, length.out = N)
#' U = rep(u, N)
#' V = rep(v, each = N)
#'
#' c.data = data.frame(U = U, V = V, C = spline.copula(U, V, R = "PE1", density=TRUE))
#' ggplot(aes(x = U, y = V), data = c.data) +
#'   geom_contour(aes(x = U, y = V, z = C, colour = after_stat(level)),
#'   data = c.data, bins=25) + xlab("u") + ylab("v")
#'

spline.copula = function(u, v, R = "PE1", mat = FALSE, density = FALSE){
  
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
  if(mat == TRUE){A}else{diag(A)}
}
