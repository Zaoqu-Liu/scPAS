#' Sparse Matrix Row Scaling (Internal Function)
#'
#' @description
#' Scale rows of a sparse matrix by centering and/or scaling.
#' Note: centering requires dense conversion as it destroys sparsity.
#'
#' @param x A sparse matrix (dgCMatrix or similar) or dense matrix
#' @param center Logical. Should rows be centered? (default: TRUE)
#' @param scale Logical. Should rows be scaled? (default: TRUE)
#'
#' @return A matrix with scaled rows (dense if centered, sparse if scale-only)
#'
#' @keywords internal
sparse_row_scale <- function(x, center = TRUE, scale = TRUE) {
  # Ensure we have a matrix
  if (!is.matrix(x) && !methods::is(x, "sparseMatrix")) {
    x <- as.matrix(x)
  }
  
  # Calculate row statistics using efficient sparse methods
  n_cols <- ncol(x)
  
  if (methods::is(x, "sparseMatrix")) {
    row_means <- Matrix::rowMeans(x)
    if (scale) {
      # Use sample variance formula: sum((x - mean)^2) / (n-1)
      # = (sum(x^2) - n * mean^2) / (n-1)
      # = n/(n-1) * (mean(x^2) - mean^2)
      row_sq_means <- Matrix::rowMeans(x^2)
      row_vars <- (row_sq_means - row_means^2) * n_cols / (n_cols - 1)
      # Handle numerical precision issues
      row_vars[row_vars < 0] <- 0
      row_sds <- sqrt(row_vars)
      # Avoid division by zero
      row_sds[row_sds == 0 | !is.finite(row_sds)] <- 1
    }
  } else {
    row_means <- rowMeans(x, na.rm = TRUE)
    if (scale) {
      row_sds <- apply(x, 1, stats::sd, na.rm = TRUE)
      row_sds[row_sds == 0 | is.na(row_sds) | !is.finite(row_sds)] <- 1
    }
  }
  
  if (!center && !scale) {
    return(x)
  }
  
  # Centering requires dense conversion (centering destroys sparsity)
  if (center) {
    if (methods::is(x, "sparseMatrix")) {
      # Convert to dense for centering
      x_scaled <- as.matrix(x)
    } else {
      x_scaled <- x
    }
    # Center rows (subtract row means)
    x_scaled <- x_scaled - row_means
  } else {
    x_scaled <- x
  }
  
  # Scale rows (divide by row sds)
  if (scale) {
    if (methods::is(x_scaled, "sparseMatrix")) {
      # Scale sparse matrix row-wise using diagonal multiplication
      D_inv <- Matrix::Diagonal(x = 1 / row_sds)
      x_scaled <- D_inv %*% x_scaled
    } else {
    x_scaled <- x_scaled / row_sds
    }
  }
  
  return(x_scaled)
}


#' Permutation Test (Internal Function)
#'
#' @description
#' Perform permutation test for statistical significance.
#' Supports parallel processing via the future package.
#'
#' @param scaled_exp Scaled expression matrix
#' @param Coefs Coefficient vector
#' @param permutation_times Number of permutations
#' @param n_cores Number of cores to use (default: 1 for sequential)
#' @param seed Random seed for reproducibility
#'
#' @return List containing risk scores and background distribution
#'
#' @keywords internal
parallel_permutation_test <- function(scaled_exp, Coefs, permutation_times, 
                                     n_cores = 1, seed = 12345) {
  
  # Calculate observed risk score
  risk_score <- crossprod(scaled_exp, Coefs)
  
  if (n_cores > 1) {
    # Use parallel processing
    if (!requireNamespace("future", quietly = TRUE)) {
      message("Package 'future' not available. Falling back to sequential processing.")
      n_cores <- 1
    } else if (!requireNamespace("future.apply", quietly = TRUE)) {
      message("Package 'future.apply' not available. Falling back to sequential processing.")
      n_cores <- 1
    }
  }
  
  if (n_cores > 1) {
    # Set up parallel backend
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    
    future::plan(future::multisession, workers = n_cores)
    
    # Parallel permutation
    set.seed(seed)
    randomPermutation <- future.apply::future_sapply(
      1:permutation_times,
      FUN = function(x) {
        set.seed(seed + x)
        sample(Coefs, length(Coefs), replace = FALSE)
      },
      future.seed = TRUE
    )
  } else {
    # Sequential permutation (original code)
    set.seed(seed)
    randomPermutation <- sapply(1:permutation_times, FUN = function(x){
      set.seed(seed + x)
      sample(Coefs, length(Coefs), replace = FALSE)
    })
  }
  
  # Ensure sparse format
  randomPermutation <- methods::as(randomPermutation, "sparseMatrix")
  
  # Calculate background distribution
  risk_score.background <- crossprod(scaled_exp, randomPermutation)
  
  return(list(
    risk_score = risk_score,
    risk_score.background = risk_score.background
  ))
}

