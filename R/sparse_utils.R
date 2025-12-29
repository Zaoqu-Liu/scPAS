#' Sparse Matrix Row Scaling (Internal Function)
#'
#' @description
#' Efficiently scale rows of a sparse matrix while preserving sparsity.
#' This function is optimized for large sparse matrices and avoids converting
#' to dense format.
#'
#' @param x A sparse matrix (dgCMatrix or similar)
#' @param center Logical. Should rows be centered? (default: TRUE)
#' @param scale Logical. Should rows be scaled? (default: TRUE)
#'
#' @return A sparse matrix with scaled rows
#'
#' @keywords internal
sparse_row_scale <- function(x, center = TRUE, scale = TRUE) {
  if (!methods::is(x, "sparseMatrix")) {
    x <- methods::as(x, "sparseMatrix")
  }
  
  # Calculate row statistics
  if (center) {
    row_means <- Matrix::rowMeans(x)
  } else {
    row_means <- rep(0, nrow(x))
  }
  
  if (scale) {
    # Calculate row standard deviations efficiently
    row_vars <- Matrix::rowMeans(x^2) - row_means^2
    row_sds <- sqrt(pmax(row_vars, 0))  # Avoid negative values from floating point errors
    row_sds[row_sds == 0] <- 1  # Avoid division by zero
  } else {
    row_sds <- rep(1, nrow(x))
  }
  
  # Scale: (x - mean) / sd
  # Use sparse matrix operations to preserve sparsity
  x_scaled <- x
  
  # Center
  if (center) {
    x_scaled <- x_scaled - row_means
  }
  
  # Scale
  if (scale) {
    x_scaled <- x_scaled / row_sds
  }
  
  # Ensure it's still sparse
  x_scaled <- methods::as(x_scaled, "sparseMatrix")
  
  return(x_scaled)
}


#' Parallel Permutation Test (Internal Function)
#'
#' @description
#' Perform permutation test in parallel to speed up computation.
#' This function uses the future package for parallel processing.
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

