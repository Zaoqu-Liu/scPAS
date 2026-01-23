# ============================================================================
# scPAS C++ 数值稳定性测试
# ============================================================================

cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("scPAS v1.0.3 C++ 数值稳定性测试\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

# 加载包
library(scPAS)

# 访问内部函数
APML0 <- scPAS:::APML0

tests_passed <- 0
tests_failed <- 0

test_that <- function(description, expr) {
  cat(sprintf("测试: %s ... ", description))
  tryCatch({
    result <- eval(expr)
    if (isTRUE(result)) {
      cat("✓ 通过\n")
      tests_passed <<- tests_passed + 1
    } else {
      cat("✗ 失败\n")
      tests_failed <<- tests_failed + 1
    }
  }, error = function(e) {
    cat(sprintf("✗ 错误: %s\n", e$message))
    tests_failed <<- tests_failed + 1
  })
}

# ============================================================================
# 测试 C++ 函数的数值稳定性
# ============================================================================
cat("\n[1] C++ APML0 数值稳定性测试\n")
cat(paste(rep("-", 50), collapse=""), "\n")

# 测试1: 常规数据
test_that("常规数据 - Gaussian回归", {
  set.seed(123)
  n <- 50
  p <- 20
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1:3] %*% c(1, 2, 3) + rnorm(n)
  
  result <- tryCatch({
    fit <- APML0(X, y, family = "gaussian", penalty = "Enet")
    all(is.finite(fit$Beta))
  }, error = function(e) FALSE)
  
  result
})

# 测试2: 具有零方差列的数据
test_that("零方差列处理 - Gaussian回归", {
  set.seed(123)
  n <- 50
  p <- 20
  X <- matrix(rnorm(n * p), n, p)
  X[, 5] <- 5  # 常量列
  y <- X[, 1:3] %*% c(1, 2, 3) + rnorm(n)
  
  result <- tryCatch({
    fit <- APML0(X, y, family = "gaussian", penalty = "Enet")
    all(is.finite(fit$Beta))
  }, error = function(e) FALSE)
  
  result
})

# 测试3: 稀疏数据
test_that("稀疏数据处理 - Gaussian回归", {
  set.seed(123)
  n <- 100
  p <- 50
  X <- matrix(rnorm(n * p), n, p)
  X[X < 1] <- 0  # 稀疏化
  y <- rnorm(n)
  
  result <- tryCatch({
    fit <- APML0(X, y, family = "gaussian", penalty = "Enet")
    all(is.finite(fit$Beta))
  }, error = function(e) FALSE)
  
  result
})

# 测试4: 二分类数据
test_that("二分类数据 - Logistic回归", {
  set.seed(123)
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), n, p)
  prob <- plogis(X[, 1] + X[, 2])
  y <- rbinom(n, 1, prob)
  
  result <- tryCatch({
    fit <- APML0(X, y, family = "binomial", penalty = "Enet")
    all(is.finite(fit$Beta))
  }, error = function(e) FALSE)
  
  result
})

# 测试5: 极端值处理
test_that("极端值处理 - Gaussian回归", {
  set.seed(123)
  n <- 50
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  X[1, 1] <- 1e6  # 添加极端值
  y <- rnorm(n)
  
  result <- tryCatch({
    fit <- APML0(X, y, family = "gaussian", penalty = "Enet")
    all(is.finite(fit$Beta))
  }, error = function(e) FALSE)
  
  result
})

# 测试6: 网络惩罚
cat("\n[2] 网络惩罚测试\n")
cat(paste(rep("-", 50), collapse=""), "\n")

test_that("网络惩罚 - Gaussian回归", {
  set.seed(123)
  n <- 50
  p <- 20
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1:3] %*% c(1, 2, 3) + rnorm(n)
  # 创建邻接矩阵
  Omega <- diag(p)
  Omega[1, 2] <- Omega[2, 1] <- 0.5
  Omega[2, 3] <- Omega[3, 2] <- 0.5
  
  result <- tryCatch({
    fit <- APML0(X, y, family = "gaussian", penalty = "Net", Omega = Omega)
    all(is.finite(fit$Beta))
  }, error = function(e) FALSE)
  
  result
})

test_that("网络惩罚 - Logistic回归", {
  set.seed(123)
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), n, p)
  y <- rbinom(n, 1, 0.5)
  Omega <- diag(p)
  
  result <- tryCatch({
    fit <- APML0(X, y, family = "binomial", penalty = "Net", Omega = Omega)
    all(is.finite(fit$Beta))
  }, error = function(e) FALSE)
  
  result
})

# 测试7: Cox回归
cat("\n[3] Cox回归测试\n")
cat(paste(rep("-", 50), collapse=""), "\n")

test_that("Cox回归 - 基本功能", {
  set.seed(123)
  n <- 100
  p <- 20
  X <- matrix(rnorm(n * p), n, p)
  time <- rexp(n, 0.1)
  status <- sample(0:1, n, replace = TRUE)
  y <- cbind(time, status)
  
  result <- tryCatch({
    fit <- APML0(X, y, family = "cox", penalty = "Enet")
    all(is.finite(fit$Beta))
  }, error = function(e) FALSE)
  
  result
})

test_that("Cox回归 - 网络惩罚", {
  set.seed(123)
  n <- 100
  p <- 15
  X <- matrix(rnorm(n * p), n, p)
  time <- rexp(n, 0.1)
  status <- sample(0:1, n, replace = TRUE)
  y <- cbind(time, status)
  Omega <- diag(p)
  
  result <- tryCatch({
    fit <- APML0(X, y, family = "cox", penalty = "Net", Omega = Omega)
    all(is.finite(fit$Beta))
  }, error = function(e) FALSE)
  
  result
})

# 测试8: 多个alpha值
cat("\n[4] 多参数测试\n")
cat(paste(rep("-", 50), collapse=""), "\n")

test_that("多个alpha值", {
  set.seed(123)
  n <- 50
  p <- 20
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  
  all_pass <- TRUE
  for (alpha in c(0, 0.25, 0.5, 0.75, 1)) {
    result <- tryCatch({
      fit <- APML0(X, y, family = "gaussian", penalty = "Enet", alpha = alpha)
      all(is.finite(fit$Beta))
    }, error = function(e) FALSE)
    if (!result) all_pass <- FALSE
  }
  
  all_pass
})

test_that("多个nlambda值", {
  set.seed(123)
  n <- 50
  p <- 20
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  
  all_pass <- TRUE
  for (nlam in c(10, 50, 100)) {
    result <- tryCatch({
      fit <- APML0(X, y, family = "gaussian", penalty = "Enet", nlambda = nlam)
      # lambda存储在fit$fit$lambda中
      length(fit$fit$lambda) == nlam && all(is.finite(fit$Beta))
    }, error = function(e) FALSE)
    if (!result) all_pass <- FALSE
  }
  
  all_pass
})

# ============================================================================
# 测试总结
# ============================================================================
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("C++ 数值稳定性测试总结\n")
cat(paste(rep("=", 70), collapse=""), "\n")

total_tests <- tests_passed + tests_failed
cat(sprintf("\n总测试数: %d\n", total_tests))
cat(sprintf("通过: %d (%.1f%%)\n", tests_passed, 100 * tests_passed / total_tests))
cat(sprintf("失败: %d (%.1f%%)\n", tests_failed, 100 * tests_failed / total_tests))

if (tests_failed == 0) {
  cat("\n✅ 所有C++数值稳定性测试通过！\n")
} else {
  cat("\n⚠️ 有测试失败\n")
}

cat("\n")
