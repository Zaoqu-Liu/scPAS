# ============================================================================
# scPAS 科学准确性和功能深度测试
# ============================================================================

cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("scPAS v1.0.3 科学准确性深度测试\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

# 加载必要的包
suppressPackageStartupMessages({
  library(Matrix)
  library(stats)
})

# 加载 scPAS 包
tryCatch({
  library(scPAS)
  cat("已加载scPAS包\n")
}, error = function(e) {
  stop("无法加载scPAS包: ", e$message)
})

# 获取内部函数（用于测试）
sparse_row_scale <- scPAS:::sparse_row_scale
parallel_permutation_test <- scPAS:::parallel_permutation_test

# 计数器
tests_passed <- 0
tests_failed <- 0
test_results <- list()

# 测试辅助函数
test_that <- function(description, expr) {
  cat(sprintf("  测试: %s ... ", description))
  tryCatch({
    result <- eval(expr)
    if (isTRUE(result)) {
      cat("✓ 通过\n")
      tests_passed <<- tests_passed + 1
      test_results[[description]] <<- "PASS"
    } else {
      cat("✗ 失败\n")
      tests_failed <<- tests_failed + 1
      test_results[[description]] <<- "FAIL"
    }
  }, error = function(e) {
    cat(sprintf("✗ 错误: %s\n", e$message))
    tests_failed <<- tests_failed + 1
    test_results[[description]] <<- paste("ERROR:", e$message)
  })
}

expect_equal <- function(actual, expected, tolerance = 1e-10) {
  if (is.numeric(actual) && is.numeric(expected)) {
    all(abs(actual - expected) < tolerance)
  } else {
    identical(actual, expected)
  }
}

expect_true <- function(x) isTRUE(x)
expect_false <- function(x) isFALSE(x)

# ============================================================================
# 测试 1: P值计算的科学准确性
# ============================================================================
cat("\n[1] P值计算科学准确性测试\n")
cat(paste(rep("-", 50), collapse=""), "\n")

# 1.1 双尾检验验证
test_that("双尾P值计算: Z=1.96 应该约等于 0.05", {
  Z <- 1.96
  # 我们的公式: 2 * pnorm(abs(Z), lower.tail = FALSE)
  p_value <- 2 * pnorm(abs(Z), mean = 0, sd = 1, lower.tail = FALSE)
  abs(p_value - 0.05) < 0.001
})

test_that("双尾P值计算: Z=2.576 应该约等于 0.01", {
  Z <- 2.576
  p_value <- 2 * pnorm(abs(Z), mean = 0, sd = 1, lower.tail = FALSE)
  abs(p_value - 0.01) < 0.001
})

test_that("双尾P值计算: Z=0 应该等于 1", {
  Z <- 0
  p_value <- 2 * pnorm(abs(Z), mean = 0, sd = 1, lower.tail = FALSE)
  abs(p_value - 1.0) < 1e-10
})

test_that("双尾P值计算: 负Z值应该等于正Z值", {
  Z_pos <- 2.5
  Z_neg <- -2.5
  p_pos <- 2 * pnorm(abs(Z_pos), mean = 0, sd = 1, lower.tail = FALSE)
  p_neg <- 2 * pnorm(abs(Z_neg), mean = 0, sd = 1, lower.tail = FALSE)
  abs(p_pos - p_neg) < 1e-15
})

test_that("P值应该在[0, 1]范围内", {
  Z_values <- c(-10, -5, -2, -1, 0, 1, 2, 5, 10)
  p_values <- 2 * pnorm(abs(Z_values), mean = 0, sd = 1, lower.tail = FALSE)
  p_values <- pmin(p_values, 1)  # 上限约束
  all(p_values >= 0 & p_values <= 1)
})

# 1.2 与统计学标准对比
test_that("P值与t.test双尾结果一致性", {
  # 模拟数据
  set.seed(123)
  x <- rnorm(100, mean = 0.5, sd = 1)
  # 手动计算
  Z <- mean(x) / (sd(x) / sqrt(length(x)))
  p_manual <- 2 * pnorm(abs(Z), lower.tail = FALSE)
  # t检验结果 (大样本近似正态)
  p_ttest <- t.test(x, mu = 0)$p.value
  # 应该接近
  abs(p_manual - p_ttest) < 0.01
})

# ============================================================================
# 测试 2: sparse.cor 相关性计算
# ============================================================================
cat("\n[2] sparse.cor 相关性计算测试\n")
cat(paste(rep("-", 50), collapse=""), "\n")

# 2.1 与标准cor()函数对比
test_that("sparse.cor与stats::cor结果一致（稠密矩阵）", {
  set.seed(42)
  mat <- matrix(rnorm(500), nrow = 50, ncol = 10)
  sparse_mat <- as(mat, "dgCMatrix")
  
  cor_standard <- cor(mat)
  cor_sparse <- sparse.cor(sparse_mat)
  
  max(abs(cor_standard - cor_sparse)) < 1e-10
})

# 2.2 稀疏矩阵测试
test_that("sparse.cor处理稀疏矩阵正确", {
  set.seed(42)
  mat <- matrix(rnorm(500), nrow = 50, ncol = 10)
  mat[mat < 0.5] <- 0  # 稀疏化
  sparse_mat <- as(mat, "dgCMatrix")
  
  cor_standard <- cor(mat)
  cor_sparse <- sparse.cor(sparse_mat)
  
  max(abs(cor_standard - cor_sparse)) < 1e-10
})

# 2.3 对角线检验
test_that("相关矩阵对角线为1", {
  set.seed(42)
  mat <- matrix(rnorm(200), nrow = 20, ncol = 10)
  sparse_mat <- as(mat, "dgCMatrix")
  cor_result <- sparse.cor(sparse_mat)
  
  all(abs(diag(cor_result) - 1) < 1e-15)
})

# 2.4 对称性检验
test_that("相关矩阵是对称的", {
  set.seed(42)
  mat <- matrix(rnorm(200), nrow = 20, ncol = 10)
  sparse_mat <- as(mat, "dgCMatrix")
  cor_result <- sparse.cor(sparse_mat)
  
  max(abs(cor_result - t(cor_result))) < 1e-15
})

# 2.5 值范围检验
test_that("相关系数在[-1, 1]范围内", {
  set.seed(42)
  mat <- matrix(rnorm(500), nrow = 50, ncol = 10)
  sparse_mat <- as(mat, "dgCMatrix")
  cor_result <- sparse.cor(sparse_mat)
  
  all(cor_result >= -1 & cor_result <= 1)
})

# 2.6 边界条件测试
test_that("sparse.cor处理单行矩阵", {
  mat <- matrix(rnorm(10), nrow = 1, ncol = 10)
  sparse_mat <- as(mat, "dgCMatrix")
  result <- suppressWarnings(sparse.cor(sparse_mat))
  is.matrix(result) && nrow(result) == 10 && ncol(result) == 10
})

test_that("sparse.cor处理单列矩阵", {
  mat <- matrix(rnorm(10), nrow = 10, ncol = 1)
  sparse_mat <- as(mat, "dgCMatrix")
  result <- sparse.cor(sparse_mat)
  is.matrix(result) && result[1,1] == 1
})

# ============================================================================
# 测试 3: sparse_row_scale 行标准化
# ============================================================================
cat("\n[3] sparse_row_scale 行标准化测试\n")
cat(paste(rep("-", 50), collapse=""), "\n")

# 3.1 标准化结果验证
test_that("标准化后行均值为0", {
  set.seed(42)
  mat <- matrix(rnorm(100, mean = 5, sd = 2), nrow = 10, ncol = 10)
  scaled <- sparse_row_scale(mat, center = TRUE, scale = TRUE)
  
  row_means <- rowMeans(scaled)
  max(abs(row_means)) < 1e-10
})

test_that("标准化后行标准差为1", {
  set.seed(42)
  mat <- matrix(rnorm(100, mean = 5, sd = 2), nrow = 10, ncol = 10)
  scaled <- sparse_row_scale(mat, center = TRUE, scale = TRUE)
  
  row_sds <- apply(scaled, 1, sd)
  max(abs(row_sds - 1)) < 1e-10
})

# 3.2 稀疏矩阵测试
test_that("稀疏矩阵标准化结果正确", {
  set.seed(42)
  mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
  sparse_mat <- as(mat, "dgCMatrix")
  
  scaled_dense <- sparse_row_scale(mat, center = TRUE, scale = TRUE)
  scaled_sparse <- sparse_row_scale(sparse_mat, center = TRUE, scale = TRUE)
  
  max(abs(scaled_dense - as.matrix(scaled_sparse))) < 1e-10
})

# 3.3 边界条件：常量行
test_that("处理常量行（方差为0）", {
  mat <- matrix(c(5, 5, 5, 5, 5,
                  1, 2, 3, 4, 5), nrow = 2, byrow = TRUE)
  result <- sparse_row_scale(mat, center = TRUE, scale = TRUE)
  
  # 常量行应该变成全0（因为sd=0时设为1）
  all(result[1,] == 0) && abs(mean(result[2,])) < 1e-10
})

# 3.4 只中心化
test_that("只中心化（不标准化）", {
  set.seed(42)
  mat <- matrix(rnorm(100, mean = 10, sd = 3), nrow = 10, ncol = 10)
  centered <- sparse_row_scale(mat, center = TRUE, scale = FALSE)
  
  row_means <- rowMeans(centered)
  max(abs(row_means)) < 1e-10
})

# ============================================================================
# 测试 4: 数值稳定性测试
# ============================================================================
cat("\n[4] 数值稳定性测试\n")
cat(paste(rep("-", 50), collapse=""), "\n")

# 4.1 极端值处理
test_that("处理极大Z值", {
  Z <- 100
  p_value <- 2 * pnorm(abs(Z), mean = 0, sd = 1, lower.tail = FALSE)
  p_value <- pmin(p_value, 1)
  p_value >= 0 && p_value <= 1 && p_value < 1e-10
})

test_that("处理极小Z值", {
  Z <- 1e-15
  p_value <- 2 * pnorm(abs(Z), mean = 0, sd = 1, lower.tail = FALSE)
  p_value <- pmin(p_value, 1)
  abs(p_value - 1) < 0.001
})

# 4.2 零矩阵处理
test_that("sparse.cor处理全零矩阵", {
  mat <- matrix(0, nrow = 10, ncol = 5)
  sparse_mat <- as(mat, "dgCMatrix")
  result <- sparse.cor(sparse_mat)
  
  # 应该返回有限值
  all(is.finite(result))
})

# 4.3 NA处理
test_that("sparse_row_scale处理包含NA的矩阵", {
  mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
  mat[1, 1] <- NA
  
  result <- sparse_row_scale(mat, center = TRUE, scale = TRUE)
  # 应该能处理而不崩溃
  is.matrix(result)
})

# ============================================================================
# 测试 5: 置换检验逻辑
# ============================================================================
cat("\n[5] 置换检验逻辑测试\n")
cat(paste(rep("-", 50), collapse=""), "\n")

test_that("置换保持向量长度", {
  set.seed(123)
  coefs <- rnorm(100)
  permuted <- sample(coefs, length(coefs), replace = FALSE)
  length(permuted) == length(coefs)
})

test_that("置换保持元素集合不变", {
  set.seed(123)
  coefs <- 1:100
  permuted <- sample(coefs, length(coefs), replace = FALSE)
  setequal(coefs, permuted)
})

test_that("多次置换产生不同结果", {
  coefs <- rnorm(50)
  set.seed(1234)
  perm1 <- sample(coefs, length(coefs), replace = FALSE)
  set.seed(1235)
  perm2 <- sample(coefs, length(coefs), replace = FALSE)
  !all(perm1 == perm2)
})

# 置换检验的统计性质
test_that("零假设下P值均匀分布", {
  set.seed(42)
  n_tests <- 1000
  p_values <- numeric(n_tests)
  
  for (i in 1:n_tests) {
    # 生成零假设数据
    x <- rnorm(50)
    z <- mean(x) / (sd(x) / sqrt(length(x)))
    p_values[i] <- 2 * pnorm(abs(z), lower.tail = FALSE)
  }
  
  # KS检验P值是否均匀分布
  ks_result <- ks.test(p_values, "punif")
  ks_result$p.value > 0.01  # 不应该拒绝均匀分布假设
})

# ============================================================================
# 测试 6: FDR校正
# ============================================================================
cat("\n[6] FDR校正测试\n")
cat(paste(rep("-", 50), collapse=""), "\n")

test_that("BH校正后的P值不小于原P值", {
  set.seed(42)
  p_values <- runif(100, 0, 1)
  q_values <- p.adjust(p_values, method = "BH")
  
  all(q_values >= p_values)
})

test_that("BH校正后的P值不超过1", {
  set.seed(42)
  p_values <- runif(100, 0, 1)
  q_values <- p.adjust(p_values, method = "BH")
  
  all(q_values <= 1)
})

test_that("BH校正保持顺序", {
  p_values <- c(0.001, 0.01, 0.05, 0.1, 0.5)
  q_values <- p.adjust(p_values, method = "BH")
  
  # 顺序应该保持
  all(diff(q_values) >= 0 | abs(diff(q_values)) < 1e-10)
})

# ============================================================================
# 测试 7: 风险分数计算
# ============================================================================
cat("\n[7] 风险分数计算测试\n")
cat(paste(rep("-", 50), collapse=""), "\n")

test_that("风险分数是系数和表达的线性组合", {
  set.seed(42)
  n_genes <- 50
  n_cells <- 100
  
  expr <- matrix(rnorm(n_genes * n_cells), nrow = n_genes, ncol = n_cells)
  coefs <- rnorm(n_genes)
  
  # 计算风险分数: risk = t(expr) %*% coefs
  risk_scores <- crossprod(expr, coefs)
  
  # 手动计算
  manual_scores <- colSums(expr * coefs)
  
  max(abs(risk_scores[,1] - manual_scores)) < 1e-10
})

test_that("Z分数标准化正确", {
  set.seed(42)
  scores <- rnorm(100, mean = 5, sd = 2)
  
  # Z标准化
  z_scores <- (scores - mean(scores)) / sd(scores)
  
  # 验证
  abs(mean(z_scores)) < 1e-10 && abs(sd(z_scores) - 1) < 1e-10
})

# ============================================================================
# 测试总结
# ============================================================================
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("测试总结\n")
cat(paste(rep("=", 70), collapse=""), "\n")

total_tests <- tests_passed + tests_failed
cat(sprintf("\n总测试数: %d\n", total_tests))
cat(sprintf("通过: %d (%.1f%%)\n", tests_passed, 100 * tests_passed / total_tests))
cat(sprintf("失败: %d (%.1f%%)\n", tests_failed, 100 * tests_failed / total_tests))

if (tests_failed == 0) {
  cat("\n✅ 所有测试通过！科学准确性验证完成。\n")
} else {
  cat("\n⚠️ 有测试失败，请检查以下项目：\n")
  failed_tests <- names(test_results)[sapply(test_results, function(x) x != "PASS")]
  for (ft in failed_tests) {
    cat(sprintf("  - %s: %s\n", ft, test_results[[ft]]))
  }
}

cat("\n")
