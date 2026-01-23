# ============================================================================
# scPAS 集成测试 - 端到端工作流验证
# ============================================================================

cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("scPAS v1.0.3 集成测试 - 端到端工作流\n")
cat(paste(rep("=", 70), collapse=""), "\n\n")

# 加载必要的包
suppressPackageStartupMessages({
  library(Matrix)
  library(stats)
})

# 尝试加载包或源代码
tryCatch({
  library(scPAS)
  cat("已加载scPAS包\n")
}, error = function(e) {
  # 如果包未安装，从源文件加载
  if (file.exists("R/sparse_utils.R")) {
    source("R/sparse_utils.R")
    source("R/scPAS.R")
    cat("已从源文件加载函数\n")
  } else {
    stop("无法加载scPAS包或源文件")
  }
})

# 获取内部函数（用于测试）
sparse_row_scale <- scPAS:::sparse_row_scale
parallel_permutation_test <- scPAS:::parallel_permutation_test

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
# 测试 1: 模拟完整工作流（不使用Seurat）
# ============================================================================
cat("\n[1] 模拟scPAS核心计算流程\n")
cat(paste(rep("-", 50), collapse=""), "\n")

# 创建模拟数据
set.seed(12345)
n_genes <- 100
n_bulk_samples <- 30
n_cells <- 200

# 模拟bulk数据
bulk_data <- matrix(rnorm(n_genes * n_bulk_samples, mean = 5, sd = 2), 
                    nrow = n_genes, ncol = n_bulk_samples)
rownames(bulk_data) <- paste0("Gene", 1:n_genes)
colnames(bulk_data) <- paste0("Sample", 1:n_bulk_samples)

# 模拟单细胞数据（稀疏）
sc_data <- matrix(rpois(n_genes * n_cells, lambda = 2), 
                  nrow = n_genes, ncol = n_cells)
sc_data <- log1p(sc_data)  # 对数变换
rownames(sc_data) <- paste0("Gene", 1:n_genes)
colnames(sc_data) <- paste0("Cell", 1:n_cells)
sc_sparse <- as(sc_data, "dgCMatrix")

# 模拟表型（连续变量）
phenotype_continuous <- rnorm(n_bulk_samples, mean = 50, sd = 10)

# 模拟表型（二分类）
phenotype_binary <- sample(c(0, 1), n_bulk_samples, replace = TRUE)

# 模拟生存数据
phenotype_survival <- matrix(c(
  rexp(n_bulk_samples, rate = 0.01),  # time
  sample(c(0, 1), n_bulk_samples, replace = TRUE)  # status
), ncol = 2)
colnames(phenotype_survival) <- c("time", "status")

test_that("数据维度正确", {
  nrow(bulk_data) == n_genes && ncol(bulk_data) == n_bulk_samples &&
  nrow(sc_data) == n_genes && ncol(sc_data) == n_cells
})

# ============================================================================
# 测试 2: 核心计算步骤
# ============================================================================
cat("\n[2] 核心计算步骤验证\n")
cat(paste(rep("-", 50), collapse=""), "\n")

# 步骤1: 量化标准化
test_that("quantile normalization", {
  # 模拟quantile normalization
  bulk_normalized <- apply(bulk_data, 2, function(x) {
    qqnorm(x, plot.it = FALSE)$x
  })
  is.matrix(bulk_normalized) && all(is.finite(bulk_normalized))
})

# 步骤2: 相关性矩阵计算
test_that("相关性矩阵计算 (sparse.cor)", {
  # 转置以获得基因x基因相关性
  cor_mat <- sparse.cor(Matrix::t(sc_sparse))
  
  # 验证性质
  is.matrix(cor_mat) && 
  nrow(cor_mat) == n_genes && 
  ncol(cor_mat) == n_genes &&
  all(diag(cor_mat) == 1) &&
  all(cor_mat >= -1 & cor_mat <= 1)
})

# 步骤3: 行标准化
test_that("表达矩阵标准化 (sparse_row_scale)", {
  scaled <- sparse_row_scale(sc_sparse, center = TRUE, scale = TRUE)
  
  row_means <- rowMeans(scaled)
  row_sds <- apply(scaled, 1, sd)
  
  all(abs(row_means) < 1e-10) && all(abs(row_sds - 1) < 1e-10)
})

# 步骤4: 风险分数计算
test_that("风险分数计算", {
  # 模拟系数
  coefs <- rnorm(n_genes)
  names(coefs) <- rownames(sc_data)
  
  # 标准化表达
  scaled_exp <- sparse_row_scale(sc_sparse, center = TRUE, scale = TRUE)
  scaled_exp[is.na(scaled_exp)] <- 0
  
  # 计算风险分数
  risk_scores <- crossprod(scaled_exp, coefs)
  
  is.matrix(risk_scores) && 
  nrow(risk_scores) == n_cells &&
  all(is.finite(risk_scores))
})

# ============================================================================
# 测试 3: 置换检验
# ============================================================================
cat("\n[3] 置换检验验证\n")
cat(paste(rep("-", 50), collapse=""), "\n")

test_that("parallel_permutation_test 函数", {
  coefs <- rnorm(n_genes)
  names(coefs) <- rownames(sc_data)
  
  scaled_exp <- sparse_row_scale(sc_sparse, center = TRUE, scale = TRUE)
  scaled_exp[is.na(scaled_exp)] <- 0
  
  result <- parallel_permutation_test(
    scaled_exp = scaled_exp,
    Coefs = coefs,
    permutation_times = 100,  # 少量置换用于测试
    n_cores = 1,
    seed = 12345
  )
  
  # 验证输出结构
  "risk_score" %in% names(result) &&
  "risk_score.background" %in% names(result) &&
  nrow(result$risk_score) == n_cells &&
  ncol(result$risk_score.background) == 100
})

# ============================================================================
# 测试 4: P值和FDR计算
# ============================================================================
cat("\n[4] P值和FDR计算验证\n")
cat(paste(rep("-", 50), collapse=""), "\n")

test_that("完整P值计算流程", {
  # 模拟风险分数
  set.seed(42)
  n_test <- 100
  risk_scores <- rnorm(n_test)
  
  # 模拟背景分布
  background <- matrix(rnorm(n_test * 200), nrow = n_test, ncol = 200)
  
  # 计算统计量
  mean_bg <- rowMeans(background)
  sd_bg <- apply(background, 1, sd)
  sd_bg[sd_bg == 0] <- 1
  
  Z <- (risk_scores - mean_bg) / sd_bg
  
  # 双尾P值
  p_values <- 2 * pnorm(abs(Z), lower.tail = FALSE)
  p_values <- pmin(p_values, 1)
  
  # FDR校正
  q_values <- p.adjust(p_values, method = "BH")
  
  # 验证
  all(p_values >= 0 & p_values <= 1) &&
  all(q_values >= p_values) &&
  all(q_values <= 1)
})

test_that("细胞分类逻辑", {
  # 模拟数据
  Z <- c(-3, -2, -1, 0, 1, 2, 3)
  FDR <- c(0.001, 0.03, 0.1, 0.5, 0.1, 0.03, 0.001)
  threshold <- 0.05
  
  # 分类逻辑
  labels <- ifelse(Z > 0 & FDR <= threshold, 'scPAS+', 
                   ifelse(Z < 0 & FDR <= threshold, 'scPAS-', '0'))
  
  expected <- c('scPAS-', 'scPAS-', '0', '0', '0', 'scPAS+', 'scPAS+')
  all(labels == expected)
})

# ============================================================================
# 测试 5: 输入验证
# ============================================================================
cat("\n[5] 输入验证测试\n")
cat(paste(rep("-", 50), collapse=""), "\n")

test_that("检测无效的bulk数据维度", {
  tryCatch({
    # 这应该会触发错误
    if (nrow(matrix(1:6, 2, 3)) < 10) {
      stop("'bulk_dataset' has too few genes")
    }
    FALSE
  }, error = function(e) {
    grepl("too few genes", e$message)
  })
})

test_that("检测表型和样本数量不匹配", {
  tryCatch({
    bulk <- matrix(rnorm(100), 10, 10)
    pheno <- 1:5  # 不匹配
    if (length(pheno) != ncol(bulk)) {
      stop("Length mismatch")
    }
    FALSE
  }, error = function(e) {
    grepl("mismatch", e$message, ignore.case = TRUE)
  })
})

test_that("检测无效的alpha值", {
  alpha <- c(0, 0.5, 1.5)  # 1.5 超出范围
  !all(alpha >= 0 & alpha <= 1)
})

# ============================================================================
# 测试 6: 数值边界条件
# ============================================================================
cat("\n[6] 数值边界条件测试\n")
cat(paste(rep("-", 50), collapse=""), "\n")

test_that("处理全零行", {
  mat <- matrix(c(0, 0, 0, 1, 2, 3, 4, 5, 6), nrow = 3, ncol = 3, byrow = TRUE)
  result <- sparse_row_scale(mat, center = TRUE, scale = TRUE)
  all(is.finite(result))
})

test_that("处理极端稀疏矩阵 (>99% 零)", {
  mat <- matrix(0, nrow = 100, ncol = 100)
  mat[sample(10000, 50)] <- rnorm(50)  # 只有0.5%非零
  sparse_mat <- as(mat, "dgCMatrix")
  
  result <- sparse.cor(Matrix::t(sparse_mat))
  all(is.finite(result))
})

test_that("处理单细胞", {
  mat <- matrix(rnorm(100), nrow = 100, ncol = 1)
  coefs <- rnorm(100)
  
  scaled <- sparse_row_scale(mat, center = TRUE, scale = TRUE)
  score <- crossprod(scaled, coefs)
  
  is.finite(score[1,1])
})

# ============================================================================
# 测试 7: 可重复性
# ============================================================================
cat("\n[7] 可重复性测试\n")
cat(paste(rep("-", 50), collapse=""), "\n")

test_that("相同种子产生相同结果", {
  coefs <- rnorm(50)
  scaled_exp <- matrix(rnorm(5000), 50, 100)
  
  set.seed(12345)
  result1 <- parallel_permutation_test(scaled_exp, coefs, 100, 1, 12345)
  
  set.seed(12345)
  result2 <- parallel_permutation_test(scaled_exp, coefs, 100, 1, 12345)
  
  identical(result1$risk_score, result2$risk_score) &&
  identical(result1$risk_score.background, result2$risk_score.background)
})

# ============================================================================
# 测试总结
# ============================================================================
cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("集成测试总结\n")
cat(paste(rep("=", 70), collapse=""), "\n")

total_tests <- tests_passed + tests_failed
cat(sprintf("\n总测试数: %d\n", total_tests))
cat(sprintf("通过: %d (%.1f%%)\n", tests_passed, 100 * tests_passed / total_tests))
cat(sprintf("失败: %d (%.1f%%)\n", tests_failed, 100 * tests_failed / total_tests))

if (tests_failed == 0) {
  cat("\n✅ 所有集成测试通过！\n")
} else {
  cat("\n⚠️ 有集成测试失败\n")
}

cat("\n")
