# ============ Internal helper functions for Seurat 4/5 compatibility ============
# Detect Seurat version once at load time for efficiency
.getSeuratVersion <- function() {

  utils::packageVersion("SeuratObject")
}

# GetAssayData: Seurat 5 uses 'layer' instead of 'slot'
.getAssayData <- function(object, assay, slot_name = 'data') {
  if (.getSeuratVersion() >= "5.0.0") {
    Seurat::GetAssayData(object = object, assay = assay, layer = slot_name)
  } else {
    Seurat::GetAssayData(object = object, assay = assay, slot = slot_name)
  }
}

# Get graphs from Seurat object (compatible with both v4 and v5)
.getGraphs <- function(object, name) {
  # Both Seurat 4 and 5 support Graphs() function from SeuratObject
  graphs <- SeuratObject::Graphs(object)
  if (name %in% names(graphs)) {
    return(object[[name]])
  }
  # Fallback to direct slot access if needed
  if (!is.null(object@graphs) && name %in% names(object@graphs)) {
    return(object@graphs[[name]])
  }
  stop(paste0("Graph '", name, "' not found. Please run FindNeighbors() first."))
}

# Get misc data from Seurat object
.getMisc <- function(object, name = NULL) {
  misc <- Seurat::Misc(object)
  if (is.null(name)) {
    return(misc)
  }
  return(misc[[name]])
}

# Set misc data in Seurat object
.setMisc <- function(object, name, value) {
  Seurat::Misc(object, slot = name) <- value
  return(object)
}

# Create Assay object (compatible with both v4 and v5)
.createAssayObject <- function(data) {
  if (.getSeuratVersion() >= "5.0.0") {
    # Seurat 5: CreateAssayObject still works but we use it carefully
    # The 'data' parameter in v5 is interpreted as normalized data
    Seurat::CreateAssayObject(data = data)
  } else {
    # Seurat 4: standard usage
    Seurat::CreateAssayObject(data = data)
  }
}

#' scPAS : A tool for identifying Phenotype-Associated cell Subpopulations from single-cell sequencing data by integrating bulk data
#'
#' @param bulk_dataset Matrix. Bulk expression matrix of related disease. Each row represents a gene and each column represents a sample. The input expression values are continuous, such as microarray fluorescent units in logarithmic scale, RNA-seq log-CPMs, log-RPKMs or log-TPMs.
#' @param sc_dataset Matrix or seurat object. Single-cell RNA-seq expression matrix of related disease. Each row represents a gene and each column represents a sample. A Seurat object that contains the preprocessed data and constructed network is preferred. Otherwise, a cell-cell similarity network is constructed based on the input matrix.Otherwise, the raw count expression matrix will be processed by using Seurat's default parameters. See run_Seurat for details.
#' @param phenotype Phenotype annotation of each bulk sample. It can be a continuous dependent variable,
#' binary group indicator vector, or clinical survival data:
#'   \itemize{
#'   \item Continuous dependent variable. Should be a quantitative vector for \code{family = gaussian}.
#'   \item Binary group indicator vector. Should be either a 0-1 encoded vector or a factor with two levels for \code{family = binomial}.
#'   \item Clinical survival data. Should be a two-column matrix with columns named 'time' and 'status'. The latter is a binary variable,
#'   with '1' indicating event (e.g.recurrence of cancer or death), and '0' indicating right censored.
#'   The function \code{Surv()} in package survival produces such a matrix.
#'   }
#' @param assay Name of Assay to get.
#' @param tag Names for each phenotypic group. Used for logistic regressions only.
#' @param nfeature Numeric. The Number of features to select as top variable features in sc_dataset. Top variable features will be used to intersect with the features of bulk_dataset. Default is NULL.All features will be used.
#' @param do_imputation Logical. Whether to perform imputation on single-cell data (default: TRUE).
#' @param imputation_method Character. Name of alternative method for imputation.
#' @param alpha Numeric. Parameter used to balance the effect of the l1 norm and the network-based penalties. It can be a number or a searching vector.
#' If \code{alpha = NULL}, a default searching vector is used. The range of alpha is in \code{[0,1]}. A larger alpha lays more emphasis on the l1 norm.
#' @param network_class  The source of feature-feature similarity network. By default this is set to \code{sc} and the other one is \code{bulk}.
#' @param family Character. Response type for the regression model.  It depends on the type of the given phenotype and
#' can be \code{family = gaussian} for linear regression, \code{family = binomial} for classification, or \code{family = cox} for Cox regression.
#' @param permutation_times Integer. Number of permutation iterations for statistical
#'   significance testing (default: 2000). Higher values increase accuracy but also
#'   computation time. Recommended: 1000-5000. For faster testing, use 500-1000.
#' @param FDR.threshold Numeric. FDR value threshold for identifying phenotype-associated cells.
#'   The default is 0.05.
#' @param n_cores Integer. Number of CPU cores to use for parallel permutation test
#'   (default: 1 for sequential processing). Setting n_cores > 1 enables parallel
#'   computing which can significantly speed up the analysis (2-4x faster with 4 cores).
#'   Requires 'future' and 'future.apply' packages. If these packages are not available,
#'   the function will automatically fall back to sequential processing.
#' @param independent Logical. The background distribution of risk scores is constructed
#'   independently of each cell.
#'
#' @return This function returns a Seurat object with the following components added to:
#'   \item{scPAS_para}{A list contains the final model parameters added to misc.}
#'   \item{PAS result}{A data frame containing risk scores (scPAS_RS), normalized risk
#'     scores (scPAS_NRS), p-value (scPAS_Pvalue), adjusted p-value (scPAS_FDR), and
#'     cell classification labels (scPAS) added to metaData.}
#'
#' @import Seurat Matrix preprocessCore
#'
#' @export
scPAS <- function(bulk_dataset, sc_dataset, phenotype, assay = 'RNA', tag = NULL, 
                  nfeature = NULL, do_imputation = TRUE, imputation_method = c('KNN','ALRA'),
                  alpha = NULL, network_class = c('SC','bulk'), independent = TRUE, 
                  family = c("gaussian","binomial","cox"), permutation_times = 2000,
                    FDR.threshold = 0.05, n_cores = 1){
  
  # ============ Input Validation ============
  # Validate bulk_dataset
  if (missing(bulk_dataset) || is.null(bulk_dataset)) {
    stop("'bulk_dataset' is required.")
  }
  if (!is.matrix(bulk_dataset) && !methods::is(bulk_dataset, "Matrix")) {
    bulk_dataset <- as.matrix(bulk_dataset)
  }
  if (nrow(bulk_dataset) < 10) {
    stop("'bulk_dataset' has too few genes (< 10). Please check input data.")
  }
  if (ncol(bulk_dataset) < 3) {
    stop("'bulk_dataset' has too few samples (< 3). At least 3 samples are required.")
  }
  
  # Validate sc_dataset
  if (missing(sc_dataset) || is.null(sc_dataset)) {
    stop("'sc_dataset' is required.")
  }
  
  # Validate phenotype
  if (missing(phenotype) || is.null(phenotype)) {
    stop("'phenotype' is required.")
  }
  
  # Match arguments
  network_class <- match.arg(network_class)
  family <- match.arg(family)
  imputation_method <- match.arg(imputation_method)
  
  # Validate phenotype based on family
  if (family == "gaussian") {
    if (!is.numeric(phenotype)) {
      stop("For family='gaussian', phenotype must be a numeric vector.")
    }
    if (length(phenotype) != ncol(bulk_dataset)) {
      stop(sprintf("Length of phenotype (%d) must match number of bulk samples (%d).", 
                   length(phenotype), ncol(bulk_dataset)))
    }
  } else if (family == "binomial") {
    if (is.factor(phenotype)) {
      if (nlevels(phenotype) != 2) {
        stop("For family='binomial', phenotype factor must have exactly 2 levels.")
      }
      phenotype <- as.numeric(phenotype) - 1
    } else {
      phenotype <- as.numeric(phenotype)
      if (!all(phenotype %in% c(0, 1))) {
        stop("For family='binomial', phenotype must be 0/1 or a two-level factor.")
      }
    }
    if (length(phenotype) != ncol(bulk_dataset)) {
      stop(sprintf("Length of phenotype (%d) must match number of bulk samples (%d).", 
                   length(phenotype), ncol(bulk_dataset)))
    }
  } else if (family == "cox") {
    if (!is.matrix(phenotype) || ncol(phenotype) != 2) {
      stop("For family='cox', phenotype must be a 2-column matrix (time, status). Use survival::Surv().")
    }
    if (!all(c("time", "status") %in% colnames(phenotype))) {
      colnames(phenotype) <- c("time", "status")
    }
    if (nrow(phenotype) != ncol(bulk_dataset)) {
      stop(sprintf("Number of phenotype rows (%d) must match number of bulk samples (%d).", 
                   nrow(phenotype), ncol(bulk_dataset)))
    }
  }
  
  # Validate numeric parameters
  if (!is.numeric(permutation_times) || permutation_times < 100) {
    stop("'permutation_times' must be a numeric value >= 100.")
  }
  if (!is.numeric(FDR.threshold) || FDR.threshold <= 0 || FDR.threshold >= 1) {
    stop("'FDR.threshold' must be between 0 and 1 (exclusive).")
  }
  if (!is.numeric(n_cores) || n_cores < 1) {
    stop("'n_cores' must be a positive integer.")
  }
  n_cores <- as.integer(n_cores)
  
  # Validate alpha if provided
  if (!is.null(alpha)) {
    if (!is.numeric(alpha) || any(alpha < 0) || any(alpha > 1)) {
      stop("'alpha' must be numeric values between 0 and 1.")
    }
  }
  
  # ============ Main Processing ============
  # Set default assay for Seurat object
  if (inherits(sc_dataset, 'Seurat')) {
  Seurat::DefaultAssay(sc_dataset) <- assay
  }

  if(inherits(sc_dataset, 'Seurat')){

    if(is.null(nfeature)){
      common <- intersect(rownames(bulk_dataset), rownames(sc_dataset))
    }
    else if(is.numeric(nfeature) & length(nfeature)==1){
      sc_dataset <- Seurat::FindVariableFeatures(sc_dataset, selection.method = "vst", verbose = FALSE, nfeatures = nfeature)

      common <- intersect(rownames(bulk_dataset), Seurat::VariableFeatures(sc_dataset))
      common <- common[!grepl(pattern ='^RP[LS]',common)]
      common <- common[!grepl(pattern ='^MT-',common)]
    }else if(is.character(nfeature) & length(nfeature)>1){
      common <- intersect(rownames(bulk_dataset), rownames(sc_dataset))
      common <- intersect(common, nfeature)
      common <- common[!grepl(pattern = '^RP[LS]', common)]
      common <- common[!grepl(pattern = '^MT-', common)]
    }
  }else{
    message("Step 0: The single-cell data is not a Seurat object, and a default Seurat pipeline will be run.")
    sc_dataset <- run_Seurat(sc_dataset)
    if(is.null(nfeature)){
      common <- intersect(rownames(bulk_dataset), rownames(sc_dataset))
      common <- common[!grepl(pattern = '^RP[LS]', common)]
      common <- common[!grepl(pattern = '^MT-', common)]
    }
    else if(is.numeric(nfeature)){
      sc_dataset <- Seurat::FindVariableFeatures(sc_dataset, selection.method = "vst", verbose = FALSE, nfeatures = nfeature)
      common <- intersect(rownames(bulk_dataset), Seurat::VariableFeatures(sc_dataset))
    }else if(is.character(nfeature) & length(nfeature)>1){
      common <- intersect(rownames(bulk_dataset), rownames(sc_dataset))
      common <- intersect(common, nfeature)
      common <- common[!grepl(pattern = '^RP[LS]', common)]
      common <- common[!grepl(pattern = '^MT-', common)]
    }
  }
  if (length(common) == 0) {
    stop("There are no common genes between the given single-cell and bulk samples.")
  }

  message("Step 1: Quantile normalization of bulk data.")
  Expression_bulk <- preprocessCore::normalize.quantiles(as.matrix(bulk_dataset))
  rownames(Expression_bulk) <- rownames(bulk_dataset)
  colnames(Expression_bulk) <- colnames(bulk_dataset)
  Expression_bulk <- Expression_bulk[common,]


  if(do_imputation){
    sc_dataset <- imputation(sc_dataset, assay = assay, method = imputation_method)
    assay <- Seurat::DefaultAssay(sc_dataset)
  }

  message("Step 2: Extracting single-cell expression profiles....")
  sc_exprs <- .getAssayData(object = sc_dataset, assay = assay, slot_name = 'data')
  #Expression_cell <- as(preprocessCore::normalize.quantiles(as.matrix(sc_exprs)), "dgCMatrix")
  Expression_cell <- sc_exprs
  rownames(Expression_cell) <- rownames(sc_exprs)
  colnames(Expression_cell) <- colnames(sc_exprs)
  Expression_cell <- Expression_cell[common,]
  sc_exprs <- NULL
  bulk_dataset <- NULL

  x <- t(Expression_bulk)


  # Construct a gene network
  if(network_class == 'bulk'){
    message("Step 3: Constructing a gene-gene similarity by bulk data....")
    cor.m <- stats::cor(x)
  }else{
    message("Step 3: Constructing a gene-gene similarity by single cell data....")
    # Ensure sparse matrix format for efficient transpose
    if (!methods::is(Expression_cell, "sparseMatrix")) {
      Expression_cell <- methods::as(Expression_cell, "sparseMatrix")
    }
    Expression_cell_t <- Matrix::t(Expression_cell)
    cor.m <- sparse.cor(Expression_cell_t)
  }
  # Set negative correlations to zero (only positive correlations are used)
  cor.m[cor.m < 0] <- 0
  # Ensure row/column names for FindNeighbors
  if (is.null(rownames(cor.m))) {
    rownames(cor.m) <- colnames(cor.m) <- rownames(Expression_cell)
  }
  SNN <- Seurat::FindNeighbors(1 - cor.m, distance.matrix = TRUE)
  Network <- as.matrix(SNN$snn)
  diag(Network) <- 0
  # Replace NA values and apply threshold
  Network[is.na(Network)] <- 0
  Network[Network > 0.2] <- 1
  Network[Network <= 0.2] <- 0

  message("Step 4: Optimizing the network-regularized sparse regression model....")
  if (family == "binomial"){
    y <- as.numeric(phenotype)
    z <- table(y)
    if (!is.null(tag) && length(tag) >= 2) {
    message(sprintf("Current phenotype contains %d %s and %d %s samples.", z[1], tag[1], z[2], tag[2]))
    } else {
      message(sprintf("Current phenotype contains %d class-0 and %d class-1 samples.", z[1], z[2]))
    }
    message("Performing logistic regression on the given phenotypes...")
  }
  if (family == "gaussian"){
    y <- as.numeric(phenotype)
    message("Performing linear regression on the given phenotypes...")
  }
  if (family == "cox"){
    y <- as.matrix(phenotype)
    if (ncol(y) != 2){
      stop("The size of survival data is wrong. Please check inputs and selected regression type.")
    }else{
      message("Performing cox regression on the given clinical outcomes...")
    }
  }

  if (is.null(alpha)){
    alpha <- c(0.001,0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

  }
  lambda <- c()
  for (i in 1:length(alpha)){
    set.seed(123)
    fit0 <- APML0(x = x, y = y, family = family, penalty = 'Net', Omega = Network, alpha = alpha[i], nlambda = 100, nfolds = min(10, nrow(x)))
    fit1 <- APML0(x = x, y = y, family = family, penalty = 'Net', Omega = Network, alpha = alpha[i], lambda = fit0$lambda.min)
    lambda <- c(lambda, fit0$lambda.min)
    if (family == "binomial"){
      Coefs <- as.numeric(fit1$Beta[2:(ncol(x)+1)])
    }else{
      Coefs <- as.numeric(fit1$Beta)
    }
    Feature1 <- colnames(x)[which(Coefs > 0)]
    Feature2 <- colnames(x)[which(Coefs < 0)]
    percentage <- (length(Feature1) + length(Feature2)) / ncol(x)
    message(sprintf("alpha = %s", alpha[i]))
    message(sprintf("lambda = %s", fit0$lambda.min))
    message(sprintf("scPAS identified %d risk+ features and %d risk- features.", length(Feature1), length(Feature2)))
    message(sprintf("The percentage of selected features is: %s%%", formatC(percentage*100, format = 'f', digits = 3)))
    cat("\n")
  }
  message("|**************************************************|")

  message("Step 5: Calculating quantified risk scores....")
  names(Coefs) <- colnames(x)
  
  # Scale expression data by row (gene)
  scaled_exp <- sparse_row_scale(Expression_cell, center = TRUE, scale = TRUE)
  colnames(scaled_exp) <- colnames(Expression_cell)
  rownames(scaled_exp) <- rownames(Expression_cell)
  # Replace NA values with 0
  scaled_exp[is.na(scaled_exp)] <- 0
  
  # Permutation test for significance
  if (n_cores > 1) {
    message(paste0("Step 6: Qualitative identification by permutation test (", 
                  permutation_times, " permutations, ", n_cores, " cores)"))
  } else {
    message(paste0("Step 6: Qualitative identification by permutation test (", 
                  permutation_times, " permutations, sequential)"))
  }
  
  perm_results <- parallel_permutation_test(
    scaled_exp = scaled_exp,
    Coefs = Coefs,
    permutation_times = permutation_times,
    n_cores = n_cores,
    seed = 12345
  )
  
  risk_score <- perm_results$risk_score
  risk_score.background <- perm_results$risk_score.background
  
  # Ensure matrix format for background distribution
  if (!is.matrix(risk_score.background)) {
    risk_score.background <- as.matrix(risk_score.background)
  }
  
  if(independent){
    mean.background <- Matrix::rowMeans(risk_score.background)
    sd.background <- apply(risk_score.background, 1, stats::sd)
    # Handle zero standard deviation (constant background)
    sd.background[sd.background == 0 | !is.finite(sd.background)] <- 1
  }else{
    mean.background <- mean(as.matrix(risk_score.background))
    sd.background <- stats::sd(as.matrix(risk_score.background))
    # Handle zero standard deviation
    if (sd.background == 0 || !is.finite(sd.background)) {
      sd.background <- 1
    }
  }

  Z <- (risk_score[,1] - mean.background) / sd.background
  # Handle non-finite Z values
  Z[!is.finite(Z)] <- 0

  # Two-tailed p-value: multiply by 2 for proper statistical testing
  # This tests for both positive and negative deviations from null
  p.value <- 2 * stats::pnorm(q = abs(Z), mean = 0, sd = 1, lower.tail = FALSE)
  # Cap p-value at 1 (can exceed 1 due to floating point in extreme cases)
  p.value <- pmin(p.value, 1)
  q.value <- stats::p.adjust(p = p.value, method = 'BH')
  risk_score_data.frame <- data.frame(
    cell = colnames(Expression_cell),
    raw_score = risk_score[,1],
    Z.statistics = Z,
    p.value = p.value,
    FDR = q.value,
    stringsAsFactors = FALSE
  )

  risk_score_data.frame$cell_label <- ifelse(Z > 0 & q.value <= FDR.threshold, 'scPAS+', 
                                              ifelse(Z < 0 & q.value <= FDR.threshold, 'scPAS-', '0'))

  sc_dataset <- .setMisc(sc_dataset, "scPAS_para", list(
    alpha = alpha, 
    lambda = lambda, 
    family = family,
    Coefs = Coefs,
    bulk = x,
    phenotype = y,
    Network = Network
  ))

  sc_dataset <- Seurat::AddMetaData(sc_dataset, metadata = risk_score_data.frame$raw_score, col.name = "scPAS_RS")
  sc_dataset <- Seurat::AddMetaData(sc_dataset, metadata = risk_score_data.frame$Z.statistics, col.name = "scPAS_NRS")
  sc_dataset <- Seurat::AddMetaData(sc_dataset, metadata = risk_score_data.frame$p.value, col.name = "scPAS_Pvalue")
  sc_dataset <- Seurat::AddMetaData(sc_dataset, metadata = risk_score_data.frame$FDR, col.name = "scPAS_FDR")
  sc_dataset <- Seurat::AddMetaData(sc_dataset, metadata = risk_score_data.frame$cell_label, col.name = "scPAS")

  message("Finished.")
  return(sc_dataset)

}

#' Preprocess the single-cell raw data using functions in the \code{Seurat} package
#'
#' This function provide a simplified-version of Seurat analysis pipeline for single-cell RNA-seq data. It contains the following steps in the pipeline:
#' \itemize{
#'    \item Create a \code{Seurat} object from raw data.
#'    \item Normalize the count data present in a given assay.
#'    \item Identify the variable features.
#'    \item Scales and centers features in the dataset.
#'    \item Run a PCA dimensionality reduction.
#'    \item Constructs a Shared Nearest Neighbor (SNN) Graph for a given dataset.
#'    \item Identify clusters of cells by a shared nearest neighbor (SNN) modularity optimization based clustering algorithm.
#'    \item Run t-distributed Stochastic Neighbor Embedding (t-SNE) dimensionality reduction on selected features.
#'    \item Runs the Uniform Manifold Approximation and Projection (UMAP) dimensional reduction technique.
#' }
#'
#' @param counts A \code{matrix}-like object with unnormalized data with cells as columns and features as rows.
#' @param project Project name for the \code{Seurat} object.
#' @param min.cells Include features detected in at least this many cells. Will subset the counts matrix as well.
#' To reintroduce excluded features, create a new object with a lower cutoff.
#' @param min.features Include cells where at least this many features are detected.
#' @param normalization.method Method for normalization.
#'   \itemize{
#'   \item LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor.
#'   This is then natural-log transformed using log1p.
#'   \item CLR: Applies a centered log ratio transformation.
#'   \item RC: Relative counts. Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor.
#'   No log-transformation is applied. For counts per million (CPM) set \code{scale.factor = 1e6}.
#' }
#' @param scale.factor Sets the scale factor for cell-level normalization.
#' @param selection.method How to choose top variable features. Choose one of :
#'   \itemize{
#'   \item vst: First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess).
#'   Then standardizes the feature values using the observed mean and expected variance (given by the fitted line).
#'   Feature variance is then calculated on the standardized values after clipping to a maximum (see clip.max parameter).
#'   \item mean.var.plot (mvp): First, uses a function to calculate average expression (mean.function) and dispersion (dispersion.function)
#'   for each feature. Next, divides features into num.bin (deafult 20) bins based on their average expression, and calculates
#'   z-scores for dispersion within each bin. The purpose of this is to identify variable features while controlling for the strong
#'   relationship between variability and average expression.
#'   \item dispersion (disp): selects the genes with the highest dispersion values
#'   }
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
#' @param dims_Neighbors Dimensions of reduction to use as input.
#' @param dims_TSNE Which dimensions to use as input features for t-SNE.
#' @param dims_UMAP Which dimensions to use as input features for UMAP.
#' @param meta.data meta data of single cell data.
#' @param verbose Print output.
#'
#' @return A \code{Seurat} object containing cell-cell similarity network, t-SNE and UMAP representations.
#' @import Seurat
#' @export
run_Seurat <- function(counts, project = "Single_Cell", min.cells = 400, min.features = 200,meta.data =NULL,
                       normalization.method = "LogNormalize", scale.factor = 10000,
                       selection.method = "vst", resolution = 0.6,
                       dims_Neighbors = 1:10, dims_TSNE = 1:10, dims_UMAP = 1:10,
                       verbose = TRUE){
  # Seurat functions available via Imports
  if(inherits(counts, c('matrix', 'dgCMatrix'))){
    data <- Seurat::CreateSeuratObject(counts = counts, project = project, min.cells = min.cells, min.features = min.features, meta.data = meta.data)
  }else if(inherits(counts, 'Seurat')){
    data <- counts
  }else{
    stop("The class of scRNA-seq data is wrong. Please input a count matrix or Seurat object")
  }
  data <- Seurat::NormalizeData(object = data, normalization.method = normalization.method, scale.factor = scale.factor, verbose = verbose)
  data <- Seurat::FindVariableFeatures(object = data, selection.method = selection.method, verbose = verbose)
  data <- Seurat::ScaleData(object = data, verbose = verbose)
  data <- Seurat::RunPCA(object = data, features = Seurat::VariableFeatures(data), verbose = verbose)
  data <- Seurat::FindNeighbors(object = data, dims = dims_Neighbors, verbose = verbose)
  data <- Seurat::FindClusters(object = data, resolution = resolution, verbose = verbose)
  data <- Seurat::RunTSNE(object = data, dims = dims_TSNE)
  data <- Seurat::RunUMAP(object = data, dims = dims_UMAP, verbose = verbose)

  return(data)
}

#' The function of imputaion.
#'
#' @param obj A seurat object.
#' @param assay The assay for imputation. The default is 'RNA'.
#' @param method The method for imputation. The default is 'RNA'.
#'
#' @return  A seurat object after imputaion.
#'
#' @import Matrix
#'
#'
#'
imputation <- function(obj,assay='RNA',method=c('KNN','ALRA')){
  method=match.arg(method)
  if(method=='KNN'){
    print("Step2: Imputation of missing values in single cell RNA-sequencing data with KNN....")
    obj <- imputation_KNN(obj = obj,assay = assay)
  }else if(method=='ALRA'){
    print("Step2: Imputation of missing values in single cell RNA-sequencing data with ALRA")
    obj <- imputation_ALRA(obj = obj,assay = assay)
  }else{
    warning(paste0("The '", method, "' method does not exist, so imputation is invalid!"))
  }
  return(obj)
}


#' A method for imputation of missing values in single cell RNA-sequencing data based on ALRA.
#'
#' @param obj A seurat object.
#' @param assay The assay for imputation. The default is 'RNA'.
#'
#' @return  A seurat object after imputaion.
#'
#' @import Matrix
#'
#'
#'
imputation_ALRA <- function(obj,assay='RNA'){
  # Check if ALRA is available
  if (!requireNamespace("ALRA", quietly = TRUE)) {
    stop("Package 'ALRA' is required for ALRA imputation. Please install it with: install.packages('ALRA')")
  }
  data <- .getAssayData(object = obj, assay = assay, slot_name = 'data')
  data_alra <- t(ALRA::alra(t(as.matrix(data)))[[3]])
  colnames(data_alra) <- colnames(data)
  data_alra <- Matrix::Matrix(data_alra, sparse = TRUE)

  obj[["imputation"]] <- .createAssayObject(data = data_alra)
  Seurat::DefaultAssay(obj) <- "imputation"
  return(obj)
}

#' A method for imputation of missing values in single cell RNA-sequencing data based on the average expression value of nearest neighbor cells.
#'
#' @param obj A seurat object.
#' @param assay The assay for imputation. The default is 'RNA'.
#' @param LogNormalized Whether the data is LogNormalized.
#'
#' @return A seurat object after imputaion.
#'
#' @import Matrix
#'
#'
#'
imputation_KNN <- function (obj,assay='RNA', LogNormalized = TRUE)
{
  # Matrix functions available via Imports
  exp_sc <- .getAssayData(object = obj, assay = assay, slot_name = 'data')
  nn_network <- .getGraphs(obj, paste0(assay, "_nn"))
  
  if (!methods::is(object = exp_sc, class2 = "sparseMatrix")) {
    exp_sc <- methods::as(exp_sc, "sparseMatrix")
  }
  if (!methods::is(object = nn_network, class2 = "sparseMatrix")) {
    nn_network <- methods::as(nn_network, "sparseMatrix")
  }
  if (LogNormalized) {
    exp_sc <- methods::as(exp(exp_sc) - 1, "sparseMatrix")
  }
  network_count <- methods::as(Matrix::Diagonal(x = 1/rowSums(nn_network)), "sparseMatrix")
  exp_sc_mean <- Matrix::tcrossprod(x = Matrix::tcrossprod(x = exp_sc, y = nn_network), y = network_count)
  
  if (LogNormalized) {
    exp_sc_mean <- log1p(exp_sc_mean)
  }
  colnames(exp_sc_mean) <- colnames(exp_sc)
  obj[["imputation"]] <- .createAssayObject(data = exp_sc_mean)
  Seurat::DefaultAssay(obj) <- "imputation"
  return(obj)
}


#' Compute correlation matrix for a sparse matrix
#'
#' @description
#' Computes the Pearson correlation matrix for a sparse matrix
#' without full dense conversion.
#'
#' @param x Matrix. Normalized single cell expression profile (cells x genes).
#'   Each row is a cell, each column is a gene.
#'
#' @return A correlation matrix (genes x genes).
#'
#' @details
#' This function handles sparse matrices by:
#' 1. Processing non-zero elements separately
#' 2. Handling numerical precision issues
#' 3. Ensuring proper correlation matrix properties (diagonal = 1, range [-1, 1])
#'
#' @keywords internal
sparse.cor <- function(x){
  # Ensure x is a sparse matrix
  if(!methods::is(x, "sparseMatrix")){
    x <- methods::as(x, "sparseMatrix")
  }
  
  n <- nrow(x)  # number of observations (cells)
  m <- ncol(x)  # number of variables (genes)
  
  # Handle edge cases
  if (n < 2) {
    warning("Less than 2 observations, returning identity matrix")
    return(diag(m))
  }
  if (m < 2) {
    return(matrix(1, 1, 1))
  }
  
  # Column means (gene means across cells)
  Ex <- as.vector(Matrix::colMeans(x))
  
  # Get indices of rows with at least one non-zero element
  ii <- unique(x@i) + 1  # 1-based indices
  
  # If all rows have non-zero elements, use standard correlation
  if (length(ii) == n) {
    # Convert to dense and use standard cor
    return(stats::cor(as.matrix(x)))
  }
  
  # Compute covariance matrix efficiently
  # Cov(X,Y) = E[XY] - E[X]E[Y]
  # For sparse matrices, we compute E[XY] from non-zero rows
  # and adjust for zero rows
  
  # Extract non-zero rows as dense matrix
  x_nonzero <- as.matrix(x[ii, , drop = FALSE])
  
  # Center the non-zero portion
  x_centered <- x_nonzero - matrix(Ex, nrow = length(ii), ncol = m, byrow = TRUE)
  
  # Contribution from non-zero rows
  covmat_nonzero <- crossprod(x_centered)
  
  # Contribution from zero rows (they contribute -Ex * -Ex = Ex^2)
  n_zero_rows <- n - length(ii)
  if (n_zero_rows > 0) {
    covmat_zero <- tcrossprod(Ex) * n_zero_rows
    covmat <- (covmat_nonzero + covmat_zero) / (n - 1)
  } else {
    covmat <- covmat_nonzero / (n - 1)
  }
  
  # Compute standard deviations
  sdvec <- sqrt(pmax(diag(covmat), 0))  # Ensure non-negative
  
  # Handle zero variance (constant columns)
  sdvec[sdvec == 0 | !is.finite(sdvec)] <- 1
  
  # Compute correlation matrix
  cormat <- covmat / tcrossprod(sdvec)
  
  # Ensure diagonal is exactly 1 and values are in [-1, 1]
  diag(cormat) <- 1
  cormat[cormat > 1] <- 1
  cormat[cormat < -1] <- -1
  
  # Handle any remaining NA/NaN values
  cormat[!is.finite(cormat)] <- 0
  
  return(cormat)
}



#' scPAS.prediction: A function that uses the scPAS model to make predictions on independent data
#'
#' @param model Seurat object. A Seurat object containing the scPAS model (from running scPAS()).
#' @param test.data Matrix or Seurat object. Single-cell RNA-seq expression matrix of related disease. 
#'   Each row represents a gene and each column represents a sample. A Seurat object that contains 
#'   the preprocessed data and constructed network is preferred.
#' @param assay Name of Assay to get.
#' @param FDR.threshold Numeric. FDR value threshold for identifying phenotype-associated cells.
#'   The default is 0.05.
#' @param do_imputation Logical. Whether to perform imputation on the test data (default: FALSE).
#' @param imputation_method Character. Imputation method: "KNN" or "ALRA".
#' @param independent Logical. Whether to compute background distribution independently for each cell.
#' @param permutation_times Integer. Number of permutations for significance testing (default: 2000).
#' @param n_cores Integer. Number of CPU cores for parallel processing (default: 1).
#'
#' @return A Seurat object or data frame containing the prediction results.
#'
#' @export
scPAS.prediction <- function(model, test.data, assay = 'RNA', FDR.threshold = 0.05, 
                             do_imputation = FALSE, imputation_method = 'KNN', 
                             independent = TRUE, permutation_times = 2000, n_cores = 1){

  # Validate inputs
  if (!inherits(model, 'Seurat')) {
    stop("'model' must be a Seurat object returned by scPAS()")
  }
  if (is.null(.getMisc(model, "scPAS_para"))) {
    stop("The model does not contain scPAS parameters. Please run scPAS() first.")
  }
  
  model_params <- .getMisc(model, "scPAS_para")

  if(inherits(test.data, 'Seurat')){
    if(do_imputation){
      test.data <- imputation(test.data, assay = assay, method = imputation_method)
      assay <- Seurat::DefaultAssay(test.data)
    }
    test.exp <- .getAssayData(object = test.data, assay = assay, slot_name = 'data')
    Expression_cell <- test.exp
    rownames(Expression_cell) <- rownames(test.exp)
    colnames(Expression_cell) <- colnames(test.exp)
  }else{
    test.exp <- as.matrix(test.data)
    Expression_cell <- methods::as(preprocessCore::normalize.quantiles(as.matrix(test.exp)), 'dgCMatrix')
    rownames(Expression_cell) <- rownames(test.exp)
    colnames(Expression_cell) <- colnames(test.exp)
  }

  Coefs <- model_params$Coefs
  common <- intersect(names(Coefs), rownames(Expression_cell))

  # Check feature overlap
  n_nonzero_coefs <- sum(Coefs[common] != 0)
  if(n_nonzero_coefs < 20){
    warning(sprintf("Only %d non-zero features overlap between model and test data. Results may be unreliable.", n_nonzero_coefs))
    if(n_nonzero_coefs == 0){
      stop("No valid features found. The test data may not be suitable for this model.")
  }
  }
  
  message(sprintf("Using %d common features (%d with non-zero coefficients)", 
                  length(common), n_nonzero_coefs))

  Coefs <- Coefs[common]
  Expression_cell <- Expression_cell[common, ]
  
  # Scale expression data
  scaled_exp <- sparse_row_scale(Expression_cell, center = TRUE, scale = TRUE)
  colnames(scaled_exp) <- colnames(Expression_cell)
  rownames(scaled_exp) <- rownames(Expression_cell)
  scaled_exp[is.na(scaled_exp)] <- 0
  
  # Permutation test for significance
  if (n_cores > 1) {
    message(sprintf("Running permutation test with %d permutations using %d cores...", 
                    permutation_times, n_cores))
  } else {
    message(sprintf("Running permutation test with %d permutations (sequential)...", 
                    permutation_times))
  }
  
  perm_results <- parallel_permutation_test(
    scaled_exp = scaled_exp,
    Coefs = Coefs,
    permutation_times = permutation_times,
    n_cores = n_cores,
    seed = 12345
  )
  
  risk_score <- perm_results$risk_score
  risk_score.background <- perm_results$risk_score.background
  
  # Ensure risk_score.background is a matrix
  if (!is.matrix(risk_score.background)) {
    risk_score.background <- as.matrix(risk_score.background)
  }

  if(independent){
    mean.background <- Matrix::rowMeans(risk_score.background)
    sd.background <- apply(risk_score.background, 1, stats::sd)
    # Handle zero standard deviation
    sd.background[sd.background == 0 | !is.finite(sd.background)] <- 1
  }else{
    mean.background <- mean(as.matrix(risk_score.background))
    sd.background <- stats::sd(as.matrix(risk_score.background))
    if (sd.background == 0 || !is.finite(sd.background)) {
      sd.background <- 1
    }
  }

  Z <- (risk_score[,1] - mean.background) / sd.background

  # Two-tailed p-value for proper statistical testing
  p.value <- 2 * stats::pnorm(q = abs(Z), mean = 0, sd = 1, lower.tail = FALSE)
  p.value <- pmin(p.value, 1)  # Cap at 1
  q.value <- stats::p.adjust(p = p.value, method = 'BH')

  risk_score_data.frame <- data.frame(
    sample = colnames(Expression_cell),
    scPAS_RS = risk_score[,1],
    scPAS_NRS = Z,
    scPAS_Pvalue = p.value,
    scPAS_FDR = q.value,
    stringsAsFactors = FALSE
  )
  risk_score_data.frame$scPAS <- ifelse(Z > 0 & q.value <= FDR.threshold, 'scPAS+', 
                                        ifelse(Z < 0 & q.value <= FDR.threshold, 'scPAS-', '0'))

  if(inherits(test.data, 'Seurat')){
    test.data <- Seurat::AddMetaData(test.data, metadata = risk_score_data.frame$scPAS_RS, col.name = "scPAS_RS")
    test.data <- Seurat::AddMetaData(test.data, metadata = risk_score_data.frame$scPAS_NRS, col.name = "scPAS_NRS")
    test.data <- Seurat::AddMetaData(test.data, metadata = risk_score_data.frame$scPAS_Pvalue, col.name = "scPAS_Pvalue")
    test.data <- Seurat::AddMetaData(test.data, metadata = risk_score_data.frame$scPAS_FDR, col.name = "scPAS_FDR")
    test.data <- Seurat::AddMetaData(test.data, metadata = risk_score_data.frame$scPAS, col.name = "scPAS")
    
    # Report summary
    n_pos <- sum(risk_score_data.frame$scPAS == "scPAS+")
    n_neg <- sum(risk_score_data.frame$scPAS == "scPAS-")
    message(sprintf("Prediction complete: %d scPAS+ cells, %d scPAS- cells (FDR < %.2f)", 
                    n_pos, n_neg, FDR.threshold))
    return(test.data)
  }else{
    return(risk_score_data.frame)
  }
}
