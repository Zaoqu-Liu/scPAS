# scPAS: Single-Cell Phenotype-Associated Subpopulation Identifier

<!-- badges: start -->
[![R-universe](https://zaoqu-liu.r-universe.dev/badges/scPAS)](https://zaoqu-liu.r-universe.dev/scPAS)
[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen.svg)](https://github.com/Zaoqu-Liu/scPAS)
[![R Version](https://img.shields.io/badge/R-%E2%89%A5%204.0.0-blue.svg)](https://www.r-project.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://img.shields.io/badge/DOI-10.1093%2Fbib%2Fbbae655-blue.svg)](https://doi.org/10.1093/bib/bbae655)
[![pkgdown](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://zaoqu-liu.github.io/scPAS/)
<!-- badges: end -->

## Overview

**scPAS** is a computational framework for identifying phenotype-associated cell subpopulations from single-cell RNA sequencing (scRNA-seq) data through the integration of bulk transcriptomic profiles and clinical phenotypes. The methodology employs network-regularized sparse regression to quantify the strength of association between individual cells and phenotypic outcomes, enabling both quantitative scoring and statistical inference at single-cell resolution.

<p align="center">
  <img src="man/figures/Flow_diagram.png" alt="scPAS Workflow" width="85%"/>
</p>

## Methodological Framework

### Statistical Model

scPAS implements the Augmented and Penalized Minimization with L0 (APML0) algorithm, which optimizes the following objective function:

$$\hat{\boldsymbol{\beta}} = \arg\min_{\boldsymbol{\beta}} \left\{ L(\boldsymbol{\beta}; \mathbf{X}, \mathbf{y}) + \lambda_1 \|\boldsymbol{\beta}\|_1 + \lambda_2 \boldsymbol{\beta}^T \mathbf{L} \boldsymbol{\beta} \right\}$$

where:
- $L(\boldsymbol{\beta})$ denotes the loss function (Gaussian, binomial, or Cox partial likelihood)
- $\lambda_1 \|\boldsymbol{\beta}\|_1$ enforces sparsity through L1 regularization
- $\lambda_2 \boldsymbol{\beta}^T \mathbf{L} \boldsymbol{\beta}$ incorporates network structure via Laplacian regularization

### Supported Regression Families

| Family | Phenotype Type | Application |
|--------|----------------|-------------|
| **Gaussian** | Continuous | Age, BMI, gene expression levels |
| **Binomial** | Binary | Case-control, treatment response |
| **Cox** | Time-to-event | Overall survival, progression-free survival |

### Statistical Inference

Significance is assessed through permutation testing with false discovery rate (FDR) correction using the Benjamini-Hochberg procedure. Cells are classified as:

- **scPAS+**: Positively associated (risk score > 0, FDR < threshold)
- **scPAS−**: Negatively associated (risk score < 0, FDR < threshold)  
- **Non-significant**: FDR ≥ threshold

## Installation

### From R-universe (Recommended)

```r
install.packages("scPAS", repos = c(
  "https://zaoqu-liu.r-universe.dev",
  "https://cloud.r-project.org"
))
```

### From GitHub

```r
# Install devtools if necessary
if (!require("devtools")) install.packages("devtools")

# Install scPAS
devtools::install_github("Zaoqu-Liu/scPAS")
```

### Dependencies

```r
# Bioconductor dependency
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("preprocessCore")

# Optional: parallel computing
install.packages(c("future", "future.apply"))
```

## Documentation

Comprehensive documentation is available at the [pkgdown website](https://zaoqu-liu.github.io/scPAS/).

| Vignette | Description |
|----------|-------------|
| [Quick Start](https://zaoqu-liu.github.io/scPAS/articles/quick-start.html) | Basic usage and workflow |
| [Algorithm](https://zaoqu-liu.github.io/scPAS/articles/algorithm.html) | Mathematical methodology |
| [Visualization](https://zaoqu-liu.github.io/scPAS/articles/visualization.html) | Publication-quality figures |
| [Survival Analysis](https://zaoqu-liu.github.io/scPAS/articles/case-survival.html) | Cox regression application |
| [Binary Classification](https://zaoqu-liu.github.io/scPAS/articles/case-binary.html) | Treatment response prediction |

## Usage

### Basic Example

```r
library(scPAS)
library(Seurat)

# Execute scPAS analysis
result <- scPAS(
  bulk_dataset = bulk_expression,    # Matrix: genes × samples
  sc_dataset = seurat_object,        # Seurat object
  phenotype = phenotype_vector,      # Phenotypic data

  family = "gaussian",               # Regression family
  nfeature = 3000,                   # Variable features
  permutation_times = 1000,          # Permutation iterations
  n_cores = 4                        # Parallel cores
)

# Extract significant cells
significant_cells <- subset(result, subset = scPAS_FDR < 0.05)
table(significant_cells$scPAS)
```

### Survival Analysis (Cox Regression)

```r
library(survival)

# Define survival phenotype
surv_phenotype <- Surv(time = clinical_data$time, 
                       event = clinical_data$status)

# Run Cox regression analysis
result <- scPAS(
  bulk_dataset = bulk_expression,
  sc_dataset = seurat_object,
  phenotype = surv_phenotype,
  family = "cox"
)
```

## Output Structure

scPAS appends the following columns to the Seurat object metadata:

| Column | Description |
|--------|-------------|
| `scPAS_RS` | Raw risk score |
| `scPAS_NRS` | Normalized risk score (Z-statistic) |
| `scPAS_Pvalue` | Permutation-based p-value |
| `scPAS_FDR` | Benjamini-Hochberg adjusted p-value |
| `scPAS` | Cell classification (scPAS+/scPAS−/0) |

## Key Features

- **Multi-modal Integration**: Bridges bulk and single-cell transcriptomics
- **Network Regularization**: Incorporates gene-gene co-expression structure
- **Flexible Phenotypes**: Supports continuous, binary, and survival outcomes
- **Scalable Computation**: Parallel processing for large-scale datasets
- **Seurat Integration**: Native support for Seurat v4 objects

## Citation

If you use scPAS in your research, please cite:

> Xie A, Wang H, Zhao J, Wang Z, Xu J, Xu Y. **scPAS: single-cell phenotype-associated subpopulation identifier.** *Briefings in Bioinformatics*. 2024;26(1):bbae655. doi: [10.1093/bib/bbae655](https://doi.org/10.1093/bib/bbae655)

```bibtex
@article{xie2024scpas,
  title={scPAS: single-cell phenotype-associated subpopulation identifier},
  author={Xie, Aimin and Wang, Hao and Zhao, Jianqiang and Wang, Zhe and Xu, Jing and Xu, Yang},
  journal={Briefings in Bioinformatics},
  volume={26},
  number={1},
  pages={bbae655},
  year={2024},
  publisher={Oxford University Press},
  doi={10.1093/bib/bbae655}
}
```

## Authors

**Original Author**
- Aimin Xie (aiminyy1993@gmail.com)

**Maintainer**
- [Zaoqu Liu](https://orcid.org/0000-0002-0452-742X) (liuzaoqu@163.com)
  - Department of Interventional Radiology, The First Affiliated Hospital of Zhengzhou University

## License

This package is distributed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html).

## Acknowledgments

- Original algorithm development by Aimin Xie et al.
- Package maintenance and optimization by Zaoqu Liu
- Computational infrastructure support from Zhengzhou University
