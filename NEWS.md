# scPAS 1.0.4

## API Compatibility

* Full compatibility with both Seurat v4 (4.0.0+) and Seurat v5 (5.0.0+)
* Added internal helper functions for version-aware API calls:
  - `.getAssayData()`: Auto-detects `slot` vs `layer` parameter
  - `.getGraphs()`: Compatible graph data access
  - `.getMisc()` / `.setMisc()`: Compatible misc slot access
  - `.createAssayObject()`: Compatible assay creation
* Updated all vignettes to use Seurat-version-aware code patterns

---

# scPAS 1.0.3

## Bug Fixes

* Fixed p-value calculation to use proper two-tailed test (multiplied by 2)
* Improved numerical stability in `sparse.cor()` function with handling of zero-variance columns
* Fixed division-by-zero protection in C++ code for zero standard deviation columns
* Fixed parameter name conflict (`imputation` → `do_imputation`) to avoid collision with function name
* Fixed `tag` parameter handling when NULL in binomial family
* Fixed loop variable scope issues in main scPAS function

## Improvements

* Added parallel computing support via `n_cores` parameter
* Enhanced input validation with informative error messages
* Improved sparse matrix operations in `sparse_row_scale()`
* Cleaned up code comments and naming conventions
* Added comprehensive test suite (57 tests)

## Documentation

* New vignettes:
  - Quick Start Guide (`quick-start.Rmd`)
  - Algorithm and Methodology (`algorithm.Rmd`)
  - Visualization Gallery (`visualization.Rmd`)
  - Case Study: Cancer Survival Analysis (`case-survival.Rmd`)
  - Case Study: Binary Classification (`case-binary.Rmd`)
* Added pkgdown website support
* Updated README with detailed usage examples

---

# scPAS 1.0.2

## Bug Fixes

* Fixed sparse matrix transpose issues in `sparse.cor()` function
* Fixed `rowMeans()` and `colMeans()` handling for sparse matrices
* Fixed `FindNeighbors()` rownames requirement
* Improved logical indexing to avoid `which()` errors with sparse matrices
* Enhanced NA handling in correlation calculations

## Improvements

* Added parallel computing support for permutation tests
* Improved documentation and examples

---
  
# scPAS 1.0.1

## Bug Fixes

* Fixed cross-platform compatibility issues (Windows, macOS, Linux)
* Fixed DLL/SO loading on different platforms

---

# scPAS 1.0.0

* Initial release
* Network-regularized sparse regression for phenotype-associated cell identification
* Support for Gaussian, binomial, and Cox regression families
* Integration with Seurat v4
* KNN and ALRA imputation methods
