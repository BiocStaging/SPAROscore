# SPAROscore

**SParsity-Adaptive RObust Gene Signature Scoring** for Transcriptomics Data

SPAROscore is a robust gene signature scoring method for bulk, single-cell, and spatial transcriptomics data. It is designed to perform consistently across diverse gene expression datasets by adapting to varying levels of data sparsity, enabling efficient and reliable signature scoring across multiple transcriptomic modalities.

SPAROscore quantifies the relative expression of a gene signature with respect to the background expression of each sample, cell, or spatial domain. This method makes the resulting scores easy to interpret biologically and facilitates comparisons across heterogeneous datasets.

## Installation

``` r
# Install from GitHub
devtools::install_github("MangiolaLaboratory/SPAROscore")
```

## Quick Start

``` r
library(SPAROscore)

# Basic usage with a matrix
signature <- c("gene1", "gene2", "gene3")
scores <- sparoscore(data = expression_matrix, signatures = signature)

# With multiple signatures
signatures <- list(
  pathway1 = c("geneA", "geneB", "geneC"),
  pathway2 = c("geneD", "geneE", "geneF")
)
scores <- sparoscore(data = expression_matrix, signatures = signatures)

# With Seurat objects
library(Seurat)
seurat_obj <- sparoscore(
  data = seurat_obj,
  signatures = signature,
  assay = "RNA",
  layer = "count"
)

# With Bioconductor objects (SingleCellExperiment, SpatialExperiment, etc.)
library(SingleCellExperiment)
sce <- sparoscore(data = sce, signatures = signature, assay = "counts")
```

## Documentation

## Supported Data Types

-   **Matrices**: `matrix`, `sparseMatrix`, `DelayedMatrix`, `data.frame`
-   **Seurat**: Full integration with Seurat objects
-   **Bioconductor**: Full integration with `SummarizedExperiment` and its derived objects like `SingleCellExperiment`, `SpatialExperiment`
-   **Signatures**: Character vectors, named lists, `GeneSet`, `GeneSetCollection`

## License

GPL-3.0 - See LICENSE file

\<!-- \## Citation

```         
Kamaraj, V. (2024). SPAROscore: Sparsity-Adaptive Robust Gene Signature Scoring. 
R package version 0.1.0.
```

## Contact

**Maintainer**: Venkatesh Kamaraj\
**Email**: [venkatesh.kamaraj\@adelaide.edu.au](mailto:venkatesh.kamaraj@adelaide.edu.au){.email}\
**Affiliation**: University of Adelaide -- !\>
