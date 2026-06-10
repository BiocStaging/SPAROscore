
# set the generic for get_sparoscores() methods
# Should work for seurat, sce, spe, sme, matrix, spmatrix, delayedarrays, dataframe
#' Title
#'
#' @param data
#' @param signatures
#' @param assay
#' @param layer
#' @param is_assay_ranks
#' @param store_ranks
#' @param count_caps
#' @param rank_caps
#' @param handle_ties
#' @param handle_missing_genes
#'
#' @returns
#' @export
#'
#' @examples
setGeneric("sparoscore",
           function(data,
                    signatures,
                    assay = "RNA",
                    layer = "count",
                    is_assay_ranks = FALSE,
                    store_ranks = FALSE,
                    count_caps = NULL,
                    rank_caps = NULL,
                    handle_ties = "min",
                    handle_missing_genes = "skip")
               standardGeneric("sparoscore"))


# set the method for sparoscore() where data is a
# Seurat object
setMethod("sparoscore",
          signature('Seurat','ANY','ANY','ANY','ANY',
                    'ANY','ANY','ANY','ANY', "ANY"),
          function(data,
                   signatures,
                   assay = "RNA",
                   layer = "count",
                   is_assay_ranks = FALSE,
                   store_ranks = FALSE,
                   count_caps = NULL,
                   rank_caps = NULL,
                   handle_ties = "min",
                   handle_missing_genes = "skip"){

              #call the helper function
              seurat_object <- augment_sparoscores_seurat(
                  seurat_object = data,
                  signatures = signatures,
                  assay = assay,
                  layer = layer,
                  is_assay_ranks = is_assay_ranks,
                  store_ranks = store_ranks,
                  count_caps = count_caps,
                  rank_caps = rank_caps,
                  handle_ties = handle_ties,
                  handle_missing_genes = handle_missing_genes)

              return(seurat_object)
          }
)



# set the method for sparoscore() where data is a
# SummarizedExperiment object
setMethod("sparoscore",
          signature('SummarizedExperiment','ANY','ANY','ANY','ANY',
                    'ANY','ANY','ANY','ANY', "ANY"),
          function(data,
                   signatures,
                   assay = "RNA",
                   layer = "count",
                   is_assay_ranks = FALSE,
                   store_ranks = FALSE,
                   count_caps = NULL,
                   rank_caps = NULL,
                   handle_ties = "min",
                   handle_missing_genes = "skip"){

              #call the helper function
              sce_object <- augment_sparoscores_seurat(
                  seurat_object = data,
                  signatures = signatures,
                  assay = assay,
                  is_assay_ranks = is_assay_ranks,
                  store_ranks = store_ranks,
                  count_caps = count_caps,
                  rank_caps = rank_caps,
                  handle_ties = handle_ties,
                  handle_missing_genes = handle_missing_genes)

              return(sce_object)
          }
)
