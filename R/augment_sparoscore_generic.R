
# set the generic for get_sparoscores() methods
setGeneric("augment_sparoscores",
           function(RNA_object,
                    signature_genes,
                    assay = "RNA",
                    layer = "count",
                    use_existing_sparoranks = FALSE,
                    store_sparoranks = FALSE,
                    count_caps = NULL,
                    rank_caps = NULL,
                    handle_ties = "min",
                    handle_missing_genes = "skip")
               standardGeneric("augment_sparoscores"))


# set the method for augment_sparoscores() where RNA_object is a
# Seurat object
setMethod("augment_sparoscores",
          signature('Seurat','ANY','ANY','ANY','ANY',
                    'ANY','ANY','ANY','ANY', "ANY"),
          function(RNA_object,
                   signature_genes,
                   assay = "RNA",
                   layer = "count",
                   use_existing_sparoranks = FALSE,
                   store_sparoranks = FALSE,
                   count_caps = NULL,
                   rank_caps = NULL,
                   handle_ties = "min",
                   handle_missing_genes = "skip"){

              #call the helper function
              seurat_object <- augment_sparoscores_seurat(
                  seurat_object = RNA_object,
                  signature_genes = signature_genes,
                  assay = assay,
                  layer = layer,
                  use_existing_sparoranks = use_existing_sparoranks,
                  store_sparoranks = store_sparoranks,
                  count_caps = count_caps,
                  rank_caps = rank_caps,
                  handle_ties = handle_ties,
                  handle_missing_genes = handle_missing_genes)

              return(seurat_object)
          }
)



# set the method for augment_sparoscores() where RNA_object is a
# SummarizedExperiment object
setMethod("augment_sparoscores",
          signature('SummarizedExperiment','ANY','ANY','ANY','ANY',
                    'ANY','ANY','ANY','ANY', "ANY"),
          function(RNA_object,
                   signature_genes,
                   assay = "RNA",
                   layer = "count",
                   use_existing_sparoranks = FALSE,
                   store_sparoranks = FALSE,
                   count_caps = NULL,
                   rank_caps = NULL,
                   handle_ties = "min",
                   handle_missing_genes = "skip"){

              #call the helper function
              sce_object <- augment_sparoscores_seurat(
                  seurat_object = RNA_object,
                  signature_genes = signature_genes,
                  assay = assay,
                  use_existing_sparoranks = use_existing_sparoranks,
                  store_sparoranks = store_sparoranks,
                  count_caps = count_caps,
                  rank_caps = rank_caps,
                  handle_ties = handle_ties,
                  handle_missing_genes = handle_missing_genes)

              return(sce_object)
          }
)
