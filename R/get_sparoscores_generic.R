#'Compute SPAROscore: sparsity aware robust gene signature scores
#'
#'get_sparoscores() computes signature score for each column of the data.
#'
#'For optimal use case, requires the output from get_sparoranks()
#'
#'Scores are based on Spearman Footrule Distance of the signature gene ranks
#'from the rank cap (default to the rank of geometric average expression)
#'for every sample/cell/spot
#'
#' @param sparoranks matrix named `sparoranks` outputted by get_sparoranks()
#'
#' @param rank_caps numeric vector named `rank_caps` outputted by
#' get_sparoranks() for optimal results.
#' Each of values will be used as the rank caps for the respective columns
#' In case of a custom input, the names should match that of sparoranks
#'
#' @param signature_genes for a single gene signature,
#' simple character vector of signature gene names/ids or
#' a GeneSet S4 object
#'
#' for multiple gene signatures
#' a named list of names/ids where the names denote the signature themselves or
#' a GeneSetCollection S4 object
#'
#' Note: naming scheme of the genes should match that of sparoranks
#'
#'
#' @param handle_missing_genes character string specifying how signature genes
#' missing in the imput dataset are handled.
#' Defaults to "skip", removing those genes from the analyses
#' "impute" adds zero expression values to all these genes
#'
#'
#' @returns A numeric vector for single gene signatures and a matrix for
#' multiple signatures
#'
#'
#' @export
#'
#'
#'
#' @examples
#' sparorank_output <- get_sparoranks(counts_data)
#'
#' sparoscores <- get_sparoscores(sparoranks = sparorank_output$sparoranks,
#'                                 rank_caps = sparorank_output$rank_caps
#'                                 signature_genes = c("geneA", "geneC"))
#'
#'
#'sparoscores <- get_sparoscores(sparoranks = sparorank_output$sparoranks,
#'                                 rank_caps = sparorank_output$rank_caps
#'                                 signature_genes = list(
#'                                  "sigA" = c("geneA", "geneC"),
#'                                  "sigB" = c("geneB", "geneD")))
#'


# set the generic for get_sparoscores() methods
setGeneric("get_sparoscores",
           function(sparoranks,
                    rank_caps,
                    signature_genes,
                    handle_missing_genes = "skip")
               standardGeneric("get_sparoscores"))




# set the method for get_sparoscores() where signature_genes is a
# character vector
setMethod("get_sparoscores",
          signature('ANY','ANY', 'character', "ANY"),
          function(sparoranks,
                   rank_caps,
                   signature_genes,
                   handle_missing_genes = "skip"){

              #call the helper function
              sparoscores <- compute_sparoscores(
                  sparoranks,
                  rank_caps,
                  signature_genes,
                  handle_missing_genes = "skip")
              return(sparoscores)
          }
)


# set the method for get_sparoscores() where signature_genes is a
# named list
setMethod("get_sparoscores",
          signature('ANY','ANY', 'list', "ANY"),
          function(sparoranks,
                   rank_caps,
                   signature_genes,
                   handle_missing_genes = "skip"){

              #call the helper function with vapply
              sparoscores <- vapply(X = signature_genes,
                                   FUN = function(x){
                                              compute_sparoscores(
                                              sparoranks,
                                              rank_caps,
                                              x,
                                              handle_missing_genes = "skip")
                                          },
                                   FUN.VALUE = numeric(ncol(sparoranks)))
              return(sparoscores)
          }
)



# set the method for get_sparoscores() where signature_genes is a
# GeneSet S4 class
setMethod("get_sparoscores",
          signature('ANY','ANY', 'GeneSet', "ANY"),
          function(sparoranks,
                   rank_caps,
                   signature_genes,
                   handle_missing_genes = "skip"){

              #call the helper function
              sparoscores <- compute_sparoscores(
                  sparoranks,
                  rank_caps,
                  GSEABase::geneIds(signature_genes),
                  handle_missing_genes = "skip")
              colnames(sparoscores) <- names(signature_genes)
              return(sparoscores)
          }
)


# set the method for get_sparoscores() where signature_genes is a
# GeneSetCollection S4 class
setMethod("get_sparoscores",
          signature('ANY','ANY', 'GeneSetCollection', "ANY"),
          function(sparoranks,
                   rank_caps,
                   signature_genes,
                   handle_missing_genes = "skip"){

              #call the helper function with vapp
              sparoscores <- vapply(X = GSEABase::geneIds(signature_genes),
                                    FUN = function(x){
                                        compute_sparoscores(
                                            sparoranks,
                                            rank_caps,
                                            x,
                                            handle_missing_genes = "skip")
                                    },
                                    FUN.VALUE = numeric(ncol(sparoranks)))
              colnames(sparoscores) <- names(signature_genes)
              return(sparoscores)
          }
)
