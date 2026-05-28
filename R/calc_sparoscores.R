#' Compute SPAROscore: sparsity aware robust gene signature scores
#'
#'calc_sparoscores() computes signature score for each column of the input data.
#'Scores are based on Spearman Footrule Distance of the signature gene ranks
#'from the rank cap (default to the rank of geometric average expression)
#'for every sample/cell/spot
#'
#'
#' @param sparoranks matrix outputted by get_sparoranks()
#'
#'
#' @param rank_caps numeric vector containing column-wise rank caps to be used
#'
#'
#' @param signature_genes character vector of gene names/ids. Note that the
#' naming scheme should match the original input data.
#' (i.e) rownames(counts_data)
#'
#'
#' @param handle_missing_genes character string specifying how signature genes
#' missing in the imput dataset are handled.
#' Defaults to "skip", removing those genes from the analyses
#' "impute" adds zero expression values to all these genes in the input dataset
#'
#'
#' @returns A numeric with the scores for each sample/cell/spot
#'
#'
#' @export
#'
#' @examples
#' sparorank_output <- get_sparoranks(counts_data)
#' sparoscores <- calc_sparoscores(sparoranks = sparorank_output$sparoranks,
#'                                 rank_caps = sparorank_output$rank_caps
#'                                 signature_genes = c("geneA", "geneC"))
#'
#'
calc_sparoscores <- function(sparoranks,
                             rank_caps,
                             signature_genes,
                             handle_missing_genes = "skip"){


    # Get all available genes names from the sparoranks
    available_genes <- rownames(sparoranks)

    # validate the signature_genes
    valid_gene_signature <- validate_signature(signature_genes, available_genes)

    # subset the ranks for only the valid_genes from signature_genes
    signature_rank_matrix <- sparoranks[valid_gene_signature$valid_genes, ]

    # get the rank cap values
    if(length(rank_caps) == ncol(sparoranks)){
            rank_caps <- rank_caps
    }
    else{
        stop("SPAROscore says: Invalid number of ranks provided")
    }


    # Calculate sparoscores for each column
    sparoscores <- setNames(rep(0.0,
                            ncol(signature_rank_matrix)),
                            colnames(signature_rank_matrix))


    for (cell in names(sparoscores)) {
        sparoscores[cell] <- compute_sparoscore(
            signature_ranks_vector = signature_rank_matrix[,cell],
            rank_cap = rank_caps[cell],
            handle_missing_genes = handle_missing_genes,
            missing_geneset = valid_gene_signature$missing_genes)
    }

    return(sparoscores)
}
