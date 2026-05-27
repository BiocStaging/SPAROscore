#' validate_geneset() is a function to identify/report the genes from the
#' user defined signature gene set that can/can't be used in scoring
#'
#' @param geneset A character vector of genes that the user wants to score for
#' @param all_genes A character vector of all the genes with ranks available for
#'
#'
#' @returns A list of two character vectors. Genes for which ranks are available
#' followed by genes for which ranks are not available
#' @export
#'
#' @examples validate_geneset(gene_signature, rownames(sparo_ranks))
validate_geneset <- function(geneset, all_genes){
    # Function to
    # Parameters:
    # geneset:
    # all_genes:
    # Returns:
    # List of valid_genes and missing_genes
    # valid_genes: a list of genes from the geneset that the pre-calculated ranks are available for
    # missing_genes: a list of genes from the geneset that the pre-calculated ranks are not available for

    valid_geneset <- geneset[geneset %in% all_genes]

    invalid_genes <- geneset[!(geneset %in% all_genes)]

    print(paste0("Warning: The pre-calculate ranks are not available for the following ", length(invalid_genes), " genes"))
    print(invalid_genes)
    print("These genes expression are imputed with mean gene expression of the sample/cell by default (or they can be skipped from the analyses)")

    return(list(valid_genes = valid_geneset, missing_genes = invalid_genes))
}
