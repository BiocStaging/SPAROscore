calc_sparoscore <- function(rank_matrix, gene_signature, cap_metric = "empirical_cutoff_gene_expression", handle_missing_genes = "impute"){
    # function to score samples/cells based on spearman footrule distance metric
    # Parameters
    # rank_matrix: Matrix of original gene ranks of all cells/samples for all genes. Calculated with counts_to_ranks()
    # gene_signature: A vector of gene names. Naming scheme should matching the rownames of rank_matrix
    # cap_metric: A character to mention which cap metric to use. Cap metrics are calculated in counts_to_ranks(). C
    #             can be "mean_gene_expression", "empirical_cutoff_gene_expression", "mean_nz_gene_expression", "median_nz_gene_expression", "ga_nz_gene_expression", "ga_1_gene_expression"
    # handle_missing_genes: The way to handle missing genes from the gene_signature. "impute" gives zero expression to the genes. "skip" removes those genes from analyses

    # Get all genes from the rank_matrix
    all_genes <- rownames(rank_matrix)

    # validate the gene_signature
    valid_gene_signature <- validate_geneset(gene_signature, all_genes)

    # subset the ranks for only the valid_genes from gene_signature
    signature_rank_matrix <- rank_matrix[valid_gene_signature$valid_genes, ]

    # calculate the maximum ranks (of zero expressing gene) for each cell
    rank_caps <- rank_matrix[cap_metric, ]

    # Calculate the metrics for the new algorithm
    new_algo_score <- setNames(rep(0.0, ncol(signature_rank_matrix)), colnames(signature_rank_matrix))

    for (cell in names(new_algo_score)) {
        new_algo_score[cell] <- calc_normalised_spearmans_footrule_distance(signature_rank_matrix[,cell], rank_caps[cell],
                                                                            handle_missing_genes = handle_missing_genes, missing_geneset = valid_gene_signature$missing_genes)
    }

    return(new_algo_score)
}
