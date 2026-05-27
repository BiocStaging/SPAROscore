# library(matrixStats)



# counts_to_ranks <- function(counts_data, ties.method = "min"){
#     # Function to take counts data as input and return column wise ascending ranks matrix. Ties are handles by taking average
#     # Parameters:
#     # counts_data: A matrix of counts where columns are cells/samples and rows are genes
#     # ties.method: Average ranks are taken if counts are tied
#
#     empirical_cutoff <- function(x){
#         sparsity <- sum(x == 0)/length(x)
#         cutoff <- (1 + sparsity)*(mean(x, na.rm = TRUE))
#         return(cutoff)
#     }
#
#     # get the mean of gene expression of each cell
#     mean_gene_expression <- apply(counts_data, 2, function(x){mean(x, na.rm = TRUE)})
#
#     # get the empirical cut-off gene expression of each cell
#     empirical_cutoff_gene_expression <- apply(counts_data, 2, empirical_cutoff)
#
#     # get the mean of non-zero gene expression of each cell
#     mean_nz_gene_expression <- apply(counts_data, 2, function(x){mean(x[x > 0], na.rm = TRUE)})
#
#     # get the median of non-zero gene expression of each cell
#     median_nz_gene_expression <- apply(counts_data, 2, function(x){median(x[x > 0], na.rm = TRUE)})
#
#     # get the geometric mean of non-zero gene expression of each cell
#     ga_nz_gene_expression <- apply(counts_data, 2, function(x){exp(mean(log(x[x > 0])))})
#
#     # get the geometric mean of 1 + gene_expression of each cell
#     ga_1_gene_expression <- apply(counts_data, 2, function(x){exp(mean(log(1+ x)))})
#
#
#     # append mean gene expression of cells as the last row before ranking
#     counts_data <- rbind(counts_data, mean_gene_expression, empirical_cutoff_gene_expression, mean_nz_gene_expression, median_nz_gene_expression, ga_nz_gene_expression, ga_1_gene_expression)
#
#     # Rank the gene expressions for each cell along with the average expression value
#     # rank_data <- apply(counts_data, 2, function(x) rank(-x, ties.method = ties.method))
#     rank_data <- matrixStats::colRanks(-counts_data, ties.method = ties.method, preserveShape = TRUE)
#
#     rownames(rank_data) <- rownames(counts_data)
#     colnames(rank_data) <- colnames(counts_data)
#
#     return(rank_data)
# }


# validate_geneset <- function(geneset, all_genes){
#     # Function to identify and report the genes from the geneset, that the files are available for
#     # Parameters:
#     # geneset: A character vector: set of genes the user may be interested in scoring
#     # all_genes: A character vector: All the genes that the ranks are available for
#     # Returns:
#     # List of valid_genes and missing_genes
#     # valid_genes: a list of genes from the geneset that the pre-calculated ranks are available for
#     # missing_genes: a list of genes from the geneset that the pre-calculated ranks are not available for
#
#     valid_geneset <- geneset[geneset %in% all_genes]
#
#     invalid_genes <- geneset[!(geneset %in% all_genes)]
#
#     print(paste0("Warning: The pre-calculate ranks are not available for the following ", length(invalid_genes), " genes"))
#     print(invalid_genes)
#     print("These genes expression are imputed with mean gene expression of the sample/cell by default (or they can be skipped from the analyses)")
#
#     return(list(valid_genes = valid_geneset, missing_genes = invalid_genes))
# }


# impute_missing_gene_ranks <- function(incomplete_gene_ranks, rank_cap, missing_geneset){
#     # Function to impute the ranks of missing genes with rank_cap. As expression is imputed as 0, ranks by default become rank_cap
#     # incomplete_gene_ranks: VECTOR of original gene ranks of signature genes in a particular cell/sample
#     # rank_cap: NUMERIC of rank cap calculated for the sample/cell
#     # missing_geneset: the list of genes from gene_signature not present in the sample/cell
#
#     # Create the additional rows of imputed ranks as a separate matrix
#     imputed_ranks <- rep(rank_cap, length(missing_geneset))
#
#     # Append the imputed matrix to in incomplete gene ranks
#     full_gene_ranks <- c(incomplete_gene_ranks, imputed_ranks)
#
#     return(full_gene_ranks)
# }


# calc_normalised_spearmans_footrule_distance <- function(signature_ranks_vector, rank_cap, handle_missing_genes = "impute", missing_geneset){
#     # function to calculate spearman footrule distance for a gene ranking of a cell/sample
#     # Note: Distance is normalized with theoretically calculated maximum possible spearman foortrule distance.
#     # Note: The ranks are capped at r_cap before the calculating distance
#     # Parameters
#     # signature_ranks_vector: VECTOR of original gene ranks of the signature genes of a particular cell/sample
#     # rank_cap: Rank of the gene with at least nth_max level of expression for the cell/sample
#     # handle_missing_genes: The way to handle missing genes from the gene_signature. "impute" gives zero expression to the genes. "skip" removes those genes from analyses
#     # missing_geneset: the list of genes from gene_signature not present in the sample/cell
#
#
#     # Cap ranks at rank_cap, if they are more than rank_cap
#     signature_ranks_vector[signature_ranks_vector > rank_cap] <- rank_cap
#
#     # Handle missing genes in gene_signature
#     if(handle_missing_genes == "impute"){
#         if(length(missing_geneset) > 0){
#             signature_ranks_vector <- impute_missing_gene_ranks(signature_ranks_vector, rank_cap, missing_geneset)
#         }
#     }
#
#     # worst case scenario of ranks for the geneset is all of them being rank_cap
#     # Calculate Spearman's footrule distance from the worst case ranks
#     spearman_footrule <- sum(abs(rank_cap - signature_ranks_vector))
#
#     # print(paste0("numerator: ", spearman_footrule))
#
#     # calculate the total length of the gene set
#     n <- length(signature_ranks_vector)
#
#     # Calculate theoretically possible maximum possible Spearman footrule distance
#     max_spearman_footrule <- sum(abs(rank_cap - seq_along(signature_ranks_vector)))
#
#     # print(paste0("denominator: ", max_spearman_footrule))
#
#     # Calculate normalised spearman footrule distance
#     normalised_spearmans_footrule <- spearman_footrule/max_spearman_footrule
#
#     return(normalised_spearmans_footrule)
# }


calc_new_algorithm <- function(rank_matrix, gene_signature, cap_metric = "empirical_cutoff_gene_expression", handle_missing_genes = "impute"){
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


