#' validate_signature() is a function to identify and report the genes from the
#' user defined signature gene set for which pre-calculated are ranks available
#' and can be used in signature scoring
#'
#' @param signature A character vector of genes for which the user wants to
#' compute signature score
#'
#' @param all_genes A character vector of all the genes for which
#' pre-calculated ranks are available
#'
#'
#' @returns A list of two character vectors. Genes for which ranks are available
#' followed by genes for which ranks are not available
#'
#' @examples
#' validate_signature(c("CCR7", "CD62L"), c("CCR7", "CD27",  "CD28"))
#'
#' validate_signature(gene_signature, rownames(sparo_ranks))
validate_signature <- function(signature, all_genes){

    valid_signature <- signature[signature %in% all_genes]

    invalid_genes <- signature[!(signature %in% all_genes)]

    if(length(invalid_genes) > 0){
        warning(paste0("The following ", length(invalid_genes),
                       " signature genes are missing in the input dataset: ",
                       paste0(invalid_genes, collapse = ", ")))
    }

    if(length(valid_signature) == 0){
        stop("No signature genes are present in the input dataset")
    }

    return(list(valid_genes = valid_signature, missing_genes = invalid_genes))
}




#' impute_missing_gene_ranks() is a function to impute the ranks of
#' missing genes with the rank_cap.
#'
#' @param incomplete_ranks A numeric vector of the gene ranks  in a cell/sample
#'
#' @param rank_cap A numeric value of the rank_cap for the particular cell/spot
#'
#' @param missing_genes A character vector of the missing genes for which ranks
#' should be imputed
#'
#' @returns A numeric vector of the gene ranks with the missing genes
#'
#' @examples
#' impute_missing_gene_ranks(c(1,2,3,4,5,6,7,8), 5, c("CD27", "CD28"))
impute_missing_gene_ranks <- function(incomplete_ranks,rank_cap, missing_genes){

    # Create the additional rows of imputed ranks as a separate matrix
    imputed_ranks <- rep(rank_cap, length(missing_geneset))

    # Append the imputed matrix to in incomplete gene ranks
    full_gene_ranks <- c(incomplete_gene_ranks, imputed_ranks)

    return(full_gene_ranks)
}





#' get_counts_to_ranks() is a function that take counts matrix data as input
#' and returns column wise ascending ranks matrix. Uses matrixStats::colRanks()
#' Ties are handles by minimum by default
#'
#' @param counts_data A matrixilike object of counts where
#' columns are spots/cells/samples and rows are genes.
#' Can handle matrix, sparsematrix, delayedmatrix
#'
#' @param ties.method A character string specifying how ties are treated.
#' "min" by default. Can be "min", "max", "average", or "random".
#'
#' @param rank_cap_metric A character string specifying the type of rank cap.
#' Default: "geometric_average" cap is the geometric average of (1 + count)
#' "arithmetic_average" cap is the arithmetic average of counts
#' "manual_override" cap is a user mention constant expression value.
#' "manual_override" requires manual_expr_caps to be NOT NULL
#'
#' @param manual_expr_caps used only if rank_cap_metric = "manual_override"
#' A numeric vector of length equals total columns in count_data
#'
#' @returns A matrix of type integer. Has one additional row with the rank caps
#' If ties.method = "average" then it is a matrix of type numeric.
#'
#' @examples
#' get_counts_to_ranks
get_counts_to_ranks <- function(counts_data,
                                ties.method = "min",
                                rank_cap_metric = "geometric_average",
                                manual_expr_caps = NULL){


    # check validity of input
    if(!(rank_cap_metric %in%
         c("geometric_average", "arithmetic_average", "manual_override"))){
        stop("Invalid value provided for rank_cap_metric
             Limit to using 'geometric_average',
             'arithmetic_average', or 'manual_override'")
    }

    # Compute rank caps' expression values as per user input
    if(rank_cap_metric == "geometric_average"){
        # Geometric mean of (1 + count) of each cell as the expression cap
        expression_caps <- apply(counts_data, 2,
                                 function(x){exp(mean(log(1+ x),
                                                      na.rm = TRUE))})
    }

    else if(rank_cap_metric == "arithmetic_average"){
        # Arithmetic mean of (count) of each cell as the expression cap
        expression_caps <- apply(counts_data, 2,
                                 function(x){mean(x, na.rm = TRUE)})
    }

    else if(rank_cap_metric == "manual_override"){
        # A user defined expression value as the cap
        if(is.null(manual_expr_caps)){
            stop("manual_expr_caps missing.
                 It is required as rank_cap_metric is set to manual_override")
        }
        if(length(manual_expr_caps) != ncol(counts_data)){
            stop("User provided manual expression caps do not match the total
                 columns in the expression matrix")
        }
        if(length(manual_expr_caps) == ncol(counts_data)){
            expression_caps <- manual_expr_caps
        }

    }


    # append expression caps as the last row to count_data  before ranking
    counts_data <- rbind(counts_data, expression_caps)

    # Rank each column of gene expressions along with the expression caps
    rank_data <- MatrixGenerics::colRanks(-counts_data,
                                          ties.method = ties.method,
                                          preserveShape = TRUE)

    # port column names and rownames from count_data to rank_data
    rownames(rank_data) <- rownames(counts_data)
    colnames(rank_data) <- colnames(counts_data)

    return(rank_data)
}






#' compute_sparoscore() is a function to compute SPAROscore for a single column
#' based on the gene rankings of signature genes and rank cap
#' SPAROscore is Spearman Footrule Distance metric of signature ranks
#' from the rank cap. They are then normalized with the
#' theoretically maximum possible Spearman Footrule Distance.
#'
#' @param signature_ranks_vector integer vector of original gene ranks of the
#' signature genes for a single spot/cell/sample
#'
#' @param rank_cap numeric at which signature ranks are capped before computing
#' sparoscores. Defaults to rank of geometric average expression of cell/sample
#'
#' @param handle_missing_genes character string specifying how signature genes
#' missing in the imput dataset are handled.
#' Defaults to "skip", removing those genes from the analyses
#' "impute" adds zero expression values to all these genes in the input dataset
#'
#' @param missing_geneset character vector of signature genes missing from input
#'
#' @returns A numeric of sparoscore for the particular sample/cell/spot
#'
#' @examples
compute_sparoscore <- function(signature_ranks_vector, rank_cap,
                               handle_missing_genes = "skip", missing_geneset){


    # check validity of input
    if(!(handle_missing_genes %in% c("skip", "impute"))){
        stop("Invalid value provided for handle_missing_genes.
             Limit to using 'skip' or 'impute'")
    }


    # Cap ranks at rank_cap, if they are more than rank_cap
    signature_ranks_vector[signature_ranks_vector > rank_cap] <- rank_cap



    # Handle missing genes in gene_signature
    if(handle_missing_genes == "impute"){
        if(length(missing_geneset) > 0){
            signature_ranks_vector <- impute_missing_gene_ranks(
                signature_ranks_vector, rank_cap, missing_geneset)
        }
    }

    # worst case scenario of ranks for the geneset is all of them being rank_cap
    # Calculate Spearman's footrule distance from the worst case ranks
    spearman_footrule <- sum(abs(rank_cap - signature_ranks_vector))



    # Calculate theoretically possible max possible Spearman footrule distance
    max_spearman_footrule <- sum(abs(
        rank_cap - seq_along(signature_ranks_vector)))


    # Calculate normalised spearman footrule distance
    sparo_score <- spearman_footrule/max_spearman_footrule

    return(sparo_score)
}
