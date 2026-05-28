#' Identify valid and missing signature genes
#'
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
#' validate_signature(gene_signature, rownames(sparoranks))
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




#' Impute missing genes ranks with cap values
#'
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


#' Append numeric vector to matrix like objects. Geeneralised rbind()
#'
#' append_to_matrix_like_object() function converts a numeric vector to
#' matrix, or sparseMatrix, or DelayedArrays and then appends it to
#' the respective matrix object. This way, original matrix type is conserved.
#'
#'
#' @param matrix_like_object a matrix, or sparseMatrix, or DelayedArrays
#' @param numeric_vector numeric vector
#'
#' @returns A matrix-like object with an additional row
#'
#' @examples
#' append_to_matrix_like_object(counts_data, cap_values)
#'
append_to_matrix_like_object <- function(matrix_like_object, numeric_vector){
    if(ncol(matrix_like_object) != length(numeric_vector)){
        stop("Cannot append cap values. Length not equal to columns")
    }

    if(inherits(matrix_like_object, "DelayedMatrix")){
        #convert numeric vector to delayedmatrix
        delayed_vector <- DelayedArray::DelayedArray(
            matrix(numeric_vector, nrow = 1))
        return(DelayedArray::rbind(matrix_like_object, delayed_vector))
    }
    else if(inherits(matrix_like_object, "sparseMatrix")){
        #convert numeric vector to sparsematrix
        sparse_vector <- Matrix::Matrix(numeric_vector, nrow = 1, sparse = TRUE)
        return(rbind(matrix_like_object, sparse_vector))
    }
    else if(is.matrix(matrix_like_object)){
        return(rbind(matrix_like_object, numeric_vector))
    }
    else{
        stop("Unsupported file type for counts/caps")
    }
}


#' Compute column-wise ranks with the respective rank caps
#'
#' get_sparoranks_from_counts() takes counts matrix data as input
#' and returns column wise ascending ranks matrix.
#' Uses MatrixGenerics::colRanks()
#' Ties are handles by minimum by default
#'
#' @param counts_data A matrixilike object of counts where
#' columns are spots/cells/samples and rows are genes.
#' Can handle matrix, sparsematrix, delayedmatrix
#'
#' @param handle_ties A character string specifying how ties are treated.
#' "min" by default. Can be "min", "max", "average", or "random".
#'
#' @param rank_cap_metric A character string specifying the type of rank cap.
#' Default: "geometric_average" cap is the geometric average of (1 + count)
#' "manual_override" cap is user mentioned constant expression values.
#' "manual_override" requires manual_expr_caps to be NOT NULL
#'
#' @param manual_expr_caps used only if rank_cap_metric = "manual_override"
#' A numeric vector of length equals total columns in count_data
#'
#' @returns A matrix of type integer. Has one additional row with the rank caps
#' If handle_ties = "average" then it is a matrix of type numeric.
#'
#' @examples
#' get_sparoranks_from_counts(counts_data,
#' handle_ties = "min", rank_cap_metric = "geometric_average")
get_sparoranks_from_counts <- function(counts_data,
                                handle_ties = "min",
                                rank_cap_metric = "geometric_average",
                                manual_expr_caps = NULL){


    # check validity of input
    if(!(rank_cap_metric %in%
         c("geometric_average", "manual_override"))){
        stop("Invalid value provided for rank_cap_metric
             Limit to using 'geometric_average' or 'manual_override'")
    }

    # check validity of input
    if(!(handle_ties %in%
         c("min", "max", "average", "random"))){
        stop("Invalid value provided for handle_ties
             Limit to using 'min', 'max', 'average' or 'random'")
    }

    # Compute rank caps' expression values as per user input
    if(rank_cap_metric == "geometric_average"){
        # Geometric mean of (1 + count) of each cell as the expression cap
        cap_values <- apply(counts_data, 2,
                                 function(x){exp(mean(log(1+ x),
                                                      na.rm = TRUE))})
    }

    # else if(rank_cap_metric == "arithmetic_average"){
    #     # Arithmetic mean of (count) of each cell as the expression cap
    #     cap_values <- apply(counts_data, 2,
    #                              function(x){mean(x, na.rm = TRUE)})
    # }

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
            cap_values <- manual_expr_caps
        }

    }



    # append expression caps as the last row to count_data  before ranking
    counts_data <- append_to_matrix_like_object(counts_data, cap_values)

    # Rank each column of gene expressions along with the expression caps
    rank_data <- MatrixGenerics::colRanks(-counts_data,
                                          ties.method = handle_ties,
                                          preserveShape = TRUE,
                                          useNames =  FALSE)

    # port column names and rownames from count_data to rank_data
    rownames(rank_data) <- rownames(counts_data)
    colnames(rank_data) <- colnames(counts_data)

    return(rank_data)
}






#' Compute SPAROscore for a single sample/cell/spot for a single signature
#'
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
#' # validate gene signatures
#' num_genes <- nrow(counts_data)
#' valid_gene_signature <- validate_signature(signature_genes,
#' rownames(counts_data))
#'
#' valid_genes <- valid_gene_signature$valid_genes
#' missing_genes <- valid_gene_signature$missing_genes
#'
#' # get sparoroanks
#' sparoranks <- get_sparoranks_from_counts(counts_data)
#'
#'
#'
#' compute_sparoscore(
#' signature_ranks_vector = sparoranks[signature_genes, cell_id][1:num_genes],
#' rank_cap = sparoranks[signature_genes, cell_id][num_genes + 1],
#' handle_missing_genes = "skip",
#' missing_geneset = missing_genes)
#'
#'
#'
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
    final_sparoscore <- spearman_footrule/max_spearman_footrule

    return(final_sparoscore)
}
