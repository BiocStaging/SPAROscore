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
        stop("SPAROscore says:
             No signature genes are present in the input dataset")
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
    imputed_ranks <- rep(rank_cap, length(missing_genes))

    # Append the imputed matrix to in incomplete gene ranks
    full_gene_ranks <- c(incomplete_ranks, imputed_ranks)

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
        stop("SPAROscore says:
             Cannot append cap values. Length not equal to columns")
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
        stop("SPAROscore says: Unsupported file type for counts/caps")
    }
}


#' Compute geometric average of counts
#'
#'compute_geometric_average() computes geomwtric average expression for
#' each sample/cell/spot to be used to find rank caps
#'
#' @param counts_data A matrixilike object of counts where
#' columns are spots/cells/samples and rows are genes.
#' Can handle matrix, sparsematrix, delayedmatrix
#'
#' @returns A numeric vector of column-wise geometric averages
#'
#' @examples
#' compute_geometric_average(counts_data)
#'
#'
compute_geometric_average <-function(counts_data){
    cap_values <- apply(counts_data, 2,
                        function(x){exp(mean(log(1+ x),
                                             na.rm = TRUE))})
    return(cap_values)
}


#' Compute column-wise ranks with the respective rank caps
#'
#' get_sparoranks_from_counts() takes counts matrix data as input
#' and returns column wise ascending ranks matrix.
#' Uses MatrixGenerics::colRanks()
#' Ties are handles by minimum by default
#'
#' @param counts_data A matrix-like object of counts where
#' columns are spots/cells/samples and rows are genes.
#' Can handle matrix, sparsematrix, delayedmatrix
#'
#'
#' @param count_caps A optional numeric vector comprising the expresison values
#' for each column to be used as capping value. By default the geometric average
#' counts of each column will be used as count_caps
#'
#'
#' @param handle_ties A character string specifying how ties are treated.
#' "min" by default. Can be "min", "max", "average", or "random".
#'
#'
#' @returns A matrix of type integer. Has one additional row with the rank caps
#' If handle_ties = "average" then it is a matrix of type numeric.
#'
#'
#' @examples
#' get_sparoranks_from_counts(counts_data,
#' handle_ties = "min", rank_cap_metric = "geometric_average")
get_sparoranks_from_counts <- function(counts_data,
                                count_caps =
                                    compute_geometric_average(counts_data),
                                handle_ties = "min"){



    # check validity of input
    if(!(handle_ties %in%
         c("min", "max", "average", "random"))){
        stop("SPAROscore says: Invalid value provided for handle_ties
             Limit to using 'min', 'max', 'average' or 'random'")
    }

    # Compute rank caps' expression values or use user input
    if(length(count_caps) == ncol(counts_data)){
        cap_values <- count_caps
    }
    else{
        stop("SPAROscore says: Length of user entered count_caps are
        not matching the column count od counts_data")
    }



    # get original rownames and column names fo count_data
    count_data_rnames <- dimnames(counts_data)[[1]]
    count_data_cnames <- dimnames(counts_data)[[2]]

    # append expression caps as the last row to count_data  before ranking
    counts_data <- append_to_matrix_like_object(counts_data, cap_values)


    # print message to user as to what object is being used
    if(inherits(counts_data, "DelayedMatrix")){
        message("SPAROscore says: Ranking a DelayedMatrix object")
    }
    else if(inherits(counts_data, "sparseMatrix")){
        message("SPAROscore says: Ranking a sparseMatrix object")
    }
    else if(is.matrix(counts_data)){
        message("SPAROscore says: Ranking a matrix object")
    }



    # Rank each column of gene expressions along with the expression caps
    full_rank_data <- MatrixGenerics::colRanks(-counts_data,
                                          ties.method = handle_ties,
                                          preserveShape = TRUE,
                                          useNames =  FALSE)

    # port column names and rownames from count_data to full_rank_data
    dimnames(full_rank_data)[[1]] <- c(count_data_rnames, "rank_caps")
    dimnames(full_rank_data)[[2]] <- count_data_cnames

    # separate the ranks and rank caps and return

    return(list(sparoranks = full_rank_data[count_data_rnames, ],
                rank_caps = full_rank_data["rank_caps", ]))
}






#' Compute SPAROscore for a single sample/cell/spot for a single signature
#'
#' compute_sparoscore_per_cell() is a function to compute SPAROscore for a
#' single column based on the gene rankings of signature genes and rank cap
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
#' compute_sparoscore_per_cell(
#' signature_ranks_vector = sparoranks[signature_genes, cell_id][1:num_genes],
#' rank_cap = sparoranks[signature_genes, cell_id][num_genes + 1],
#' handle_missing_genes = "skip",
#' missing_geneset = missing_genes)
#'
#'
#'
compute_sparoscore_per_cell <- function(signature_ranks_vector, rank_cap,
                               handle_missing_genes = "skip", missing_geneset){


    # check validity of input
    if(!(handle_missing_genes %in% c("skip", "impute"))){
        stop("SPAROscore says: Invalid value provided for handle_missing_genes.
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



#' Compute SPAROscore: sparsity aware robust gene signature scores
#'
#'compute_sparoscores() computes signature score for each column of the data.
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
#'
#' @examples
#' sparorank_output <- get_sparoranks(counts_data)
#' sparoscores <- compute_sparoscores(sparoranks = sparorank_output$sparoranks,
#'                                 rank_caps = sparorank_output$rank_caps
#'                                 signature_genes = c("geneA", "geneC"))
#'
#'
compute_sparoscores <- function(sparoranks,
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


    sparoscores <- vapply(
        X = names(sparoscores),
        FUN = function(cell) {
            compute_sparoscore_per_cell(
                signature_ranks_vector = signature_rank_matrix[, cell],
                rank_cap = rank_caps[cell],
                handle_missing_genes = handle_missing_genes,
                missing_geneset = valid_gene_signature$missing_genes
            )
        },
        FUN.VALUE = numeric(1)
    )
    return(sparoscores)
}




#' Title
#'
#' @param seurat_object
#' @param signature_genes
#' @param assay
#' @param layer
#' @param use_existing_sparoranks
#' @param store_sparoranks
#' @param count_caps
#' @param rank_caps
#' @param handle_ties
#' @param handle_missing_genes
#'
#' @returns
#' @export
#'
#' @examples
augment_sparoscores_seurat <- function(seurat_object,
                                       signature_genes,
                                       assay = NULL,
                                       layer = NULL,
                                       use_existing_sparoranks = FALSE,
                                       store_sparoranks = FALSE,
                                       count_caps = NULL,
                                       rank_caps = NULL,
                                       handle_ties = "min",
                                       handle_missing_genes = "skip"){

    # assign the default assay as RNA unless provided
    assay <- ifelse(is.null(assay), "RNA", assay)

    # assign the default layer as counts unless provided
    layer <- ifelse(is.null(layer), "count", layer)

    # if user ask to use pre-calculated ranks
    if(use_existing_sparoranks){
        # extract precalculated sparoranks and rank_caps
        sparoranks <- Seurat::GetAssayData(seurat_object,
                                   assay = assay,
                                   layer = "sparoranks")

        # set rank_caps to pre-calculated values unless provided by the user
        if(is.null(rank_caps)){
            rank_caps <- as.numeric(unlist(seurat_object[["SPARO_rank_caps"]]))
            names(rank_caps) <- rownames(seurat_object[[]])
        }
        else{
            rank_caps <- rank_caps
        }


        sparoscores <- get_sparoscores(sparoranks =  sparoranks,
                                       rank_caps = rank_caps,
                                       signature_genes =  signature_genes,
                                       handle_missing_genes =
                                           handle_missing_genes)
    }
    else{
        #extract count matrix from data
        counts_data <- Seurat::GetAssayData(seurat_object,
                                        assay = assay,
                                        layer = layer)

        # set count caps to geometric averages unless provided by the user
        if(is.null(count_caps)){
            count_caps <- compute_geometric_average(counts_data)
        }
        else{
            count_caps <- count_caps
        }


        # get outputs for sparoranking function
        get_sparoranks_output <- get_sparoranks(counts_data = counts_data,
                                                count_caps = count_caps,
                                                handle_ties = handle_ties)

        # set sparoranks for future input
        sparoranks <- get_sparoranks_output$sparoranks


        # set rank_caps to geometric averages unless provided by the user
        if(is.null(rank_caps)){
            rank_caps <- get_sparoranks_output$rank_caps
        }
        else{
            rank_caps <- rank_caps
        }


        sparoscores <- get_sparoscores(sparoranks =  sparoranks,
                                       rank_caps = rank_caps,
                                       signature_genes =  signature_genes,
                                       handle_missing_genes =
                                           handle_missing_genes)
    }



    # store the ranks in sparoranks layer of the assay if user asks
    if(store_sparoranks){
        seurat_object <- Seurat::SetAssayData(object = seurat_object,
                                              assay = assay,
                                              layer = "sparoranks",
                                              new.data = sparoranks)
    }


    # append the rank_caps and sparoscores to metadata
    sparoscores_df <- as.data.frame(sparoscores,
                                    row.names = rownames(sparoscores))

    sparoscores_df$SPARO_rank_caps <- rank_caps

    seurat_object <- Seurat::AddMetaData(object = seurat_object,
                                         metadata = sparoscores_df)

    return(seurat_object)
}
