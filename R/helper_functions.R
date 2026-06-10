#' Validate a gene signature against available genes
#'
#' Checks whether the genes in a user-supplied signature are present in a
#' reference set of genes with available pre-computed ranks. The function
#' returns the subset of signature genes that can be used for scoring and the
#' subset that are missing.
#'
#' A warning is issued if one or more signature genes are not found in
#' `all_genes`. An error is raised if none of the supplied signature genes are
#' present.
#'
#'
#' @param signature Character vector containing the genes that define the
#'   signature to be scored.
#'
#'
#' @param all_genes Character vector of genes for which pre-computed ranks are
#'   available.
#'
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{valid_genes}{Character vector of signature genes found in
#'   `all_genes` and therefore available for scoring.}
#'   \item{missing_genes}{Character vector of signature genes not found in
#'   `all_genes`.}
#' }
#'
#'
#' @details
#' The order of genes in `valid_genes` matches their order in the input
#' `signature`. Missing genes are reported in a warning message. If no
#' signature genes are found in `all_genes`, the function stops with an error.
#'
#'
#' @examples
#' validate_signature(
#'   signature = c("CCR7", "CD62L"),
#'   all_genes = c("CCR7", "CD27", "CD28")
#' )
#'
#' validate_signature(
#'   signature = gene_signature,
#'   all_genes = rownames(ranks)
#' )
#'
#'
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




#' Impute ranks for missing signature genes
#'
#' Replaces missing signature genes by assigning each of them a rank equal to
#' `rank_cap`. The imputed ranks are appended to the supplied vector of
#' available gene ranks, producing a complete rank vector for downstream
#' signature scoring.
#'
#'
#' @param incomplete_ranks Numeric vector containing the ranks of signature
#'   genes that are present in a cell, sample, or spatial location.
#'
#' @param rank_cap Numeric value used to impute ranks for missing genes.
#'   Typically represents the maximum rank considered during scoring.
#'
#'
#' @param missing_genes Character vector of signature genes that are absent
#'   from the input data and therefore require rank imputation.
#'
#'
#' @return A numeric vector containing the original ranks in
#'   `incomplete_ranks` followed by imputed ranks for each gene in
#'   `missing_genes`.
#'
#'
#' @details
#' One value equal to `rank_cap` is added for each gene listed in
#' `missing_genes`. If `missing_genes` is empty, the input vector is returned
#' unchanged.
#'
#'
#' @examples
#' impute_missing_gene_ranks(
#'   incomplete_ranks = c(1, 2, 3, 4, 5, 6, 7, 8),
#'   rank_cap = 5,
#'   missing_genes = c("CD27", "CD28")
#' )
#'
#'
impute_missing_gene_ranks <- function(incomplete_ranks,rank_cap, missing_genes){

    # Create the additional rows of imputed ranks as a separate matrix
    imputed_ranks <- rep(rank_cap, length(missing_genes))

    # Append the imputed matrix to in incomplete gene ranks
    full_gene_ranks <- c(incomplete_ranks, imputed_ranks)

    return(full_gene_ranks)
}


#' Append a numeric vector as a new row to a matrix-like object
#'
#' Appends a numeric vector to a matrix-like object while preserving the
#' original object class. Supported input types include base R matrices,
#' Matrix::sparseMatrix objects, and DelayedArray::DelayedMatrix objects.
#'
#' The supplied vector is converted to the appropriate one-row representation
#' before being combined with the input object.
#'
#'
#' @param matrix_like_object A matrix-like object to which a new row will be
#' appended. Supported classes are base matrices, sparse matrices from the
#' Matrix package, and delayed matrices from the DelayedArray package.
#' @param numeric_vector Numeric vector to append as a new row. Its length must
#' equal the number of columns in matrix_like_object.
#'
#'
#' @return An object of the same class as matrix_like_object, with
#' numeric_vector appended as an additional row.
#'
#'
#' @details
#' The function validates that the length of numeric_vector matches the
#' number of columns in matrix_like_object. An error is raised if the
#' dimensions are incompatible or if the supplied object class is not
#' supported.
#'
#'
#' @examples
#' # Base matrix
#' mat <- matrix(1:6, nrow = 2)
#' append_to_matrix_like_object(mat, c(7, 8, 9))
#'
#' # Sparse matrix
#' sparse_mat <- Matrix::Matrix(mat, sparse = TRUE)
#' append_to_matrix_like_object(sparse_mat, c(7, 8, 9))
#'
#' # Delayed matrix
#' delayed_mat <- DelayedArray::DelayedArray(mat)
#' append_to_matrix_like_object(delayed_mat, c(7, 8, 9))
#'
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


#' Compute column-wise geometric mean expression values
#'
#' Computes the geometric mean expression for each sample, cell, or spatial
#' location in a count matrix. These values are used to derive rank caps for
#' downstream signature scoring.
#'
#' The geometric mean is calculated as:
#' \deqn{\exp(\mathrm{mean}(\log(1 + x)))}
#' where \eqn{x} represents the expression values in a column.
#'
#'
#' @param counts A matrix-like object containing expression counts, with
#' genes in rows and samples, cells, or spatial locations in columns.
#' Supported inputs include base matrices, sparse matrices
#' (Matrix::sparseMatrix), and delayed matrices
#' (DelayedArray::DelayedMatrix).
#'
#'
#' @return A numeric vector containing the column-wise geometric mean
#' expression values. The length of the vector equals the number of columns
#' in counts.
#'
#'
#' @details
#' A pseudocount of 1 is added to all expression values before log
#' transformation to avoid undefined values for zero counts. Missing values are
#' ignored when computing the mean log expression.
#'
#'
#' @examples
#' counts <- matrix(
#' c(0, 5, 10,
#' 2, 0, 20,
#' 1, 3, 0),
#' nrow = 3
#' )
#'
#' compute_geometric_average(counts)
#'
#'
compute_geometric_average <-function(counts){
    cap_values <- apply(counts, 2,
                        function(x){exp(mean(log(1+ x),
                                             na.rm = TRUE))})
    return(cap_values)
}


#' Generate ranks from expression count data
#'
#' Computes column-wise gene expression ranks from a count matrix and derives
#' rank caps for each sample, cell, or spatial location. Ranking is performed
#' using MatrixGenerics::colRanks(), allowing efficient processing of dense,
#' sparse, and delayed matrix representations.
#'
#' Rank caps are incorporated by temporarily appending a row of cap expression
#' values to the count matrix before ranking. The resulting rank of each cap
#' value is returned separately from the gene rank matrix.
#'
#'
#' @param counts A matrix-like object containing expression counts, with
#' genes in rows and samples, cells, or spatial locations in columns.
#' Supported inputs include base matrices, sparse matrices
#' (Matrix::sparseMatrix), and delayed matrices
#' (DelayedArray::DelayedMatrix).
#'
#'
#' @param count_caps Optional numeric vector containing the expression values
#' used as rank caps for each column. The length must equal the number of
#' columns in counts. By default, column-wise geometric mean
#' expression values computed by compute_geometric_average() are used.
#'
#'
#' @param handle_ties Character string specifying how tied expression values
#' are ranked. Passed directly to MatrixGenerics::colRanks() via the
#' ties.method argument. Supported options are:
#' \describe{
#' \item{"min"}{Assign the minimum rank to tied values (default).}
#' \item{"max"}{Assign the maximum rank to tied values.}
#' \item{"average"}{Assign the average rank to tied values.}
#' \item{"random"}{Break ties at random.}
#' }
#'
#'
#' @return A named list with two elements:
#' \describe{
#' \item{ranks}{A matrix of column-wise gene ranks. Lower rank values
#' correspond to higher expression levels. The matrix is typically of type
#' integer, except when handle_ties = "average", in which case it is
#' numeric.}
#'
#' \item{rank_caps}{A numeric vector containing the rank assigned to the cap
#' value in each column.}
#' }
#'
#'
#' @details
#' Gene expression values are ranked in descending order by applying
#' MatrixGenerics::colRanks() to the negated count matrix. Consequently, the
#' most highly expressed gene in a column receives rank 1.
#'
#' To determine cap ranks, the supplied cap values are appended as an
#' additional row prior to ranking. After ranking, gene ranks and cap ranks are
#' separated and returned independently.
#'
#' Informative messages are printed indicating whether ranking is being
#' performed on a base matrix, sparse matrix, or delayed matrix object.
#'
#'
#' @examples
#' # Compute SPARO ranks using geometric mean expression as rank caps
#' rank_results <- get_ranks_from_counts(counts)
#'
#' # Specify custom cap values
#' rank_results <- get_ranks_from_counts(
#' counts,
#' count_caps = rep(10, ncol(counts))
#' )
#'
#'
get_ranks_from_counts <- function(counts,
                                count_caps =
                                    compute_geometric_average(counts),
                                handle_ties = "min"){



    # check validity of input
    if(!(handle_ties %in%
         c("min", "max", "average", "random"))){
        stop("SPAROscore says: Invalid value provided for handle_ties
             Limit to using 'min', 'max', 'average' or 'random'")
    }

    # Compute rank caps' expression values or use user input
    if(length(count_caps) == ncol(counts)){
        cap_values <- count_caps
    }
    else{
        stop("SPAROscore says: Length of user entered count_caps are
        not matching the column count od counts")
    }



    # get original rownames and column names fo count_data
    count_data_rnames <- dimnames(counts)[[1]]
    count_data_cnames <- dimnames(counts)[[2]]

    # append expression caps as the last row to count_data  before ranking
    counts <- append_to_matrix_like_object(counts, cap_values)


    # print message to user as to what object is being used
    if(inherits(counts, "DelayedMatrix")){
        message("SPAROscore says: Ranking a DelayedMatrix object")
    }
    else if(inherits(counts, "sparseMatrix")){
        message("SPAROscore says: Ranking a sparseMatrix object")
    }
    else if(is.matrix(counts)){
        message("SPAROscore says: Ranking a matrix object")
    }



    # Rank each column of gene expressions along with the expression caps
    full_rank_data <- MatrixGenerics::colRanks(-counts,
                                          ties.method = handle_ties,
                                          preserveShape = TRUE,
                                          useNames =  FALSE)

    # port column names and rownames from count_data to full_rank_data
    dimnames(full_rank_data)[[1]] <- c(count_data_rnames, "rank_caps")
    dimnames(full_rank_data)[[2]] <- count_data_cnames

    # separate the ranks and rank caps and return

    return(list(ranks = full_rank_data[count_data_rnames, ],
                rank_caps = full_rank_data["rank_caps", ]))
}






#' Compute a SPAROscore for a single sample, cell, or spatial location
#'
#' Calculates the SPAROscore for a single expression profile using the ranks of
#' genes in a signature and a corresponding rank cap. The score is based on the
#' normalized Spearman footrule distance between the observed signature gene
#' ranks and the worst-case ranking in which all genes are assigned the rank
#' cap.
#'
#' Prior to scoring, all ranks greater than rank_cap are truncated to
#' rank_cap. Missing signature genes can either be ignored or imputed with
#' the capped rank.
#'
#'
#' @param signature_ranks_vector Numeric or integer vector containing the ranks
#' of signature genes for a single sample, cell, or spatial location.
#'
#'
#' @param rank_cap Numeric value used to cap gene ranks before score
#' calculation. Typically corresponds to the rank of the geometric mean
#' expression value for the sample.
#'
#'
#' @param handle_missing_genes Character string specifying how signature genes
#' absent from the input data should be handled:
#' \itemize{
#' \item "skip" (default): exclude missing genes from the score
#' calculation.
#' \item "impute": include missing genes by assigning them the capped
#' rank value (rank_cap).
#' }
#'
#'
#' @param missing_geneset Character vector containing signature genes that are
#' not present in the input dataset.
#'
#'
#' @return A single numeric value representing the SPAROscore for the specified
#' sample, cell, or spatial location.
#'
#'
#' @details
#' The SPAROscore is computed in four steps:
#' \enumerate{
#' \item All ranks greater than rank_cap are replaced with rank_cap.
#' \item Missing genes are optionally imputed with the capped rank.
#' \item The Spearman footrule distance between the observed ranks and the
#' worst-case ranking (all genes assigned rank_cap) is calculated.
#' \item The distance is normalized by the maximum theoretically possible
#' Spearman footrule distance for a signature of the same size.
#' }
#'
#' Scores range from 0 to 1, where larger values indicate that signature genes
#' tend to have higher expression ranks relative to the cap.
#'
#'
#' @examples
#' # validate gene signatures
#' num_genes <- nrow(counts)
#' valid_gene_signature <- validate_signature(signatures,
#' rownames(counts))
#'
#' valid_genes <- valid_gene_signature$valid_genes
#' missing_genes <- valid_gene_signature$missing_genes
#'
#' # get sparoroanks
#' ranks <- get_ranks_from_counts(counts)
#'
#'
#'
#' compute_sparoscore_per_cell(
#' signature_ranks_vector = ranks$ranks[signatures, cell_id],
#' rank_cap = ranks$rank_caps[cell_id],
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



#' Compute SPAROscores for a gene signature
#'
#' Computes SPAROscores for all samples, cells, or spatial locations using
#' pre-computed SPARO ranks and rank caps. Scores quantify the enrichment of a
#' gene signature based on the normalized Spearman footrule distance between
#' observed signature gene ranks and the corresponding rank cap.
#'
#' Signature genes that are not present in the ranked dataset can either be
#' excluded from the calculation or imputed using the capped rank value.
#'
#'
#' @param ranks A matrix of gene ranks produced by
#' get_ranks_from_counts(). Rows correspond to genes and columns
#' correspond to samples, cells, or spatial locations.
#'
#'
#' @param rank_caps Numeric vector containing the rank cap for each column of
#' ranks. Typically obtained from the rank_caps element returned by
#' get_ranks_from_counts().
#'
#'
#' @param signatures Character vector containing the genes that define the
#' signature of interest. Gene identifiers must match the row names of
#' ranks.
#'
#'
#' @param handle_missing_genes Character string specifying how signature genes
#' absent from ranks should be handled:
#' \itemize{
#' \item "skip" (default): exclude missing genes from score
#' calculation.
#' \item "impute": assign missing genes the capped rank value and include
#' them in the score calculation.
#' }
#'
#'
#' @return A named numeric vector of SPAROscores, with one score per column of
#' ranks.
#'
#'
#' @details
#' The supplied signature is first validated against the genes present in
#' ranks. Missing genes are reported and handled according to
#' handle_missing_genes.
#'
#' For each column, gene ranks are extracted for the signature genes and passed
#' to compute_sparoscore_per_cell(), which computes a normalized Spearman
#' footrule distance score. Higher scores indicate that signature genes tend to
#' occupy higher expression ranks relative to the column-specific rank cap.
#'
#'
#' @examples
#' # Generate SPARO ranks and rank caps
#' rank_results <- get_ranks_from_counts(counts)
#'
#' # Compute scores for a gene signature
#' sparoscores <- compute_sparoscores(
#' ranks = rank_results$ranks,
#' rank_caps = rank_results$rank_caps,
#' signatures = c("CCR7", "IL7R", "LTB")
#' )
#'
#' # Include missing signature genes by imputing capped ranks
#' sparoscores <- compute_sparoscores(
#' ranks = rank_results$ranks,
#' rank_caps = rank_results$rank_caps,
#' signatures = c("CCR7", "IL7R", "LTB"),
#' handle_missing_genes = "impute"
#' )
#'
#' @seealso
#' \code{\link{get_ranks_from_counts}},
#' \code{\link{compute_sparoscore_per_cell}},
#' \code{\link{validate_signature}}
#'
#'
compute_sparoscores <- function(ranks,
                             rank_caps,
                             signatures,
                             handle_missing_genes = "skip"){


    # Get all available genes names from the ranks
    available_genes <- rownames(ranks)

    # validate the signatures
    valid_gene_signature <- validate_signature(signatures, available_genes)

    # subset the ranks for only the valid_genes from signatures
    signature_rank_matrix <- ranks[valid_gene_signature$valid_genes, ]

    # get the rank cap values
    if(length(rank_caps) == ncol(ranks)){
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




#' Add SPAROscores to a Seurat object
#'
#' Computes SPAROscores for a gene signature and stores the resulting scores in
#' the metadata of a Seurat object. SPARO ranks can either be calculated from
#' an expression layer or retrieved from previously stored SPARO ranks.
#'
#' Optionally, the computed SPARO ranks can be stored as an assay layer for
#' reuse in subsequent analyses.
#'
#' @param seurat_object A Seurat object containing expression data.
#'
#' @param signatures Character vector containing the genes that define the
#' signature to score. Gene identifiers must match the feature names in the
#' selected assay.
#'
#' @param assay Character string specifying the assay from which expression
#' data should be extracted. Defaults to "RNA".
#'
#' @param layer Character string specifying the assay layer containing
#' expression counts. Defaults to "counts".
#'
#' @param is_assay_ranks Logical indicating whether previously
#' calculated SPARO ranks stored in the assay layer "ranks" should be
#' used instead of recalculating ranks from expression data.
#'
#' @param store_ranks Logical indicating whether newly computed SPARO
#' ranks should be stored in the selected assay under the layer name
#' "ranks".
#'
#' @param count_caps Optional numeric vector of expression values used to
#' derive rank caps during SPARO rank calculation. If NULL, geometric mean
#' expression values are computed using compute_geometric_average().
#'
#' @param rank_caps Optional numeric vector of rank cap values used during
#' SPAROscore calculation. If NULL, rank caps returned by
#' get_ranks_from_counts() are used.
#'
#' @param handle_ties Character string specifying how tied expression values
#' are ranked. Passed to get_ranks_from_counts(). Supported values are
#' "min", "max", "average", and "random".
#'
#' @param handle_missing_genes Character string specifying how signature genes
#' absent from the dataset should be handled:
#' \itemize{
#' \item "skip" (default): exclude missing genes from score calculation.
#' \item "impute": include missing genes by assigning them the capped
#' rank value.
#' }
#'
#' @return The input Seurat object with:
#' \itemize{
#' \item A metadata column containing SPAROscores for the supplied gene
#' signature.
#' \item A metadata column named "SPARO_rank_caps" containing the rank cap
#' for each cell.
#' \item Optionally, a "ranks" layer in the selected assay if
#' store_ranks = TRUE.
#' }
#'
#' @details
#' If is_assay_ranks = FALSE, SPARO ranks are computed from the
#' selected assay layer before score calculation.
#'
#' If is_assay_ranks = TRUE, the function expects a layer named
#' "ranks" to already exist in the selected assay. Stored rank caps are
#' retrieved from the metadata column "SPARO_rank_caps" unless explicitly
#' provided through rank_caps.
#'
#' Computed SPAROscores are appended to the Seurat object's metadata and can be
#' accessed through seurat_object[[]].
#'
#' @examples
#' # Compute SPAROscores directly from RNA counts
#' seurat_object <- augment_sparoscores_seurat(
#' seurat_object = seurat_object,
#' signatures = c("CCR7", "IL7R", "LTB")
#' )
#'
#' # Store SPARO ranks for later reuse
#' seurat_object <- augment_sparoscores_seurat(
#' seurat_object = seurat_object,
#' signatures = c("CCR7", "IL7R", "LTB"),
#' store_ranks = TRUE
#' )
#'
#' # Reuse previously stored SPARO ranks
#' seurat_object <- augment_sparoscores_seurat(
#' seurat_object = seurat_object,
#' signatures = c("CCR7", "IL7R", "LTB"),
#' is_assay_ranks = TRUE
#' )
#'
#' @seealso
#' \code{\link{get_ranks_from_counts}},
#' \code{\link{compute_sparoscores}}
augment_sparoscores_seurat <- function(seurat_object,
                                       signatures,
                                       assay = "RNA",
                                       layer = "count",
                                       is_assay_ranks = FALSE,
                                       store_ranks = FALSE,
                                       count_caps = NULL,
                                       rank_caps = NULL,
                                       handle_ties = "min",
                                       handle_missing_genes = "skip"){


    if(!is.logical(is_assay_ranks)){
        stop("SPAROscore says: is_assay_ranks should be a boolean")
    }

    if(!is.logical(store_ranks)){
        stop("SPAROscore says: store_ranks should be a boolean")
    }

    # if user ask to use pre-calculated ranks
    if(is_assay_ranks){
        # extract precalculated ranks and rank_caps
        ranks <- Seurat::GetAssayData(seurat_object,
                                   assay = assay,
                                   layer = "ranks")

        # set rank_caps to pre-calculated values unless provided by the user
        if(is.null(rank_caps)){
            rank_caps <- as.numeric(unlist(seurat_object[["SPARO_rank_caps"]]))
            names(rank_caps) <- rownames(seurat_object[[]])
        }
        else{
            rank_caps <- rank_caps
        }


        sparoscores <- get_scores(ranks =  ranks,
                                       rank_caps = rank_caps,
                                       signatures =  signatures,
                                       handle_missing_genes =
                                           handle_missing_genes)
    }
    else{
        #extract count matrix from data
        counts <- Seurat::GetAssayData(seurat_object,
                                        assay = assay,
                                        layer = layer)

        # set count caps to geometric averages unless provided by the user
        if(is.null(count_caps)){
            count_caps <- compute_geometric_average(counts)
        }
        else{
            count_caps <- count_caps
        }


        # get outputs for sparoranking function
        get_ranks_output <- get_ranks(counts = counts,
                                                count_caps = count_caps,
                                                handle_ties = handle_ties)

        # set ranks for future input
        ranks <- get_ranks_output$ranks


        # set rank_caps to geometric averages unless provided by the user
        if(is.null(rank_caps)){
            rank_caps <- get_ranks_output$rank_caps
        }
        else{
            rank_caps <- rank_caps
        }


        sparoscores <- get_scores(ranks =  ranks,
                                       rank_caps = rank_caps,
                                       signatures =  signatures,
                                       handle_missing_genes =
                                           handle_missing_genes)
    }



    # store the ranks in ranks layer of the assay if user asks
    if(store_ranks){
        seurat_object <- Seurat::SetAssayData(object = seurat_object,
                                              assay = assay,
                                              layer = "ranks",
                                              new.data = ranks)
    }


    # append the rank_caps and sparoscores to metadata
    sparoscores_df <- as.data.frame(sparoscores,
                                    row.names = rownames(sparoscores))

    sparoscores_df$SPARO_rank_caps <- rank_caps

    seurat_object <- Seurat::AddMetaData(object = seurat_object,
                                         metadata = sparoscores_df)

    return(seurat_object)
}


#' Add SPAROscores to a SummarizedExperiment, SingleCellExperiment,
#' SpatialExperiment or RangedSummarizedExperiment object
#'
#' Computes SPAROscores for a gene signature and stores the resulting scores
#' in the `colData` of a SummarizedExperiment, SingleCellExperiment,
#' SpatialExperiment or RangedSummarizedExperiment object
#'
#' SPARO ranks can either be computed from an expression assay or retrieved
#' from previously stored assay data.
#'
#' Optionally, computed SPARO ranks can be stored in the assay for reuse.
#'
#' @param sce_object A SummarizedExperiment, SingleCellExperiment,
#' SpatialExperiment or RangedSummarizedExperiment object
#'
#'
#' @param signatures Character vector of genes defining the signature.
#' Gene identifiers must match `rownames()` of the selected assay.
#'
#'
#' @param assay Character string specifying the assay to use. Defaults to "RNA".
#'
#'
#' @param is_assay_ranks Logical indicating whether previously stored
#' ranks in an assay named `"ranks"` should be used.
#'
#'
#' @param store_ranks Logical indicating whether to store computed
#' ranks in the assay `"ranks"`.
#'
#'
#' @param count_caps Optional numeric vector of expression values used to
#' derive rank caps during SPARO rank calculation. If NULL, geometric mean
#' expression values are computed using compute_geometric_average().
#'
#'
#' @param rank_caps Optional numeric vector of rank cap values used during
#' SPAROscore calculation. If NULL, rank caps returned by
#' get_ranks_from_counts() are used.
#'
#'
#' @param handle_ties Character string specifying how tied expression values
#' are ranked. Passed to get_ranks_from_counts(). Supported values are
#' "min", "max", "average", and "random".
#'
#'
#' @param handle_missing_genes Character string specifying how signature genes
#' absent from the dataset should be handled:
#' \itemize{
#' \item "skip" (default): exclude missing genes from score calculation.
#' \item "impute": include missing genes by assigning them the capped
#' rank value.
#' }
#'
#'
#' @return The input object with:
#' \itemize{
#'   \item SPAROscores added to `colData(sce_object)`
#'   \item A column `"SPARO_rank_caps"` in `colData`
#'   \item Optionally, a `"ranks"` assay if `store_ranks = TRUE`
#' }
#'
#'
#' @details
#' If `is_assay_ranks = FALSE`, ranks are computed from the selected
#' assay. Otherwise, stored `"ranks"` assay values are used.
#'
#' SPAROscores are computed using the normalized Spearman footrule distance
#' between signature gene ranks and rank caps.
#'
#'
#' @examples
#' library(SingleCellExperiment)
#'
#' sce <- augment_sparoscores_sce(
#'   sce_object = sce,
#'   signatures = c("CCR7", "IL7R", "LTB")
#' )
#'
#' @seealso
#' \code{\link{get_ranks_from_counts}},
#' \code{\link{get_scores}},
#' \code{\link{compute_geometric_average}}
#'
#'
augment_sparoscores_sce <- function(sce_object,
                            signatures,
                            assay = "RNA",
                            is_assay_ranks = FALSE,
                            store_ranks = FALSE,
                            count_caps = NULL,
                            rank_caps = NULL,
                            handle_ties = "min",
                            handle_missing_genes = "skip"){

    if(!is.logical(is_assay_ranks)){
        stop("SPAROscore says: is_assay_ranks should be a boolean")
    }

    if(!is.logical(store_ranks)){
        stop("SPAROscore says: store_ranks should be a boolean")
    }

    # if user ask to use pre-calculated ranks
    if(is_assay_ranks){
        # extract precalculated ranks and rank_caps
        ranks <- SummarizedExperiment::assay(sce_object,
                                           assay = "ranks")

        # set rank_caps to pre-calculated values unless provided by the user
        if(is.null(rank_caps)){
            rank_caps <- as.numeric(unlist(colData(sce_object)[["SPARO_rank_caps"]]))
            names(rank_caps) <- rownames(colData(sce_object))
        }
        else{
            rank_caps <- rank_caps
        }


        sparoscores <- get_scores(ranks =  ranks,
                                       rank_caps = rank_caps,
                                       signatures =  signatures,
                                       handle_missing_genes =
                                           handle_missing_genes)
    }
    else{
        #extract count matrix from data
        counts <- SummarizedExperiment::assay(sce_object,
                                                   assay = assay)

        # set count caps to geometric averages unless provided by the user
        if(is.null(count_caps)){
            count_caps <- compute_geometric_average(counts)
        }
        else{
            count_caps <- count_caps
        }


        # get outputs for sparoranking function
        get_ranks_output <- get_ranks(counts = counts,
                                                count_caps = count_caps,
                                                handle_ties = handle_ties)

        # set ranks for future input
        ranks <- get_ranks_output$ranks


        # set rank_caps to geometric averages unless provided by the user
        if(is.null(rank_caps)){
            rank_caps <- get_ranks_output$rank_caps
        }
        else{
            rank_caps <- rank_caps
        }


        sparoscores <- get_scores(ranks =  ranks,
                                       rank_caps = rank_caps,
                                       signatures =  signatures,
                                       handle_missing_genes =
                                           handle_missing_genes)
    }



    # store the ranks in ranks layer of the assay if user asks
    if(store_ranks){
        SummarizedExperiment::assay(
            sce_object, "ranks") <- ranks
    }


    # append the rank_caps and sparoscores to metadata
    sparoscores_df <- as.data.frame(sparoscores,
                                    row.names = rownames(sparoscores))

    sparoscores_df$SPARO_rank_caps <- rank_caps

    colData(sce_object) <- cbind(colData(sce_object), sparoscores_df)

    return(sce_object)
}



