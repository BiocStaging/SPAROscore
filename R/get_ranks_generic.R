#' Compute gene ranks and rank caps for SPAROscore
#'
#' Generates gene rank matrices and rank caps required for SPAROscore
#' calculations.
#'
#' get_ranks() computes column-wise gene ranks from expression count
#' data and derives a rank cap for each sample, cell, or spatial location.
#' The resulting rank matrix and rank caps can be supplied directly to
#' \code{\link{compute_sparoscores}}.
#'
#' Internally, ranks are computed using
#' MatrixGenerics::colRanks(). Gene expression values are ranked in
#' descending order, such that the most highly expressed gene in a column
#' receives rank 1.
#'
#' This function is implemented as an S4 generic and supports multiple input
#' formats for counts including base matrices, sparse matrices,
#' delayed matrices, and data frames.
#'
#'
#' @param counts A matrix-like object containing gene expression counts,
#' with genes in rows and samples, cells, or spatial locations in columns.
#'
#' Supported input classes include:
#' \itemize{
#' \item matrix
#' \item Matrix::sparseMatrix
#' \item DelayedArray::DelayedMatrix
#' \item data.frame
#' }
#'
#'
#' @param count_caps Optional numeric vector containing count cap values for
#' each column. These values are used to determine the corresponding rank
#' caps. If NULL, column-wise geometric mean expression values computed by
#' \code{\link{compute_geometric_average}} are used.
#'
#'
#' @param handle_ties Character string specifying how tied expression values
#' should be ranked. Passed directly to
#' MatrixGenerics::colRanks(ties.method = ...).
#'
#' Supported values are:
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
#' \item{ranks}{
#' Matrix of gene ranks with the same dimensions as counts.
#' Ranks are integer-valued unless handle_ties = "average".
#' }
#'
#' \item{rank_caps}{
#' Numeric vector containing the rank cap associated with each column.
#' }
#' }
#'
#'
#' @details
#' A count cap value is appended to each column before ranking and its resulting
#' rank is extracted as the column-specific rank cap. By default, count caps
#' are derived from the geometric mean expression value of each column.
#'
#' The returned ranks matrix and rank_caps vector are intended for use
#' with \code{\link{get_scores}}.
#'
#' @export
#'
#' @examples
#' # Compute SPARO ranks
#' rank_results <- get_ranks(counts)
#'
#' # Extract outputs
#' ranks <- rank_results$ranks
#' rank_caps <- rank_results$rank_caps
#'
#' # Use custom tie handling
#' rank_results <- get_ranks(
#' counts,
#' handle_ties = "average"
#' )

# set the generic for get_ranks() methods
setGeneric("get_ranks",
           function(counts,
                    count_caps = compute_geometric_average(counts),
                    handle_ties = "min")
               standardGeneric("get_ranks"))


# set the method for get_ranks() where counts is a matrix
setMethod("get_ranks",
          signature('matrix','ANY', 'ANY'),
          function(counts,
                   count_caps = compute_geometric_average(counts),
                   handle_ties = "min"){

              #call the helper function
              ranks <- get_ranks_from_counts(
                  counts,
                  count_caps = compute_geometric_average(counts),
                  handle_ties = "min")
              return(ranks)
          }
)


# set the method for get_ranks() where counts is a sparseMatrix
setMethod("get_ranks",
          signature('sparseMatrix','ANY', 'ANY'),
          function(counts,
                   count_caps = compute_geometric_average(counts),
                   handle_ties = "min"){

              #call the helper function
              ranks <- get_ranks_from_counts(
                  counts,
                  count_caps = compute_geometric_average(counts),
                  handle_ties = "min")
              return(ranks)
          }
)


# set the method for get_ranks() where counts is a DelayedMatrix
setMethod("get_ranks",
          signature('DelayedMatrix','ANY', 'ANY'),
          function(counts,
                   count_caps = compute_geometric_average(counts),
                   handle_ties = "min"){

              #call the helper function
              ranks <- get_ranks_from_counts(
                  counts,
                  count_caps = compute_geometric_average(counts),
                  handle_ties = "min")
              return(ranks)
          }
)


# set the method for get_ranks() where counts is a data.frame
setMethod("get_ranks",
          signature('data.frame','ANY', 'ANY'),
          function(counts,
                   count_caps = compute_geometric_average(counts),
                   handle_ties = "min"){

              #call the helper function
              ranks <- get_ranks_from_counts(
                  as.matrix(counts),
                  count_caps = compute_geometric_average(counts),
                  handle_ties = "min")
              return(ranks)
          }
)
