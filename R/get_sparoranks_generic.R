#' Compute Gene Ranks for SPAROscore
#'
#' Generate gene rank matrices and rank caps for SPAROscore calculations.
#'
#' \emph{get_sparoranks()} computes column-wise ascending gene ranks  and the
#' rank caps from count data and returns the inputs required by
#' \emph{get_sparoscores()}. The function uses MatrixGenerics::colRanks()
#' internally and supports a variety of matrix-like input formats,
#' including data.frames, matrices, sparse matrices, and delayed matrices.
#'
#' By default, ties are resolved using the minimum rank. Rank caps are
#' determined from the geometric mean count of each column, although custom
#' count caps can also be supplied.
#'
#' @param counts_data A matrix-like object containing gene expression counts,
#' where rows correspond to genes and columns correspond to samples, cells,
#' or spatial domains.
#'
#' Supported input types include:
#' * matrix
#' * sparseMatrix
#' * DelayedMatrix
#' * data.frame
#'
#' @param count_caps An optional numeric vector specifying count cap values for
#' each column of counts_data. These values are used to determine the
#' corresponding rank caps. By default, the geometric mean count of
#' each column is used.
#'
#'
#' @param handle_ties Character string specifying how tied count values are
#' ranked. Passed to MatrixGenerics::colRanks().
#' Supported options are:
#' \describe{
#' \item{"min"}{Assign the minimum rank to tied values (default).}
#' \item{"max"}{Assign the maximum rank to tied values.}
#' \item{"average"}{Assign the average rank to tied values.}
#' \item{"random"}{Break ties at random.}
#' }
#'
#'
#' @returns A list containing:
#' \describe{
#' \item{sparoranks}{A matrix of gene ranks with the same dimensions as
#' counts_data. Integer-valued unless handle_ties = "average", in which
#' case ranks may be numeric.}
#' \item{rank_caps}{A numeric vector containing the rank cap associated
#' with each column.}
#' }
#'
#'
#' @details
#' Gene ranks are computed independently for each column using ascending order
#' of expression values. The resulting rank matrix and rank caps serve as the
#' primary inputs to \emph{get_sparoscores()}, which computes SPAROscores
#' based on deviations from the rank caps.
#'
#' @export
#'
#' @examples
#' # Compute SPARO ranks from a dense, sparse, or delayed matrix
#' sparorank_output <- get_sparoranks(counts_data)
#'
#' # Access the rank matrix and rank caps
#' sparoranks <- sparorank_output$sparoranks
#' rank_caps <- sparorank_output$rank_caps
#'
#'

# set the generic for get_sparoranks() methods
setGeneric("get_sparoranks",
           function(counts_data,
                    count_caps = compute_geometric_average(counts_data),
                    handle_ties = "min")
               standardGeneric("get_sparoranks"))


# set the method for get_sparoranks() where counts is a matrix
setMethod("get_sparoranks",
          signature('matrix','ANY', 'ANY'),
          function(counts_data,
                   count_caps = compute_geometric_average(counts_data),
                   handle_ties = "min"){

              #call the helper function
              sparoranks <- get_sparoranks_from_counts(
                  counts_data,
                  count_caps = compute_geometric_average(counts_data),
                  handle_ties = "min")
              return(sparoranks)
          }
)


# set the method for get_sparoranks() where counts is a sparseMatrix
setMethod("get_sparoranks",
          signature('sparseMatrix','ANY', 'ANY'),
          function(counts_data,
                   count_caps = compute_geometric_average(counts_data),
                   handle_ties = "min"){

              #call the helper function
              sparoranks <- get_sparoranks_from_counts(
                  counts_data,
                  count_caps = compute_geometric_average(counts_data),
                  handle_ties = "min")
              return(sparoranks)
          }
)


# set the method for get_sparoranks() where counts is a DelayedMatrix
setMethod("get_sparoranks",
          signature('DelayedMatrix','ANY', 'ANY'),
          function(counts_data,
                   count_caps = compute_geometric_average(counts_data),
                   handle_ties = "min"){

              #call the helper function
              sparoranks <- get_sparoranks_from_counts(
                  counts_data,
                  count_caps = compute_geometric_average(counts_data),
                  handle_ties = "min")
              return(sparoranks)
          }
)


# set the method for get_sparoranks() where counts is a data.frame
setMethod("get_sparoranks",
          signature('data.frame','ANY', 'ANY'),
          function(counts_data,
                   count_caps = compute_geometric_average(counts_data),
                   handle_ties = "min"){

              #call the helper function
              sparoranks <- get_sparoranks_from_counts(
                  as.matrix(counts_data),
                  count_caps = compute_geometric_average(counts_data),
                  handle_ties = "min")
              return(sparoranks)
          }
)
