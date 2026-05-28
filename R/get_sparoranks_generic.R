#' Compute ranks for SPAROscore
#'
#' get_sparoranks() is a generic function that take counts matrix data as input
#' and returns column wise ascending ranks matrix.
#' It uses the S4 generic MatrixGenerics::colRanks() and can handle a wide range
#' of input formats (see Arguments)
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
#' @export
#'
#' @examples
#'
#' # For counts data in the form of matrix, sparsematrix or delayedmatrix
#' sparoranks <- get_sparoranks(counts_data)
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
