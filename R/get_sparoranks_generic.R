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
#' Can natively handle matrix, data.frame, sparsematrix, delayedmatrix.
#' Hence, can handle count assays from Seurat objects,
#' SummarisedExperiment, SingleCellExperiment objects, SpatialExperiment objects
#' and many mare.
#'
#'
#' @param handle_ties A character string specifying how ties are treated.
#' "min" by default. Can be "min", "max", "average", or "random"
#'
#'
#' @param rank_cap_metric A character string specifying the type of rank cap.
#' Default: "geometric_average" cap is the geometric average of (1 + count)
#' "manual_override" cap is user mentioned constant expression values.
#' "manual_override" requires manual_expr_caps to be NOT NULL
#'
#'
#' @param manual_expr_caps used only if rank_cap_metric = "manual_override"
#' A numeric vector of length equals total columns in count_data. Contains rank
#' caps for each sample/cell/spot individually
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
#' sparoranks <- get_sparoranks(counts_data,
#' handle_ties = "min", rank_cap_metric = "geometric_average")
#'

# set the generic for get_sparoranks() methods
setGeneric("get_sparoranks",
           function(counts_data,
                    handle_ties = "min",
                    rank_cap_metric = "geometric_average",
                    manual_expr_caps = NULL)
               standardGeneric("get_sparoranks"))


# set the method for get_sparoranks() where counts is a matrix
setMethod("get_sparoranks",
          signature('matrix','ANY', 'ANY', 'ANY'),
          function(counts_data,
                   handle_ties = "min",
                   rank_cap_metric = "geometric_average",
                   manual_expr_caps = NULL){

              #call the helper function
              sparoranks <- get_sparoranks_from_counts(counts_data,
                            handle_ties = "min",
                            rank_cap_metric = "geometric_average",
                            manual_expr_caps = NULL)
              return(sparoranks)
          }
)


# set the method for get_sparoranks() where counts is a sparseMatrix
setMethod("get_sparoranks",
          signature('sparseMatrix','ANY', 'ANY', 'ANY'),
          function(counts_data,
                   handle_ties = "min",
                   rank_cap_metric = "geometric_average",
                   manual_expr_caps = NULL){

              #call the helper function
              sparoranks <- get_sparoranks_from_counts(counts_data,
                            handle_ties = "min",
                            rank_cap_metric = "geometric_average",
                            manual_expr_caps = NULL)
              return(sparoranks)
          }
)


# set the method for get_sparoranks() where counts is a DelayedMatrix
setMethod("get_sparoranks",
          signature('DelayedMatrix','ANY', 'ANY', 'ANY'),
          function(counts_data,
                   handle_ties = "min",
                   rank_cap_metric = "geometric_average",
                   manual_expr_caps = NULL){

              #call the helper function
              sparoranks <- get_sparoranks_from_counts(counts_data,
                            handle_ties = "min",
                            rank_cap_metric = "geometric_average",
                            manual_expr_caps = NULL)
              return(sparoranks)
          }
)


# set the method for get_sparoranks() where counts is a data.frame
setMethod("get_sparoranks",
          signature('data.frame','ANY', 'ANY', 'ANY'),
          function(counts_data,
                   handle_ties = "min",
                   rank_cap_metric = "geometric_average",
                   manual_expr_caps = NULL){

              #call the helper function
              sparoranks <- get_sparoranks_from_counts(as.matrix(counts_data),
                            handle_ties = "min",
                            rank_cap_metric = "geometric_average",
                            manual_expr_caps = NULL)
              return(sparoranks)
          }
)
