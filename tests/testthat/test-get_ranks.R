library(testthat)
library(Matrix)
library(DelayedArray)


# test general structure of the rank results
test_that("get_ranks returns expected structure", {

    counts <- matrix(
        c(
            10, 5,
            20, 1,
            30, 0
        ),
        nrow = 3,
        byrow = TRUE
    )

    rownames(counts) <- c("G1", "G2", "G3")
    colnames(counts) <- c("Cell1", "Cell2")

    result <- get_ranks(counts)

    expect_type(result, "list")
    expect_named(result, c("ranks", "rank_caps"))

    expect_equal(
        dim(result$ranks),
        dim(counts)
    )

    expect_equal(
        rownames(result$ranks),
        rownames(counts)
    )

    expect_equal(
        colnames(result$ranks),
        colnames(counts)
    )

    expect_length(
        result$rank_caps,
        ncol(counts)
    )
})


# test logic of the rank results
test_that("get_ranks correctly ranks genes in descending order", {

    counts <- matrix(
        c(
            100,
            50,
            10
        ),
        ncol = 1
    )

    rownames(counts) <- c("G1", "G2", "G3")
    colnames(counts) <- c("Cell1")

    result <- get_ranks(
        counts,
        count_caps = 20
    )

    expect_equal(
        unname(result$ranks[,1]),
        c(1,2,4)
    )

    expect_equal(
        unname(result$rank_caps),
        3
    )
})


# check if user provided custom count caps works
test_that("user supplied count_caps affects rank_caps", {

    counts <- matrix(
        c(
            10, 10,
            5, 5,
            1, 1
        ),
        nrow = 3
    )

    rownames(counts) <- c("G1", "G2", "G3")
    colnames(counts) <- c("Cell1", "Cell2")

    result1 <- get_ranks(
        counts,
        count_caps = c(7, 7)
    )

    result2 <- get_ranks(
        counts,
        count_caps = c(100, 100)
    )

    expect_false(
        identical(result1$rank_caps,
                  result2$rank_caps)
    )
})


#check if tie handling works well
test_that("invalid tie method throws error", {

    counts <- matrix(
        c(1, 2, 3),
        ncol = 1
    )

    rownames(counts) <- c("G1", "G2", "G3")
    colnames(counts) <- c("Cell1")

    expect_error(
        get_ranks(
            counts,
            handle_ties = "foobar"
        ),
        "Invalid value provided for handle_ties"
    )
})


test_that("handle_ties is passed through correctly", {

    counts <- matrix(
        c(10, 10, 5),
        ncol = 1
    )

    rownames(counts) <- c("G1", "G2", "G3")
    colnames(counts) <- "Cell1"

    result_min <- get_ranks(counts, handle_ties = "min")
    result_max <- get_ranks(counts, handle_ties = "max")

    expect_false(
        all(result_min$ranks == result_max$ranks)
    )
})


# check if method works well for dense matrices
test_that("matrix method works", {

    counts <- matrix(
        c(1, 2, 3, 4),
        nrow = 2
    )

    rownames(counts) <- c("G1","G2")
    colnames(counts) <- c("Cell1", "Cell2")

    result <- get_ranks(counts)

    expect_named(
        result,
        c("ranks", "rank_caps")
    )
})


# check if method works well for sparse matrices
test_that("sparseMatrix method works", {

    counts <- Matrix(
        matrix(
            c(
                1, 0,
                0, 2
            ),
            nrow = 2
        ),
        sparse = TRUE
    )

    rownames(counts) <- c("G1","G2")
    colnames(counts) <- c("Cell1", "Cell2")

    result <- get_ranks(counts)

    expect_named(
        result,
        c("ranks", "rank_caps")
    )
})

# check if method works well for delayed matrices
test_that("DelayedMatrix method works", {

    counts <- DelayedArray(
        matrix(
            c(
                1, 2,
                3, 4
            ),
            nrow = 2
        )
    )

    rownames(counts) <- c("G1","G2")
    colnames(counts) <- c("Cell1", "Cell2")

    result <- get_ranks(counts)

    expect_named(
        result,
        c("ranks", "rank_caps")
    )
})


#test if invalid cout caps are producing error
test_that("count_caps length must match number of columns", {

    counts <- matrix(
        c(1, 2, 3, 4),
        nrow = 2
    )

    rownames(counts) <- c("G1","G2")
    colnames(counts) <- c("Cell1", "Cell2")

    expect_error(
        get_ranks(
            counts,
            count_caps = c(1, 2, 3)
        ),
        "count_caps"
    )
})


# test scalar stability
test_that("1 gene x 1 cell matrix behaves consistently", {

    counts <- matrix(5, nrow = 1)
    rownames(counts) <- "G1"
    colnames(counts) <- "Cell1"

    result <- get_ranks(counts)

    # must still be a matrix (not vector)
    expect_true(is.matrix(result$ranks))

    expect_equal(dim(result$ranks), c(1, 1))
    expect_equal(rownames(result$ranks), "G1")
    expect_equal(colnames(result$ranks), "Cell1")

    expect_equal(result$rank_caps, 1)
})


# test single row stability
test_that("1 gene across many cells preserves structure", {

    counts <- matrix(c(10, 5, 2), nrow = 1)
    rownames(counts) <- "G1"
    colnames(counts) <- c("C1", "C2", "C3")

    result <- get_ranks(counts)

    expect_equal(dim(result$ranks), c(1, 3))
    expect_equal(rownames(result$ranks), "G1")
    expect_equal(colnames(result$ranks), c("C1", "C2", "C3"))

    expect_equal(length(result$rank_caps), 3)
})


# test single column stability
test_that("multiple genes in 1 cell preserves matrix structure", {

    counts <- matrix(c(10, 5, 1), ncol = 1)
    rownames(counts) <- c("G1", "G2", "G3")
    colnames(counts) <- "Cell1"

    result <- get_ranks(counts)

    expect_equal(dim(result$ranks), c(3, 1))
    expect_equal(rownames(result$ranks), c("G1", "G2", "G3"))
    expect_equal(colnames(result$ranks), "Cell1")

    expect_equal(length(result$rank_caps), 1)
})


# check dimnames of sparsematrix stability
test_that("sparseMatrix preserves dimnames and structure", {

    counts <- Matrix::Matrix(
        c(0, 1, 0,
          2, 0, 3),
        nrow = 2,
        sparse = TRUE
    )

    rownames(counts) <- c("G1", "G2")
    colnames(counts) <- c("C1", "C2", "C3")

    result <- get_ranks(counts)

    expect_s4_class(counts, "sparseMatrix")
    expect_equal(dim(result$ranks), c(2, 3))

    expect_equal(rownames(result$ranks), c("G1", "G2"))
    expect_equal(colnames(result$ranks), c("C1", "C2", "C3"))
})



# test dimnames of delayed matrix stability
test_that("DelayedArray 1x1 stability", {

    counts <- DelayedArray::DelayedArray(matrix(7, 1, 1))
    rownames(counts) <- "G1"
    colnames(counts) <- "C1"

    result <- get_ranks(counts)

    expect_true(is.matrix(result$ranks))
    expect_equal(dim(result$ranks), c(1, 1))
    expect_equal(rownames(result$ranks), "G1")
    expect_equal(colnames(result$ranks), "C1")
})
