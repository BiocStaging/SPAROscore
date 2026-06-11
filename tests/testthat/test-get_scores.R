library(Matrix)

# helper: small rank matrix + caps
make_test_ranks <- function() {
    ranks <- matrix(
        c(
            1, 2,
            3, 4,
            5, 6
        ),
        nrow = 3,
        byrow = TRUE
    )
    rownames(ranks) <- c("geneA", "geneB", "geneC")
    colnames(ranks) <- c("cell1", "cell2")

    rank_caps <- c(cell1 = 3, cell2 = 3)

    list(ranks = ranks, rank_caps = rank_caps)
}


test_that("get_scores works for character vector signatures", {

    obj <- make_test_ranks()

    res <- get_scores(
        ranks = obj$ranks,
        rank_caps = obj$rank_caps,
        signatures = c("geneA", "geneB")
    )

    expect_true(is.matrix(res))
    expect_equal(nrow(res), ncol(obj$ranks))
    expect_equal(colnames(res), "SPAROscore")
    expect_true(is.numeric(res))
})


test_that("get_scores works for named list signatures", {

    obj <- make_test_ranks()

    sigs <- list(
        sig1 = c("geneA", "geneB"),
        sig2 = c("geneB", "geneC")
    )

    res <- get_scores(
        ranks = obj$ranks,
        rank_caps = obj$rank_caps,
        signatures = sigs
    )

    expect_true(is.matrix(res))
    expect_equal(nrow(res), ncol(obj$ranks))
    expect_equal(colnames(res), names(sigs))
    expect_true(is.numeric(res))
})


test_that("get_scores character vs list output consistency", {

    obj <- make_test_ranks()

    char_res <- get_scores(
        ranks = obj$ranks,
        rank_caps = obj$rank_caps,
        signatures = c("geneA", "geneB")
    )

    list_res <- get_scores(
        ranks = obj$ranks,
        rank_caps = obj$rank_caps,
        signatures = list(sig = c("geneA", "geneB"))
    )

    expect_equal(
        as.numeric(char_res),
        as.numeric(list_res[, "sig"])
    )
})


# test for skip mode in handling missing genes
test_that("get_scores handles partially missing genes safely (skip mode)", {

    obj <- make_test_ranks()

    expect_no_error(
        get_scores(
            ranks = obj$ranks,
            rank_caps = obj$rank_caps,
            signatures = c("geneA", "geneX"),
            handle_missing_genes = "skip"
        )
    )
})


test_that("get_scores fails cleanly when no genes match (skip mode)", {

    obj <- make_test_ranks()

    expect_error(
        get_scores(
            ranks = obj$ranks,
            rank_caps = obj$rank_caps,
            signatures = c("geneX", "geneY"),
            handle_missing_genes = "skip"
        ),
        "No signature genes are present"
    )
})


test_that("GeneSet method works (if GSEABase available)", {

    skip_if_not_installed("GSEABase")

    library(GSEABase)

    obj <- make_test_ranks()

    gs <- GeneSet(setName = "test_set",
                  geneIds = c("geneA", "geneB"))

    res <- get_scores(
        ranks = obj$ranks,
        rank_caps = obj$rank_caps,
        signatures = gs
    )

    expect_true(is.matrix(res))
    expect_equal(nrow(res), ncol(obj$ranks))
    expect_equal(colnames(res), "SPAROscores")
})


test_that("GeneSetCollection method works (if GSEABase available)", {

    skip_if_not_installed("GSEABase")

    library(GSEABase)

    obj <- make_test_ranks()

    gs1 <- GeneSet(setName = "sig1", geneIds = c("geneA", "geneB"))
    gs2 <- GeneSet(setName = "sig2", geneIds = c("geneB", "geneC"))

    gsc <- GeneSetCollection(gs1, gs2)

    res <- get_scores(
        ranks = obj$ranks,
        rank_caps = obj$rank_caps,
        signatures = gsc
    )

    expect_true(is.matrix(res))
    expect_equal(ncol(res), 2)
    expect_equal(colnames(res), names(gsc))
    expect_true(is.numeric(res))
})


test_that("invalid rank_caps names mismatch is tolerated", {

    obj <- make_test_ranks()

    bad_caps <- c(a = 1, b = 2)

    expect_error(
        get_scores(
            ranks = obj$ranks,
            rank_caps = bad_caps,
            signatures = c("geneA", "geneB")
        )
    )
})


# test for breakage in no genes existing
test_that("BREAKS: all signature genes missing", {

    obj <- make_test_ranks()

    expect_error(
        get_scores(
            ranks = obj$ranks,
            rank_caps = obj$rank_caps,
            signatures = c("geneX", "geneY"),
            handle_missing_genes = "skip"
        ),
        regexp = "No signature genes are present"
    )
})


test_that("BREAKS: empty subset after filtering valid genes", {

    obj <- make_test_ranks()

    # valid gene exists in ranks, but not in signature after filtering logic
    expect_error(
        get_scores(
            ranks = obj$ranks,
            rank_caps = obj$rank_caps,
            signatures = setdiff(rownames(obj$ranks), rownames(obj$ranks)),
            handle_missing_genes = "skip"
        )
    )
})


#check if rank_caps have names as cells
test_that("BREAKS: unnamed rank_caps cause misalignment", {

    obj <- make_test_ranks()

    bad_caps <- c(3, 5)  # no names

    expect_error(
        get_scores(
            ranks = obj$ranks,
            rank_caps = bad_caps,
            signatures = c("geneA", "geneB")
        )
    )
})



# check for breakage when rank_caps have less entries
test_that("BREAKS: rank_caps length mismatch", {

    obj <- make_test_ranks()

    expect_error(
        get_scores(
            ranks = obj$ranks,
            rank_caps = c(cell1 = 3),  # missing column
            signatures = c("geneA", "geneB")
        )
    )
})



#signature has only one gene
test_that("BREAKS: single gene signature stability", {

    obj <- make_test_ranks()

    res <- get_scores(
        ranks = obj$ranks,
        rank_caps = obj$rank_caps,
        signatures = c("geneA")
    )

    expect_true(is.matrix(res))
    expect_equal(ncol(res), 1)
})


# check for empty geneset object
test_that("BREAKS: empty GeneSet", {

    skip_if_not_installed("GSEABase")
    library(GSEABase)

    gs <- GeneSet(setName = "empty", geneIds = character(0))

    obj <- make_test_ranks()

    expect_error(
        get_scores(
            ranks = obj$ranks,
            rank_caps = obj$rank_caps,
            signatures = gs
        )
    )
})


#check for duplicate gene signatures
test_that("BREAKS: duplicate genes inflate scoring", {

    obj <- make_test_ranks()

    res <- get_scores(
        ranks = obj$ranks,
        rank_caps = obj$rank_caps,
        signatures = c("geneA", "geneA", "geneB")
    )

    expect_true(is.matrix(res))
})


#check for NAs

test_that("BREAKS: NA values in ranks", {

    obj <- make_test_ranks()
    obj$ranks[1,1] <- NA

    expect_error(
        get_scores(
            ranks = obj$ranks,
            rank_caps = obj$rank_caps,
            signatures = c("geneA", "geneB")
        )
    )
})


# check for dim(0,0) matrices as imputs
test_that("BREAKS: empty rank matrix", {

    ranks <- matrix(numeric(0), nrow = 0, ncol = 0)
    rank_caps <- numeric(0)

    expect_error(
        get_scores(
            ranks = ranks,
            rank_caps = rank_caps,
            signatures = c("geneA")
        )
    )
})


#check for invalid input type for signature
test_that("BREAKS: wrong method dispatch fallback", {

    obj <- make_test_ranks()

    expect_error(
        get_scores(
            ranks = obj$ranks,
            rank_caps = obj$rank_caps,
            signatures = 123  # invalid type
        )
    )
})
