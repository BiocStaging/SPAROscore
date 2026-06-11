# helper testing function
make_test_counts <- function() {
    matrix(
        c(
            10, 5,  1,
            2,  8,  3,
            0,  1,  7
        ),
        nrow = 3,
        byrow = TRUE,
        dimnames = list(
            c("gene1", "gene2", "gene3"),
            c("cell1", "cell2", "cell3")
        )
    )
}


test_that("sparoscore matrix dispatch returns correct structure", {

    counts <- make_test_counts()

    res <- sparoscore(
        data = counts,
        signatures = c("gene1", "gene2")
    )

    expect_true(is.matrix(res))
    expect_equal(nrow(res), ncol(counts))
    expect_equal(colnames(res), "SPAROscore")
})


test_that("sparoscore sparseMatrix behaves like matrix", {

    skip_if_not_installed("Matrix")

    counts <- Matrix::Matrix(make_test_counts(), sparse = TRUE)

    res <- sparoscore(
        data = counts,
        signatures = c("gene1", "gene2")
    )

    expect_true(is.matrix(res))
    expect_equal(nrow(res), ncol(counts))
})


test_that("sparoscore data.frame input behaves correctly", {

    df <- as.data.frame(make_test_counts())

    res <- sparoscore(
        data = df,
        signatures = c("gene1", "gene2")
    )

    expect_true(is.matrix(res))
})


test_that("sparoscore handles missing genes safely", {

    counts <- make_test_counts()

    expect_no_error(
        res <- sparoscore(
            data = counts,
            signatures = c("gene1", "geneX")
        )
    )

    expect_true(is.matrix(res))
})


test_that("sparoscore fails when all genes missing", {

    counts <- make_test_counts()

    expect_error(
        sparoscore(
            data = counts,
            signatures = c("geneX", "geneY")
        ),
        regexp = "No signature genes are present"
    )
})


test_that("sparoscore handles multiple signatures", {

    counts <- make_test_counts()

    sigs <- list(
        A = c("gene1", "gene2"),
        B = c("gene2", "gene3")
    )

    res <- sparoscore(
        data = counts,
        signatures = sigs
    )

    expect_true(is.matrix(res))
    expect_equal(colnames(res), names(sigs))
})


test_that("sparoscore Seurat method adds metadata", {

    skip_if_not_installed("Seurat")

    seurat_object <- Seurat::CreateSeuratObject(
        counts = make_test_counts()
    )

    res <- sparoscore(
        data = seurat_object,
        signatures = c("gene1", "gene2")
    )

    meta <- Seurat::FetchData(res, vars = "rank_caps")

    expect_true("rank_caps" %in% colnames(res[[]]))
    expect_true(is.data.frame(res[[]]))
    expect_true(nrow(meta) > 0)
})



test_that("sparoscore SummarizedExperiment works", {

    skip_if_not_installed("SummarizedExperiment")

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = make_test_counts())
    )

    res <- sparoscore(
        data = se,
        signatures = c("gene1", "gene2")
    )

    expect_true("rank_caps" %in% colnames(SummarizedExperiment::colData(res)))
})


test_that("sparoscore fails on invalid assay/layer", {

    skip_if_not_installed("Seurat")

    seurat_object <- Seurat::CreateSeuratObject(
        counts = make_test_counts()
    )

    expect_error(
        sparoscore(
            data = seurat_object,
            signatures = c("gene1"),
            assay = "NON_EXISTENT"
        )
    )
})


test_that("sparoscore fails on misaligned rank_caps", {

    counts <- make_test_counts()

    expect_error(
        sparoscore(
            data = counts,
            signatures = c("gene1", "gene2"),
            rank_caps = c(1, 2)  # missing column
        )
    )
})


test_that("sparoscore fails on NA input", {

    counts <- make_test_counts()
    counts[1,1] <- NA

    expect_error(
        sparoscore(
            data = counts,
            signatures = c("gene1", "gene2")
        )
    )
})



test_that("sparoscore does not silently fallback incorrectly", {

    counts <- make_test_counts()

    # intentionally wrong class injection attempt
    class(counts) <- "fakeClass"

    expect_error(
        sparoscore(
            data = counts,
            signatures = c("gene1")
        )
    )
})





