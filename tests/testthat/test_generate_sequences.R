context("Sequences generation")

test_that("Correct motifs generation",{

  set.seed(1)
  single_motif <- generate_single_motif(1)

  expect_identical(c(1, "_", "_", 1, 1), single_motif)
})
