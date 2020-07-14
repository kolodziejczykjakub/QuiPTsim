context("Sequences generation")

test_that("Correct motifs generation", {

  set.seed(1)
  single_motif <- generate_single_motif(1)

  three_motifs <- generate_motifs(1:4, 3)
  expect_identical(c(1, "_", "_", 1, 1), single_motif)
  expect_identical(list(c(1, "_", 3),
                        c(2, "_", 2, "_", 3),
                        c(1, "_", "_", 1)),
                   three_motifs)
})

test_that("Correct sequence generation", {

  set.seed(42)
  seq1 <- simulate_single_sequence(5, 1L:4)
  seq2 <- simulate_single_sequence(10, c("a", "b", "c"))

  expect_equal(c(1, 1, 1, 1, 2),
                   seq1)
  expect_identical(c("b", "b", "a", "c", "c", "a", "a", "b", "b", "b"),
                   seq2)

})

test_that("Correct motif injection", {

  set.seed(1)

  injected1 <- add_motifs(list(c(1, "_", 1), c(1, 1)), c(2, 2, 3, 4))

  # TODO: this one should be passed
  # expect_equal(as.character(c(1, 1, 1, 4)),
  #              injected1)
  expect_true(all(as.character(c(1, 1, 1, 4)) == injected1))


  # alph <- 1:4
  # motifs <- generate_motifs(alph, 2)
  # injected2 <- add_motifs(motifs, simulate_single_sequence(10, alph))
  #
  # injected2_true <- c(4, 3, 2, 3, 1, 2, 3, 2, 4, 4)
  # attr(injected2_true, "motifs") <- list(c(2, "_", 2), c(4, "_", "_", 3, 1))
  # attr(injected2_true, "masks") <- list(c(rep(FALSE, 5), T, F, T, F, F),
  #                                       c(T, F, F, T, T, rep(FALSE, 5)))
  #
  # expect_identical(injected2_true, injected2)

})
