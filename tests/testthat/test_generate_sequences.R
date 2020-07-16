context("Sequences generation")

test_that("Correct motifs generation", {

  set.seed(1)
  single_motif <- generate_single_motif(1)
  three_motifs <- generate_motifs(1:4, 3)
  weighted_motifs <- generate_motifs(1:4, 3, n = 4, d = 1, motifProbs = c(0.7, 0.1, 0.1, 0.1))

  expect_identical(c("1", "1"), single_motif)
  expect_identical(list(c(1, 3, "_", "_",  3),
                        c(3, "_", "_", 3,"_", "_","_", 1),
                        c(2, "_", "_", "_", "_", "_", 2)),
                   three_motifs)
  expect_identical(list(c("3", "2", "1"),
                        c(1, "_", 1),
                        c(4, "_", 1)),
                   weighted_motifs)

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

test_that("Multiple sequence generation", {

  set.seed(42)
  n_seq <- 2
  len <- 10
  alph <- 1L:4
  motifs <- generate_motifs(alph, 2)
  sequences <- simulate_sequences(n_seq, len, alph, motifs, 1, seqProbs = c(0.7, 0.1, 0.1, 0.1))

  expect_true(
    all(rbind(c(3, 4, 4, 4, 1, 1, 2, 4, 2, 2),
              c(2, 1, 1, 1, 2, 1, 4, 3, 4, 1)) == sequences))
  expect_equal(attr(sequences, "motifs")[[1]][[1]],
               c(4, "_", "_", "_", 2, "_", "_", 2))
  expect_equal(attr(sequences, "masks")[[1]][[1]],
               c(F, F, T, F, F, F, T, F, F, T))
  expect_equal(attr(sequences, "target"),
               c(TRUE, FALSE))

})
test_that("Correct motif injection", {

  set.seed(1)

  injected1 <- add_motifs(list(c(1, "_", 1), c(1, 1)), c(2, 2, 3, 4))
  expect_true(all(as.character(c(1, 1, 1, 4)) == injected1))


  set.seed(1)
  alph <- 1:4
  motifs <- generate_motifs(alph, 2)
  injected2 <- add_motifs(motifs, simulate_single_sequence(10, alph))

  injected2_true <- c(1, 3, 3, 1, 3, 1, 4, 3, 2, 2)
  attr(injected2_true, "motifs") <- list(c(4, 3), c(1, 3, "_", "_", 3))
  attr(injected2_true, "masks") <- list(c(F, F, F, F, F, F, T, T, F, F),
                                        c(T, T, F, F, T, F, F, F, F, F))

  expect_true(all(injected2 == injected2_true))

  expect_true(all(attr(injected2, "motifs")[[1]] == attr(injected2_true, "motifs")[[1]]))
  expect_true(all(attr(injected2, "motifs")[[2]] == attr(injected2_true, "motifs")[[2]]))

  expect_true(all(attr(injected2, "masks")[[1]] == attr(injected2_true, "masks")[[1]]))
  expect_true(all(attr(injected2, "masks")[[2]] == attr(injected2_true, "masks")[[2]]))
})

test_that("Count n-grams assertion", {

})
