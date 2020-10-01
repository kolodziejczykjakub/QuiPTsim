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
  expect_identical(list(c("1", "_", "2", "1"), c("1", "1"), c("1", "_", "1")),
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

  expect_equal(structure(c("2", "1", "4", "1", "4", "1", "4", "1", "2", "1",
                           "2", "2", "1", "3", "4", "1", "2", "4", "4", "1"), .Dim = c(2L, 10L),
                         motifs = list(list(c("4", "_", "_", "_", "2", "_", "_","2"))),
                         masks = list(list(c(FALSE, TRUE, FALSE, FALSE, FALSE,
                                             TRUE, FALSE, FALSE, TRUE, FALSE))),
                         target = c(TRUE, FALSE)),
               sequences)

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

  expect_equal(structure(c("4", "4", "3", "1", "3", "1", "3", "3", "4", "3"),
                         motifs = list(c("4", "3"), c("1", "3", "_", "_", "3")),
                         masks = list(c(FALSE,TRUE, TRUE, FALSE, FALSE,
                                        FALSE, FALSE, FALSE, FALSE, FALSE),
                                      c(FALSE, FALSE, FALSE, TRUE, TRUE,
                                        FALSE, FALSE, TRUE, FALSE,FALSE))),
               injected2)
})

test_that("Add motifs - stop", {

  set.seed(1)
  alph <- 1:4
  motifs <- list(as.character(rep(1, 10)), as.character(rep(2, 10)))
  expect_error(add_motifs(motifs, rep(3, 10)), "Given motifs cannot be injected to a sequence!")

})

test_that("count 1-mers", {

  set.seed(1)
  n_seq <- 20
  len <- 10
  alph <- letters[1:20]
  motifs <- generate_motifs(alph, 2)
  results <- generate_sequences(n_seq, len, alph, motifs, 1, n = 1, d = 0)
  expect_equal(apply(results, 1, sum),
               c(7L, 9L, 7L, 8L, 8L, 9L, 8L, 6L, 8L, 9L, 9L, 8L, 6L, 8L, 9L,
                 8L, 7L, 8L, 8L, 9L))
  expect_equal(apply(results, 2, sum),
               c(j = 10L, r = 10L, n = 12L, f = 7L, t = 10L, k = 9L, o = 4L,
                 s = 5L, g = 12L, h = 9L, l = 9L, m = 5L, p = 8L, d = 10L, b = 6L,
                 a = 6L, c = 5L, i = 8L, q = 7L, e = 7L))

})

