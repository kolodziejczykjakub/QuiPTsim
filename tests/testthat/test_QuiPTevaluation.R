context("QuiPT evaluation")

test_that("Correct simulation details read", {
  results <- read_simulation_data("./", "example_QuiPT_simulation_data_result")

  expect_equal(dim(results), c(120L, 11L))
  expect_equal(attributes(results),
               list(names = c("X", "replication", "n_seq", "l_seq", "n_motifs",
                              "n", "d", "seqProbs", "motifProbs", "path", "md5sum"), class = "data.frame",
                    row.names = 1:120, details = list(replication = 1:10, n_seq = 600L,
                                                      l_seq = c(10L, 20L, 40L, 80L), n_motifs = 1:3, n = 4L,
                                                      d = 6L, seqProbs = NA, motifProbs = NA)))



})

test_that("Correct data read", {
  results <- read_ngram_matrix("./example_QuiPT_simulation_data_result.Rds")

  expect_equal(dim(results), c(600L, 20777L))
  expect_equal(lapply(attributes(results), length),
               list(names = 6L, class = 1L, sequences = 6000L, motifs = 300L,
                    masks = 300L, target = 600L))

})

test_that("Correct data subsampling", {

  set.seed(42)
  results <- read_ngram_matrix("./example_QuiPT_simulation_data_result.Rds",
                               n = 300,
                               fraction = 0.8)
  expect_equal(dim(results), c(300L, 20777L))
  expect_equal(lapply(attributes(results), length),
               list(names = 6L, class = 1L, sequences = 3000L, motifs = 240L,
                    masks = 240L, target = 300L))

})

