context("QuiPT simulation")

test_that("Correct result of create_simulation_data", {

  set.seed(1)
  alph <- 1L:4
  reps <- 3
  n_seq <- c(10)
  l_seq <- c(10)
  n_motifs <- 1
  path <- "./"
  results <- create_simulation_data(reps, n_seq, l_seq, n_motifs, alph,
                                    path, "SEQ", FALSE,
                                    motifProbs = c(0.7, 0.1, 0.1, 0.1),
                                    seqProbs = c(0.7, 0.1, 0.1, 0.1),
                                    n = 4, d = 4)
  results_true <- readRDS("./test_create_simulation_data_df.Rds")

  expect_equal(results_true, results)
})
