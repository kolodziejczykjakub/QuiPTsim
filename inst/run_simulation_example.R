library(QuiPTsim)

alph <- 1L:4
reps <- 3
n_seq <- c(10)
l_seq <- c(10)
n_motifs <- 1
motifProbs <- c(0.7, 0.1, 0.1, 0.1)
seqProbs = c(0.7, 0.1, 0.1, 0.1)
n <- 4
d <- 4
title <- "example_01"
path <- "./results_01/"
if (!file.exists(path)) dir.create(path)

results <- create_simulation_data(reps, n_seq, l_seq, n_motifs, alph,
                                  path, title, TRUE,
                                  motifProbs = motifProbs,
                                  seqProbs = seqProbs,
                                  n = n,
                                  d = d)
