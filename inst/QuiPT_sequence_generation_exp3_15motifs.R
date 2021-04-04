library(QuiPTsim)
library(progress)
library(tools)
# Number of sequences
n_seq <- c(1200)

# Replications
reps <- 10

# Number of motifs
motif_set_size <- 15
n_motifs <- 1:2

# Sequence lengths
l_seq <- 10

# alphabets
alph6 <- letters[1:6]
alphs <- list(alph6)

# probabilities
weights <- readRDS("./inst/encodingProbs.Rds")
weights <- weights[1]
# titles
titles <- names(weights)
# paths

paths <- lapply(titles, function(x) paste0("./exp3_reduced_alph_enc_15motifs_", x, "/"))

for (p in paths) {
  if (!file.exists(p)) dir.create(p)
}

# 4-elements with 6 gaps motifs
n <- 4
d <- 6

# simulation

# weights
# paths
# probsNames
# alphabet = alph20
# probVectors

for (i in 1:length(paths)) {
  results <- create_simulation_data_set_of_motifs(reps,
                                       n_seq,
                                       l_seq,
                                       motif_set_size,
                                       n_motifs,
                                       alphs[[i]],
                                       paths[[i]],
                                       titles[[i]],
                                       motifProbs = weights[[titles[[i]]]][["positive"]],
                                       seqProbs = weights[[titles[[i]]]][["negative"]],
                                       n = n,
                                       d = d,
                                       save_files=TRUE)
  print(paste0(rep("-", 40), collapse=""))
  print(titles[[i]])
  print(results)
  print(paste0(rep("-", 40), collapse=""))
}
