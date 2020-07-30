library(QuiPTsim)

# Number of sequences
n_seq <- c(6000)

# Replications
reps <- 100

# Number of motifs
n_motifs <- 1:3

# Sequence lengths
l_seq <- 10 * 2 ^ (0:3)

# alphabets
alph4 <- letters[1:4]
alph6 <- letters[1:6]
alphs <- list(alph6, alph4, alph4, alph6, alph4)

# probabilities
weights <- readRDS("./inst/encodingProbs.Rds")
weights[["alph6_const"]] <- list(positive = NULL, negative = NULL)
weights[["alph4_const"]] <- list(positive = NULL, negative = NULL)

# titles
titles <- names(weights)
# paths

paths <- lapply(titles, function(x) paste0("./reduced_alph_enc_", x, "/"))

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
  results <- create_simulation_data(reps, n_seq, l_seq, n_motifs, alphs[[i]],
                                    paths[[i]], titles[[i]], TRUE,
                                    motifProbs = weights[[titles[[i]]]][["positive"]],
                                    seqProbs = weights[[titles[[i]]]][["negative"]],
                                    n = n,
                                    d = d)
  print(paste0(rep("-", 40), collapse=""))
  print(titles[[i]])
  print(results)
  print(paste0(rep("-", 40), collapse=""))
}
