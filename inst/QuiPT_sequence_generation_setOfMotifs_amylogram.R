library(QuiPTsim)

# Number of sequences
n_seq <- c(10)

# Replications
reps <- 10

# Number of motifs
motif_set_size <- 15
n_motifs <- 1:3

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

paths <- lapply(titles, function(x) paste0("./exp4_reduced_alph_enc_", x, "/"))

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
  create_simulation_data_set_of_motifs(reps,
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
                                       d = d)
  print(paste0(rep("-", 40), collapse=""))
  print(titles[[i]])
  print(results)
  print(paste0(rep("-", 40), collapse=""))
  

  
}



# do.call(expand.grid, list(1:15, 1:15, 1:15))
# 
# n_motif <- 15
# validation_size <- 3
# alphabet <- letters[1:6]
# n <- 4
# d <- 6
# motifProbs <- NULL
# motifs <- lapply(1L:n_motif, function(dummy) generate_single_motif(alphabet, n, d, motifProbs))
# motifs_grid <- do.call(expand.grid, rep(list(1:n_motif), validation_size))
# possible_motifs_grid <- motifs_grid[apply(motifs_grid, 1, function(x) length(unique(unlist(x))) == 3), , drop=FALSE]
# 
# tmp <- apply(possible_motifs_grid, 1, function(x) validate_motifs(motifs[unlist(x)], 10))

