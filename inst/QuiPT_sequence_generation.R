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
alph20 <- letters[1:20]

# probabilities
weights <- readRDS("./inst/weights.Rds")

probVectors <- list(list(motifProbs = NULL,
                         seqProbs = NULL))
for (n in names(weights$positive)) {
  probVectors[[n]] <- list(motifProbs = weights[["positive"]][[n]],
                           seqProbs = weights[["negative"]][[n]])
}

# titles
probsNames <- lapply(c("const", names(probVectors[2:6])), function(x) paste0("alph20_prob_", x))
names(probVectors) <- probsNames
# paths
paths <- lapply(probsNames, function(x) paste0("./", x, "/"))

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

for (i in 1:6) {
  results <- create_simulation_data(reps, n_seq, l_seq, n_motifs, alph20,
                                    paths[[i]], probsNames[[i]], TRUE,
                                    motifProbs = probVectors[[i]][["motifProbs"]],
                                    seqProbs = probVectors[[i]][["seqProbs"]],
                                    n = n,
                                    d = d)
  print(paste0(rep("-", 40), collapse=""))
  print(probNames[[i]])
  print(results)
  print(paste0(rep("-", 40), collapse=""))
}
