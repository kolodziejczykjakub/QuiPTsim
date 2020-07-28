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
alph1 <- letters[1:4]
alph2 <- letters[1:6]
alph3 <- letters[1:8]
alph4 <- letters[1:20]

# probabilities
weights <- readRDS("./inst/weights.Rds")
length(weights$positive)

# titles
t1 <- "results_alph4elements"
t2 <- "results_alph6elements"
t3 <- "results_alph8elements"
t4 <- "results_alph20elements"

# paths
path1 <- paste0("./", t1, "/", collapse="")
path2 <- paste0("./", t2, "/", collapse="")
path3 <- paste0("./", t3, "/", collapse="")
path4 <- paste0("./", t4, "/", collapse="")

for (p in c(path1, path2, path3, path4)) {
  if (!file.exists(p)) dir.create(p)
}
# 4-elements with 6 gaps motifs
n <- 4
d <- 6

# probability vectors
motifProbs <- NULL
seqProbs <- NULL

paths <- list(path1, path2, path3, path4)
titles <- list(t1, t2, t3, t4)
alphs <- list(alph1, alph2, alph3, alph4)

for (i in 1:4) {
  results <- create_simulation_data(reps, n_seq, l_seq, n_motifs, alphs[[i]],
                                    paths[[i]], titles[[i]], TRUE,
                                    motifProbs = motifProbs,
                                    seqProbs = seqProbs,
                                    n = n,
                                    d = d)
  print(paste0(rep("-", 40), collapse=""))
  print(titles[[i]])
  print(paths[[i]])
  print(alphs[[i]])
  print(results)
  print(paste0(rep("-", 40), collapse=""))
}
