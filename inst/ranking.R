library(QuiPTsim)
library(parallel)

df <- read.csv("./reduced_alph_enc_alph4_const/alph4_const.csv")
paths <- as.character(df[df$l_seq == 10 & df$n_motifs == 3, "path"])

setups <- list(
  list(method="QuiPT"),
  list(method="Chi-squared"),
  list(method="FCBF"),
  list(method="gainratio"),
  list(method="MRMR")
)

paths <- paths[1:2]
results <- mclapply(setups, function(setup) {
  print(setup)
  create_benchmark_data(paths, setup)
}, mc.cores = 4)

