library(QuiPTsim)
library(pbapply)
library(caret)
library(pROC)

df <- read.csv("./reduced_alph_enc_alph4_const/alph4_const.csv")

paths <- df[df$l_seq == 10 & df$n_motifs == 2, "path"]

# FCBF
setup <- list(method = "FCBF")
validation_scheme <- list(type = "cv",
                          folds = 5)

m <-read_ngram_matrix(paths[1])
evaluate_filtering_results(m, filtering_results, setup, validation_scheme)

# QuiPT
setup <- list(method = "QuiPT")
validation_scheme <- list(type = "cv",
                          folds = 5,
                          n_kmers = c(5, 10, 15, 20, 25, 40, 50, 100))
m <-read_ngram_matrix(paths[1])
filtering_results <- filter_ngrams(m, setup[["method"]])
evaluate_filtering_results(m, filtering_results, setup, validation_scheme)



