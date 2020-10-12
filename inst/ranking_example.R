library(QuiPTsim)
library(pbapply)
library(caret)
library(pROC)

df <- read.csv("./reduced_alph_enc_alph4_const/alph4_const.csv")

paths <- df[df$l_seq == 10 & df$n_motifs == 2, "path"]

# FCBF
setup <- list(method = "FCBF")
validation_scheme <- list(type = "cv",
                          folds = 5,
                          cv_reps = 2)

m <-read_ngram_matrix(paths[1])
filtering_results <- filter_ngrams(m, setup[["method"]])
evaluate_filtering_results(m, filtering_results, setup, validation_scheme)

# QuiPT
setup <- list(method = "QuiPT")
validation_scheme <- list(type = "cv",
                          folds = 5,
                          n_kmers = c(2),
                          cv_reps = 2)
m <-read_ngram_matrix(paths[1])
filtering_results <- filter_ngrams(m, setup[["method"]])
evaluate_filtering_results(m, filtering_results, setup, validation_scheme)



# # Feature selection methods from praznik package
# praznik_filters <- c("MIM", "MRMR", "JMI", "JMIM", "DISR", "NJMIM", "CMIM", "CMI")
# FSelectorRcpp_measures <- c("infogain", "gainratio", "symuncert")
#
# if (!(feature_selection_method %in% c("QuiPT",
#                                       "FCBF",
#                                       "Chi-squared",
#                                       FSelectorRcpp_measures,
#                                       praznik_filters)))

x <- filter_ngrams(m, "infogain")


# QuiPT, Chi-squared
