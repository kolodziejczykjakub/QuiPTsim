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
# filtering_results <- filter_ngrams(m, setup[["method"]])
# evaluate_filtering_results(m, filtering_results, setup, validation_scheme)

# QuiPT
setup <- list(method = "QuiPT")

models_details <- list(
  list(model = "lm",
       param_name = "lambda",
       param_value = NULL),
  list(model = "knn",
       param_name = "neighbors",
       param_value = 1:2),
  list(model = "rf",
       param_name = "num.trees",
       param_value = 1:2),
  list(model = "naive bayes",
       param_name = "laplace",
       param_value = 1:2)
)
validation_scheme <- list(type = "cv",
                          folds = 5,
                          n_kmers = c(2, 5, 100),
                          cv_reps = 2,
                          models_details = models_details)
# results <- ranking_summary(paths[1:2], setup, validation_scheme)
filtering_results <- filter_ngrams(m, setup[["method"]])
res <- evaluate_filtering_results(m, filtering_results, setup, validation_scheme)


setup <- list(method = "FCBF")
filtering_results <- filter_ngrams(m, setup[["method"]])

kmers_for_nonranking_methods(filtering_results, setup[["method"]], c(0.001, 0.01, 0.05))
