library(QuiPTsim)
library(pbapply)

df <- read.csv("./reduced_alph_enc_alph4_const/alph4_const.csv")

paths <- df[df$l_seq == 10 & df$n_motifs == 2, "path"]
setup <- list(method = "FCBF")
validation_scheme <- list(type = "cv",
                          folds = 5,
                          learner = "classif.log_reg")

m <- read_ngram_matrix(paths[1])
filtering_results <- filter_ngrams(m,setup[["method"]])

X <- as.matrix(m[, filtering_results[["score"]]])
y <- attr(m, "target")
df <- data.frame(X, y)

df
cvReplications <- 1
n_folds = 5
for (i in 1:cvReplications) {

}
folds <- createFolds(y = df$y, k = n_folds)

for (foldNum in 1:n_folds) {

  df_train <- df[-folds[[foldNum]], ]
  df_test <- df[folds[[foldNum]], ]
}

# linear model
logReg <- glm(y ~ ., family = binomial(link = "logit"), data = df_train)

y_true <- df_test[["y"]]
y_proba <- predict(logReg, df_test)
y_pred <- y_proba > 0.5


evaluate_filtering_results <- function(m, filtering_results, setup, validation_scheme){

  # Validation scheme should contain:
  # + metrics
  # + CV details
  # + number of k-mers (if needed)

  # Feature selection methods
  # praznik_filters <- c("MIM", "MRMR", "JMI", "JMIM", "DISR", "NJMIM", "CMIM", "CMI")
  # FSelectorRcpp_measures <- c("infogain", "gainratio", "symuncert")
  # QuiPT
  # FCBF
  # Chi-squared


  if (setup[["method"]] == "FCBF") {
    setup[["n_kmers"]] <- sum(filtering_results[["score"]])
  }

  # Iterate over numbers of selected k-mers
  lapply(setup[["n_kmers"]], function(n_kmers) {

    #TODO: switch
    if (setup[["method"]] == "FCBF") {

      X <- as.matrix(m[, filtering_results[["score"]]])
      y <- attr(m, "target")
      df <- data.frame(X, y)

      results <- evaluate_selected_kmers(df, validation_scheme)

    }

  })


  # # one function that builds model etc.
  # results <- evaluate_selected_kmers(df, validation_scheme)

}

evaluate_selected_kmers <- function(df, validation_scheme) {

  # caret::createFolds
  # Stratified


}


ranking_summary <- function(paths, setup, validation_scheme) {

  results <- pblapply(paths, function(path) {

    # read n-gram matrix
    # TODO: fraction, mixing matrices etc.
    m <- read_ngram_matrix(path)

    # scoring k-mers
    filtering_results <- filter_ngrams(m,setup[["method"]])

    # Cross-validation, bootstraping, holdout etc.
    results <- evaluate_filtering_results(m, filtering_results, setup, validation_scheme)

    }

  })

}



