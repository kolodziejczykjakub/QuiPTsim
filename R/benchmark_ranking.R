#' Evaluation of filtered k-mers in ranking model approach
#' @param df data frame containing both selected k-mers and target variable y
#' @param validation_scheme list containing folds, n_kmers, cv_reps - validation setup
#' @export
#' @importFrom caret createFolds
#' @importFrom stats glm
evaluate_selected_kmers <- function(df, validation_scheme) {

  folds <- createFolds(y = df[["y"]],
                       k = validation_scheme[["folds"]])

  results <- lapply(folds, function(fold_indices) {

    df_train <- df[-fold_indices, ]
    df_test <- df[fold_indices, ]

    # linear model
    logReg <- glm(y ~ ., family = binomial(link = "logit"), data = df_train)

    y_true <- df_test[["y"]]
    y_proba <- predict(logReg, df_test)
    y_pred <- y_proba > 0.5


    c(compute_metrics(y_true, y_pred),
                 auc = auc(y_true, y_proba))
  })

  results <- data.frame(do.call(rbind, results))
  results[["fold"]] <- 1:validation_scheme[["folds"]]
  rownames(results) <- NULL

  results
}



#' Wrapper for `evaluate_selected_kmers` functions
#' @param m k-mer matrix
#' @param filtering_results data.frame containing ranking of top k-mers
#' @param setup filter usage setup
#' @param validation_scheme list containing folds, n_kmers, cv_reps - validation setup
#' @export
evaluate_filtering_results <- function(m, filtering_results, setup, validation_scheme){

  if (setup[["method"]] == "FCBF") {
    validation_scheme[["n_kmers"]] <- sum(filtering_results[["score"]])
  }

  # Repeated k-fold cross validation
  results <- lapply(1:validation_scheme[["cv_reps"]], function(cv_rep) {

    # CV repeated for each number of selected k-mers
    res <- lapply(validation_scheme[["n_kmers"]], function(n_kmers) {

      X <- as.matrix(m[, filtering_results[["rank"]] <= n_kmers])
      y <- attr(m, "target")

      df <- data.frame(X, y)
      res <- evaluate_selected_kmers(df, validation_scheme)
      res[["n_kmers"]] <- n_kmers

      res
    })

    do.call(rbind, res)
  })

  results
}

#' High-level wrapper of ranking evaluation tools
#' @inheritParams evaluate_filtering_results
#' @export
ranking_summary <- function(paths, setup, validation_scheme) {

  results <- pblapply(paths, function(path) {


    # TODO: fraction, mixing matrices etc.
    m <- read_ngram_matrix(path)

    # scoring k-mers
    filtering_results <- filter_ngrams(m,setup[["method"]])

    # Cross-validation, bootstraping, holdout etc.
    results <- evaluate_filtering_results(m, filtering_results, setup, validation_scheme)
  })

  results
}



