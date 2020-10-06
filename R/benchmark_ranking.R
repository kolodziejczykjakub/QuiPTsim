
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

  do.call(rbind, results)
}


evaluate_filtering_results <- function(m, filtering_results, setup, validation_scheme){

  if (setup[["method"]] == "FCBF") {
    setup[["n_kmers"]] <- sum(filtering_results[["score"]])
  }

  # Iterate over numbers of selected k-mers
  lapply(validation_scheme[["n_kmers"]], function(n_kmers) {

    #TODO: switch
    if (setup[["method"]] == "FCBF") {

      X <- as.matrix(m[, filtering_results[["score"]]])
      y <- attr(m, "target")
      df <- data.frame(X, y)

      results <- evaluate_selected_kmers(df, validation_scheme)

    }

    if (setup[["method"]] == "QuiPT") {

      X <- as.matrix(m[, filtering_results[["rank"]] <= n_kmers])
      y <- attr(m, "target")

      df <- data.frame(X, y)

      results <- evaluate_selected_kmers(df, validation_scheme)
    }

  })


  # # one function that builds model etc.
  # results <- evaluate_selected_kmers(df, validation_scheme)

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



