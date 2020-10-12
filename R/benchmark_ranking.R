#' High-level wrapper of ranking evaluation tools
#' @inheritParams evaluate_filtering_results
#' @export
#' @importFrom pbapply pblapply
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

    res[["cv.rep"]] = validation_scheme[["cv_reps"]]
    do.call(rbind, res)
  })

  do.call(rbind, results)
}

#' Evaluation of filtered k-mers in ranking model approach
#' @param df data frame containing both selected k-mers and target variable y
#' @param validation_scheme list containing folds, n_kmers, cv_reps - validation setup
#' @export
#' @importFrom caret createFolds
#' @importFrom stats glm
#' @importFrom pROC auc
evaluate_selected_kmers <- function(df, validation_scheme) {

  folds <- createFolds(y = df[["y"]],
                       k = validation_scheme[["folds"]])

  results <- lapply(folds, function(fold_indices) {

    df_train <- df[-fold_indices, ]
    df_test <- df[fold_indices, ]

    evaluate_models(df_train, df_test) #, models = c())

  })

  results <- data.frame(do.call(rbind, results))
  results[["fold"]] <- 1:validation_scheme[["folds"]]
  rownames(results) <- NULL

  results
}

#' Function trains and evaluates model on selected k-mers
#' @export
#' @importFrom glmnet glmnet
#' @importFrom ranger ranger
#' @importFrom e1071 naiveBayes
#' @importFrom class knn
evaluate_models <- function(df_train, df_test) {

  X_train <- df_train[, !names(df_train) %in% "y" ]
  y_train <- df_train[["y"]]

  X_test <- df_test[, !names(df_test) %in% "y" ]
  y_test <- df_test[["y"]]

  # Logistic regression with l1 penalty
  lambdas <- 1:5 * 0.01
  logReg <- glmnet(x = as.matrix(X_train),
                   y = as.numeric(y_train),
                   family = "binomial",
                   lambda = lambdas)

  logReg_probs <- predict(logReg, as.matrix(X_test), type = "response")

  # k-NN
  n_neighbors <- 1:10
  kNN_probs <- do.call(cbind, lapply(n_neighbors, function(neighbors) {
    kNN_classifier <- knn(X_train, X_test,
                          cl = y_train,
                          k = neighbors,
                          prob = TRUE)
    ifelse(y_test == TRUE,
           attr(kNN_classifier, "prob"),
           1 - attr(kNN_classifier, "prob"))
  }))

  # Naive bayes
  naiveBayes_classifier <- naiveBayes(X_train, y_train)
  naiveBayes_probs <- predict(naiveBayes_classifier, X_test, type = "raw")[, 2, drop=FALSE]

  # Random forest - ranger
  numTrees <- c(100, 200, 300)
  rf_probs <- do.call(cbind, lapply(numTrees, function(n) {

    rf_classifier <- ranger(y~.,
                            data = df_train,
                            probability = TRUE,
                            num.trees = n)
    rf_preds <- predict(rf_classifier, df_test)
    rf_preds$predictions[, 1]
  }))

  # Aggregation
  aggregated_results <- list(
    list(probs = logReg_probs,
         model = "lm",
         params = list(lambda = lambdas)),
    list(probs = kNN_probs,
         model = "knn",
         params = list(neighbors = n_neighbors)),
    list(probs = naiveBayes_probs,
         model = "naive bayes",
         params = list(laplace = 0)),
    list(probs = rf_probs,
         model = "rf",
         params = list(num.trees = numTrees))
  )

  results <- lapply(aggregated_results, function(results) {
    res <- apply(results[["probs"]], 2, function(y_pred) {
      c(compute_metrics(y_test, y_pred > 0.5),
        auc = auc(as.numeric(y_test), y_pred))
      })
    data.frame(model = results[["model"]],
               results[["params"]],
               do.call(rbind, res))
    })

  do.call(rbind_results, results)
}

#' functions binds rows of two data frames in a smart way
#' (non overlapping columns are filled with NAs)
#' @param df1
#' @param df2
#' @export
rbind_smart <- function(df1, df2) {
  df1[setdiff(names(df2), names(df1))] <- NA
  df2[setdiff(names(df1), names(df2))] <- NA
  rbind(df1, df2)
}

#' functions binds rows of data frames in a smart way
#' @param ... dataframes to bind
#' @export
rbind_results <- function(...) {
  dfs <- list(...)
  df <- dfs[[1]]

  if (length(dfs) > 1) {
    for (i in 2:length(dfs)) {
      df <- rbind_smart(df, dfs[[i]])
    }
  }

  df
}
