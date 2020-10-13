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

  modele <- list(
    list(model = "lm",
         param_name = "lambda",
         param_value = 1:2),
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

  models_probs <- lapply(modele, function(m)
    build_model(X_train, y_train, X_test, y_test, m[["param_value"]], m[["model"]]))

  results <- lapply(1:length(modele), function(i) {

    res <- apply(models_probs[[i]], 2, function(y_pred) {
      c(compute_metrics(y_test, y_pred > 0.5),
        auc = auc(as.numeric(y_test), y_pred))
    })

    data.frame(model = modele[[i]][["model"]],
               param = modele[[i]][["param_name"]],
               value = modele[[i]][["param_value"]],
               do.call(rbind, res),
               row.names = NULL)
  })

  do.call(rbind, results)
}

#' function builds a model and predicts out of fold probabilites
#' @export
#' @importFrom glmnet glmnet
#' @importFrom ranger ranger
#' @importFrom e1071 naiveBayes
#' @importFrom class knn
build_model <- function(X_train, y_train, X_test, y_test, param, method) {

  switch(method,
         "lm" = {

           logReg <- glmnet(x = as.matrix(X_train),
                            y = as.numeric(y_train),
                            family = "binomial",
                            lambda = param)

           predict(logReg, as.matrix(X_test), type = "response")

         },
         "knn" = {

           do.call(cbind, lapply(param, function(neighbors) {
             kNN_classifier <- knn(train = X_train, test = X_test,
                                   cl = y_train,
                                   k = neighbors,
                                   prob = TRUE)
             ifelse(y_test == TRUE,
                    attr(kNN_classifier, "prob"),
                    1 - attr(kNN_classifier, "prob"))
           }))

         },
         "rf" = {

           do.call(cbind, lapply(param, function(n) {
             rf_classifier <- ranger(y~.,
                                     data = cbind(X_train, y = y_train, row.names = NULL),
                                     probability = TRUE,
                                     num.trees = n)
             rf_preds <- predict(rf_classifier, cbind(X_test, y = y_test, row.names = NULL))
             rf_preds$predictions[, 1]
           }))

         },
         "naive bayes" = {

           do.call(cbind, lapply(param, function(lapl) {
             naiveBayes_classifier <- naiveBayes(X_train, y_train, laplace = lapl)
             predict(naiveBayes_classifier, X_test, type = "raw")[, 2, drop=FALSE]
           }))

         },
         stop("Wrong model name"))
}
