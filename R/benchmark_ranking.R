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

    evaluate_models(df_train, df_test, validation_scheme)

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
evaluate_models <- function(df_train, df_test, validation_scheme) {

  # TODO: works with only one hyperparameter
  models_details <- validation_scheme[["models_details"]]

  X_train <- df_train[, !names(df_train) %in% "y" ]
  y_train <- df_train[["y"]]

  X_test <- df_test[, !names(df_test) %in% "y" ]
  y_test <- df_test[["y"]]

  models_probs <- lapply(models_details, function(m)
    build_model(X_train, y_train, X_test, y_test, m[["param_value"]], m[["model"]]))

  results <- lapply(1:length(models_details), function(i) {

    res <- apply(models_probs[[i]], 2, function(y_pred) {
      c(compute_metrics(y_test, y_pred > 0.5),
        auc = auc(as.numeric(y_test), y_pred))
    })

    data.frame(model = models_details[[i]][["model"]],
               param = models_details[[i]][["param_name"]],
               value = models_details[[i]][["param_value"]],
               do.call(rbind, res),
               row.names = NULL)
  })

  do.call(rbind, results)
}

#' function builds a model and predicts out of fold probabilites
#' @param X_train
#' @param y_train
#' @param X_test
#' @param y_test
#' @param param vector of parameters (lm : lambda,
#'  knn: neighbors, naive bayes: laplace, rf: num.trees)
#' @param method model (one of: lm, knn, naive bayes, rf)
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

#' calculates number of k-mers to select for non-ranking methods
#' @param filtering_results data.frame containing ranking of k-mers
#' @param method feature filtering method
#' @thresholds optional parameter for statistical tests
#' @export
kmers_for_nonranking_methods <- function(filtering_results, method, thresholds) {

  switch(method,
         "QuiPT" =,
         "Chi-squared"= unlist(lapply(thresholds, function(th)
           sum(filtering_results[["p.value"]]< th))),
         "gainratio" = ,
         "infogain" = ,
         "symuncert" = cpt.var(sort(filtering_results, decreasing = TRUE),
                               class = F,
                               penalty = "BIC")[1],
         "FCBF"= sum(filtering_results[["score"]]),
         stop("Unkown filter!"))
}
