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
evaluate_selected_kmers <- function(df, validation_scheme) {

  folds <- createFolds(y = df[["y"]],
                       k = validation_scheme[["folds"]])

  results <- lapply(1:length(folds), function(i) {

    df_train <- df[-folds[[i]], ]
    df_test <- df[folds[[i]], ]

    data.frame(evaluate_models(df_train, df_test, validation_scheme),
               fold = i)

  })

  results <- data.frame(do.call(rbind, results))
  rownames(results) <- NULL

  results
}

#' Function trains and evaluates model on selected k-mers
#' @param df_train
#' @export
#' @importFrom glmnet glmnet
#' @importFrom ranger ranger
#' @importFrom e1071 naiveBayes
#' @importFrom class knn
#' @importFrom pROC auc
evaluate_models <- function(df_train, df_test, validation_scheme) {

  # TODO: works with only one hyperparameter
  models_details <- validation_scheme[["models_details"]]

  X_train <- df_train[, !names(df_train) %in% "y", drop=FALSE]
  y_train <- df_train[["y"]]

  X_test <- df_test[, !names(df_test) %in% "y", drop=FALSE]
  y_test <- df_test[["y"]]

  models_probs <- lapply(models_details, function(m) {
    build_model(X_train, y_train, X_test, y_test, m[["param_value"]], m[["model"]])
  })


  results <- lapply(1:length(models_details), function(i) {

    # in case of not specified lambdas in LASSO lm
    if (models_details[[i]][["model"]] == "lm" &
        is.null(models_details[[i]][["param_value"]])) {
      models_details[[i]][["param_name"]] <- "lambda"
      models_details[[i]][["param_value"]] <- attr(models_probs[[i]], "lambda")
    }

    res <- apply(models_probs[[i]], 2, function(y_pred) {
      c(compute_metrics(y_test, y_pred > 0.5),
        auc = suppressMessages(auc(as.numeric(y_test), y_pred)))
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
#' @importFrom kknn kknn
build_model <- function(X_train, y_train, X_test, y_test, param, method) {

  switch(method,
         "lm" = {

           # glmnet performs a check
           # this is a workaround
           if (ncol(X_train) < 2){
             X_train <- cbind(X_train, rep(1, nrow(X_train)))
             X_test <- cbind(X_test, rep(1, nrow(X_test)))
           }

           # if all variables have zero variance, random guess instead of LM
           if (sum(apply(X_train, 2, var)) == 0){
             
             ans <- sample(0:1, size = nrow(X_test), prob = table(y_train) / length(y_train), replace = TRUE)
             ans <- matrix(ans)
             attr(ans, "lambda") <- -999
             
           } else {
             # if param is not specified, default lambda parameters are computed
             if (is.null(param)) {
               logReg <- glmnet(x = as.matrix(X_train),
                                y = as.numeric(y_train),
                                family = "binomial")
             } else {
               logReg <- glmnet(x = as.matrix(X_train),
                                y = as.numeric(y_train),
                                family = "binomial",
                                lambda = param)
             }
             
             ans <- predict(logReg, as.matrix(X_test), type = "response")
             attr(ans, "lambda") <- logReg[["lambda"]]
           }
           
           ans
         },
         "knn" = {

           do.call(cbind, lapply(param, function(neighbors) {

             kknn_model <- kknn(y~.,
                                data.frame(X_train, y = as.numeric(y_train)),
                                data.frame(X_test),
                                k = neighbors,
                                kernel = "rectangular")

             kknn_model$fitted.values
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
#' @importFrom changepoint cpt.var
#' @export
kmers_for_nonranking_methods <- function(filtering_results, method, thresholds) {

  switch(method,
         "QuiPT" =,
         "Chi-squared"= unlist(lapply(thresholds, function(th)
           sum(filtering_results[["score"]]< th))),
         "gainratio" = ,
         "infogain" = ,
         "symuncert" = cpt.var(sort(filtering_results[["score"]], decreasing = TRUE),
                               class = F,
                               penalty = "BIC")[1],
         "FCBF"= sum(filtering_results[["score"]]),
         stop("Unkown filter!"))
}

#' function creates and evaluates filtering rankings
#' @param paths list of paths containing n-gram matrices
#' @param output_prefix output files directory
#' @param feature_selection_method filter (e.g. QuiPT)
#' @param n number of total sequences
#' @param fraction fraction of positive examples
#' @param validation_scheme list with ranking details
#' @export
filter_rankings <- function(paths, output_prefix, feature_selection_method, n, fraction, validation_scheme) {

  output_paths <- lapply(1:length(paths), function(i)
    paste0(output_prefix, "_",feature_selection_method, "_", i, ".Rds"))

  for (i in 1:length(paths)) {

    m <- read_ngram_matrix(paths[[i]], n, fraction)

    filtering_results <- filter_ngrams(m, feature_selection_method = feature_selection_method)
    results <- evaluate_filtering_results(m,
                                          filtering_results,
                                          list(method = feature_selection_method),
                                          validation_scheme)

    toSave <- list(filtering_results = filtering_results, results = results)

    saveRDS(toSave,
            output_paths[[i]])
    message(paste0("File ", i, " saved in directory: ", output_paths[[i]]))
  }

  output_paths
}

#' function that evaluates nonranking methods
#' @param paths list of paths containing n-gram matrices
#' @param output_prefix output files directory
#' @param feature_selection_method filter (e.g. QuiPT)
#' @param n number of total sequences
#' @param fraction fraction of positive examples
#' @param validation_scheme list with ranking details
#' @param thresholds p-value thresholds for statistical tests
#' @export
filter_nonrankings <- function(paths, output_prefix, feature_selection_method, n, fraction, validation_scheme,
                               thresholds = NULL) {

  output_paths <- lapply(1:length(paths), function(i)
    paste0(output_prefix, "_",feature_selection_method, "_nonranking_", i, ".Rds"))

  for (i in 1:length(paths)) {

    m <- read_ngram_matrix(paths[[i]], n, fraction)

    filtering_results <- filter_ngrams(m, feature_selection_method = feature_selection_method)
    n_kmers <- kmers_for_nonranking_methods(filtering_results,
                                            feature_selection_method,
                                            thresholds)

    print(paste0("Method: ", feature_selection_method, ", path: ", paths[[i]], ", number of k-mers:", n_kmers))

    # if threshold if further, trim to best 4096 k-mers
    n_kmers = min(n_kmers, 4096)

    validation_scheme[["n_kmers"]] <- n_kmers

    results <- evaluate_filtering_results(m,
                                          filtering_results,
                                          list(method = feature_selection_method),
                                          validation_scheme)

    toSave <- list(filtering_results = filtering_results, results = results)

    saveRDS(toSave,
            output_paths[[i]])
    message(paste0("File ", i, " saved in directory: ", output_paths[[i]]))

  }
  output_paths
}

#' function creates and evaluates filtering rankings
#' @param paths list of paths containing n-gram matrices
#' @param num_reps number of repetitions of an experiment
#' @param num_matrices_to_rbind number of matrices to bind in order to perform exp III correctly
#' @param output_prefix output files directory
#' @param feature_selection_method filter (e.g. QuiPT)
#' @param n number of total sequences
#' @param fraction fraction of positive examples
#' @param validation_scheme list with ranking details
#' @export
filter_rankings_exp3 <- function(paths, num_reps, num_matrices_to_rbind, output_prefix, feature_selection_method, n, fraction, validation_scheme) {
  
  output_paths <- lapply(1:length(num_reps), function(i)
    paste0(output_prefix, "_",feature_selection_method, "_", i, ".Rds"))
  
  for (i in 1:length(num_reps)) {
    
    submatrix_nrow <- as.integer(n / num_matrices_to_rbind)
    m <- do.call(rbind_ngram_matrices, lapply(sample(paths, num_matrices_to_rbind), function(x)
      read_ngram_matrix(x, submatrix_nrow, fraction)))
    
    filtering_results <- filter_ngrams(m, feature_selection_method = feature_selection_method)
    results <- evaluate_filtering_results(m,
                                          filtering_results,
                                          list(method = feature_selection_method),
                                          validation_scheme)
    
    toSave <- list(filtering_results = filtering_results, results = results)
    
    saveRDS(toSave,
            output_paths[[i]])
    message(paste0("File ", i, " saved in directory: ", output_paths[[i]]))
  }
  
  output_paths
}

#' function that evaluates nonranking methods
#' @param paths list of paths containing n-gram matrices
#' @param num_reps number of repetitions of an experiment
#' @param num_matrices_to_rbind number of matrices to bind in order to perform exp III correctly
#' @param output_prefix output files directory
#' @param feature_selection_method filter (e.g. QuiPT)
#' @param n number of total sequences
#' @param fraction fraction of positive examples
#' @param validation_scheme list with ranking details
#' @param thresholds p-value thresholds for statistical tests
#' @export
filter_nonrankings_exp3 <- function(paths, num_reps, num_matrices_to_rbind, output_prefix, feature_selection_method, n, fraction, validation_scheme,
                               thresholds = NULL) {
  
  output_paths <- lapply(1:length(num_reps), function(i)
    paste0(output_prefix, "_",feature_selection_method, "_nonranking_", i, ".Rds"))
  
  for (i in 1:length(num_reps)) {
    
    submatrix_nrow <- as.integer(n / num_matrices_to_rbind)
    m <- do.call(rbind_ngram_matrices, lapply(sample(paths, num_matrices_to_rbind), function(x)
      read_ngram_matrix(x, submatrix_nrow, fraction)))
    
    filtering_results <- filter_ngrams(m, feature_selection_method = feature_selection_method)
    n_kmers <- kmers_for_nonranking_methods(filtering_results,
                                            feature_selection_method,
                                            thresholds)
    
    print(paste0("Method: ", feature_selection_method, ", path: ", paths[[i]], ", number of k-mers:", n_kmers))
    
    # if threshold if further, trim to best 4096 k-mers
    n_kmers <- ifelse(n_kmers > 4096, 4096, n_kmers)
    
    validation_scheme[["n_kmers"]] <- n_kmers
    
    results <- evaluate_filtering_results(m,
                                          filtering_results,
                                          list(method = feature_selection_method),
                                          validation_scheme)
    
    toSave <- list(filtering_results = filtering_results, results = results)
    
    saveRDS(toSave,
            output_paths[[i]])
    message(paste0("File ", i, " saved in directory: ", output_paths[[i]]))
    
  }
  output_paths
}

