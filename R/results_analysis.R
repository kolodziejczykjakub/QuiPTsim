#' function aggregates computation times 
#' @param paths list of paths that contain filtering results
#' @return computation times
#' @importFrom pbapply pblapply
#' @export
collect_filtering_times <- function(paths) {
  ans <- pblapply(paths, function(path) {
    x <- readRDS(path)
    attr(x$filtering_results, "time")[1]
  })
  
  unlist(ans)
}

#' parsing filtering results for ranking methods
#' @inheritParams collect_filtering_results
#' @return aggregated results
#' @export
#' @importFrom pbapply pblapply 
#' @importFrom dplyr group_by summarise
parse_results <- function(paths) {
  
  x <- pblapply(paths, function(path) {
    
    results <- readRDS(path)[["results"]]
    
    # clearing row containing dummy 1-mer
    results_lm <- results[results$model == "lm" & results$model != 1, ]
    results_wo_lm <- results[results$model != "lm" & results$model != 1, ]
    
    agg_results_wo_lm <- summarise(group_by(results_wo_lm, n_kmers, model, param, value),
                                   accuracy = mean(as.numeric(accuracy)),
                                   sensitivity = mean(as.numeric(sensitivity)),
                                   specificity = mean(as.numeric(specificity)),
                                   F1score = mean(as.numeric(F1score)),
                                   precision = mean(as.numeric(precision)),
                                   recall = mean(as.numeric(recall)),
                                   auc = mean(as.numeric(auc)), .groups = 'keep')
    
    # Results for LM are aggregated by maximising metrics fold-wise
    agg_results_lm_1 <- summarise(group_by(results_lm,n_kmers, fold, model),
                                  accuracy = max(as.numeric(accuracy)),
                                  sensitivity = max(as.numeric(sensitivity)),
                                  specificity = max(as.numeric(specificity)),
                                  F1score = max(as.numeric(F1score)),
                                  precision = max(as.numeric(precision)),
                                  recall = max(as.numeric(recall)),
                                  auc = max(as.numeric(auc)), .groups = 'keep')
    
    agg_results_lm <- summarise(group_by(agg_results_lm_1, n_kmers, model),
                                accuracy = mean(as.numeric(accuracy)),
                                sensitivity = mean(as.numeric(sensitivity)),
                                specificity = mean(as.numeric(specificity)),
                                F1score = mean(as.numeric(F1score)),
                                precision = mean(as.numeric(precision)),
                                recall = mean(as.numeric(recall)),
                                auc = mean(as.numeric(auc)), .groups = 'keep')
    
    out = data.frame(path = path, rbind(agg_results_wo_lm, agg_results_lm))
    
    # For QuiPT and Chi tests
    # thresholds are added
    
    if (grepl("nonranking", path)) {
      kmer_levels <- sort(unique(out$n_kmers))
      thresholds <- c(0.01, 0.05)
      out[["threshold"]] <- NA
      if (length(kmer_levels) > 1){
        for (i in 1:length(kmer_levels)){
          out[out[["n_kmers"]] == kmer_levels[i], "threshold"] = thresholds[i]
        }
      }
    }
    
    out
  })
  
  results <- do.call(rbind, x)
  
  # pretty model names
  results[["Model"]] <- unlist(lapply(paste(results[["model"]], results[["param"]], results[["value"]]), function(x)
    switch(x,
           "knn neighbors 1" = "1-NN",
           "knn neighbors 2" = "2-NN",
           "knn neighbors 4" = "4-NN",
           "knn neighbors 8" = "8-NN",
           "knn neighbors 16" = "16-NN",
           "lm NA NA" = "LR LASSO",
           "rf num.trees 500" = "RF (500 trees)",
           "rf num.trees 1000" = "RF (1000 trees)",
           "naive bayes laplace 0" = "Naive Bayes")))
  
  results[["Model"]] <- factor(results[["Model"]], 
                               levels = c(
    "LR LASSO",
    "RF (500 trees)",
    "RF (1000 trees)",
    "1-NN",
    "2-NN",
    "4-NN",
    "8-NN",
    "16-NN",
    "Naive Bayes"
  ))

  results
}

#' aggregates results from `parse_results`
#' @param df output of `parse_results()`
#' @return aggregated results (metrics means and SEs)
#' @export
aggregate_results <- function(df, ranking=TRUE) {
  
  std <- function(x) sd(x)/sqrt(length(x))
  
  if (ranking) {
    summarise(group_by(df, n_kmers, model, param, value, Model),
              Accuracy_mean = mean(as.numeric(accuracy)),
              Sensitivity_mean = mean(as.numeric(sensitivity)),
              Specificity_mean = mean(as.numeric(specificity)),
              F1Score_mean = mean(as.numeric(F1score)),
              Precision_mean = mean(as.numeric(precision)),
              Recall_mean = mean(as.numeric(recall)),
              AUC_mean = mean(as.numeric(auc)),
              
              Accuracy_std = std(as.numeric(accuracy)),
              Sensitivity_std = std(as.numeric(sensitivity)),
              Specificity_std = std(as.numeric(specificity)),
              F1Score_std = std(as.numeric(F1score)),
              Precision_std = std(as.numeric(precision)),
              Recall_std = std(as.numeric(recall)),
              AUC_std = std(as.numeric(auc)),
              .groups = 'keep')
  } else { 
    summarise(group_by(df, model, param, value, Model, threshold),
              
              n_kmers_mean = mean(as.numeric(n_kmers)),
              Accuracy_mean = mean(as.numeric(accuracy)),
              Sensitivity_mean = mean(as.numeric(sensitivity)),
              Specificity_mean = mean(as.numeric(specificity)),
              F1Score_mean = mean(as.numeric(F1score)),
              Precision_mean = mean(as.numeric(precision)),
              Recall_mean = mean(as.numeric(recall)),
              AUC_mean = mean(as.numeric(auc)),
              
              n_kmers_std = std(as.numeric(n_kmers)),
              Accuracy_std = std(as.numeric(accuracy)),
              Sensitivity_std = std(as.numeric(sensitivity)),
              Specificity_std = std(as.numeric(specificity)),
              F1Score_std = std(as.numeric(F1score)),
              Precision_std = std(as.numeric(precision)),
              Recall_std = std(as.numeric(recall)),
              AUC_std = std(as.numeric(auc)),
              .groups = 'keep')
  }

}
