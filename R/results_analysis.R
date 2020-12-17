#' function aggregates computation times 
#' @param paths list of paths that contain filtering results
#' @return computation times
#' @export
collect_filtering_times <- function(paths) {
  lapply(paths, function(path) {
    x <- readRDS(path)
    attr(x$filtering_results, "time")[1]
  })
}

#' parsing filtering results
#' @inheritParams collect_filtering_results
#' @return aggregated results
#' @export
#' @importFrom pbapply pblapply 
#' @importFrom dplyr group_by summarise
parse_results <- function(paths) {
  
  x <- pblapply(paths, function(path) {
    
    results <- readRDS(path)[["results"]]
    
    # clearing row containing dummy 1-mer
    results_lm <- results[results$model == "lm" & results$n_kmers != 1, ]
    results_wo_lm <- results[results$model != "lm" & results$n_kmers != 1, ]
    
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
    
    data.frame(path = path, rbind(agg_results_wo_lm, agg_results_lm))
  })
  
  do.call(rbind, x)
}