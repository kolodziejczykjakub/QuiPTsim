#' Computes feature selection on n-gram matrices
#' Result is an input to the benchmark.
#' @param paths
#' @param setup
#' @export
create_benchmark_data <- function(paths, setup) {

  results <- list()
  details <- list()

  for (path in paths) {

    m <- read_ngram_matrix(path)

    res <- filter_ngrams(m, setup[["method"]])

    results <- c(results, list(res))

  }

  attr(results, "setup") <- setup

  results
}

#' Function summarizes results for a given feature selection method
#'
benchmark_summary <- function(scores, setup) {

  evaluation_metrics <- lapply(scores, function(sc) {


    y_true <- sc[["positive.ngram"]]
    y_pred <- calculate_score(sc, setup)

    list(sensitivity = sensitivity(y_true, y_pred),
         specificity = specificity(y_true, y_pred),
         F1score = F1score(y_true, y_pred),
         precision = precision(y_true, y_pred),
         recall = recall(y_true, y_pred))

  })

  evaluation_metrics
}

