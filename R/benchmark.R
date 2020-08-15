#' Computes feature selection on n-gram matrices
#' Result is an input to the benchmark.
#' @param paths
#' @param setup
#' @export
create_benchmark_data <- function(paths, setup) {

  results <- setup
  details <- list()

  for (path in paths) {

    m <- read_ngram_matrix(path)

    posNgrams <- positive_ngrams(m)

    res <- filter_ngrams(m, setup[["method"]])

    results <- c(results, res)

  }

  results
}

#' Function summarizes results for a given feature selection method
#'
benchmark_summary <- function() {

}

