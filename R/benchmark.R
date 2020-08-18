#' Computes feature selection on n-gram matrices
#' Result is an input to the benchmark.
#' @param paths
#' @param setup
#' @export
create_benchmark_data <- function(paths, setup) {

  results <- list()
  details <- list()

  totalIterations = length(paths)
  pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)",
                         total = totalIterations)
  for (i in 1:totalIterations) {

    path <- paths[i]

    m <- read_ngram_matrix(path)
    res <- filter_ngrams(m, setup[["method"]])
    results <- c(results, list(res))

    pb$tick(1)
  }

  attr(results, "setup") <- setup

  results
}

#' Function summarizes results for a given feature selection method
#'
benchmark_summary <- function(scores, setup) {

    calculate_score(scores, setup)

#
#   metrics_names <- names(evaluation_metrics[[1]])
#
#   aggregated_metrics <- lapply(lapply(metrics_names, function(metric_name)
#     lapply(evaluation_metrics, function(results)
#       results[[metric_name]])), unlist)
#
#   data.frame(list(mean = sapply(aggregated_metrics, mean),
#                   standard.dev = sapply(aggregated_metrics, sd)),
#              row.names = metrics_names)

  evaluation_metrics
}

#'
#'

QuiPTsimBenchmark <- function() {

  # create benchmark data

  # benchmark summary
}

