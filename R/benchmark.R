#' Computes feature selection on n-gram matrices
#' Result is an input to the benchmark.
#' @param paths
#' @param setup
#' @export
create_benchmark_data <- function(paths, setup) {

  if (!("shuffle_matrices" %in% names(setup))) {

    results <- pblapply(paths, function(path) {

      m <- read_ngram_matrix(path)
      filter_ngrams(m,setup[["method"]])
    })

  } else {
    listOfPaths <- lapply(1:length(paths), function(x) sample(paths, setup[["shuffle_matrics"]]))
    # TODO: no idea, how this one should be implemented.

  }


  attr(results, "setup") <- setup

  results
}

#' Function summarizes results for a given feature selection method
#'
benchmark_summary <- function(scores, setup) {

  evaluation_metrics <- calculate_score(scores, setup)

  attr(evaluation_metrics, "setup") <- setup

  evaluation_metrics
}

#'
#'

QuiPTsimBenchmark <- function() {

  # create benchmark data

  # benchmark summary
}

