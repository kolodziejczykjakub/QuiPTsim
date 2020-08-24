#' Computes feature selection on n-gram matrices
#' Result is an input to the benchmark.
#' @param paths list of n-gram matrices' paths
#' @param setup list of experiment details
#' @export

create_benchmark_data <- function(paths, setup) {

  if (!("shuffle_matrices" %in% names(setup))) {

    results <- pblapply(paths, function(path) {


      if (!("fraction" %in% names(setup) & "n" %in% names(setup))) {
        m <- read_ngram_matrix(path)
      } else {
        m <- tryCatch({read_ngram_matrix(path, n = setup[["n"]], fraction = setup[["fraction"]])},
                      error = function(c) {
                        message("Either number of sequences or fraction have not been set")
                        message("Original message:")
                        message(c)
                        })

      }

      filter_ngrams(m,setup[["method"]])
    })

  } else {

    listOfPaths <- lapply(1:length(paths), function(x) sample(paths, setup[["shuffle_matrices"]]))

    results <- pblapply(listOfPaths, function(pathsRep) {

      ngram_matrices <- lapply(pathsRep, function(x)
        read_ngram_matrix(x,
                          n = round(setup[["n"]] / setup[["shuffle_matrices"]], 0),
                          fraction = setup[["fraction"]]))

      m <- rbind_ngram_matrices(ngram_matrices)

      filter_ngrams(m,setup[["method"]])

    })

}


  attr(results, "setup") <- setup

  results
}

#' Function summarizes results for a given feature selection method
#' @inheritParams calculate_score
#' @export

benchmark_summary <- function(scores, setup) {

  evaluation_metrics <- calculate_score(scores, setup)

  attr(evaluation_metrics, "setup") <- setup

  evaluation_metrics
}


#' Function wraps up QuiPT evaluation auxiliary functions
#' @param paths list of paths of n-gram matrices
#' @param setup list of experiment details
#' @seealso `create_benchmark_data` `benchmark_summary`
#' @export

QuiPTsimBenchmark <- function(paths, setup) {

  # create benchmark data
  benchmark_data <- create_benchmark_data(paths, setup)

  # benchmark summary
  benchmark_summary(benchmark_data, attr(benchmark_data, "setup"))
}

