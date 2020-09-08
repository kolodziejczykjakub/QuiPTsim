#' Computes feature selection on n-gram matrices
#' Result is an input to the benchmark.
#' @param paths list of n-gram matrices' paths
#' @param setup list of experiment details
#' @export

create_benchmark_data <- function(paths, setup) {

  if (!("shuffle_matrices" %in% names(setup))) {

    res <- pblapply(paths, function(path) {


      if (!("fraction" %in% names(setup) & "n" %in% names(setup))) {
        m <- read_ngram_matrix(path)
        } else {
        m <- read_ngram_matrix(path, n = setup[["n"]], fraction = setup[["fraction"]])
        }

      computation.time <- system.time(res <- filter_ngrams(m, setup[["method"]]))[1]
      list(time = computation.time,
           results = res)
    })
  } else {

    listOfPaths <- lapply(1:length(paths), function(x) sample(paths, setup[["shuffle_matrices"]]))

    res <- pblapply(listOfPaths, function(pathsRep) {

      ngram_matrices <- lapply(pathsRep, function(x)
        read_ngram_matrix(x,
                          n = round(setup[["n"]] / setup[["shuffle_matrices"]], 0),
                          fraction = setup[["fraction"]]))

      m <- rbind_ngram_matrices(ngram_matrices)

      computation.time <- system.time(res <- filter_ngrams(m, setup[["method"]]))[1]
      list(time = computation.time,
           results = res)

    })
  }

  results <- lapply(res, function(x) x[["results"]])
  attr(results, "time") <- unlist(lapply(res, function(x) x[["time"]]))
  attr(results, "setup") <- setup

  results
}

#' Function summarizes results for a given feature selection method
#' @inheritParams calculate_score
#' @export

benchmark_summary <- function(scores, setup) {

  evaluation_metrics <- calculate_score(scores, setup)

  attr(evaluation_metrics, "setup") <- setup
  attr(evaluation_metrics, "time") <- attr(scores, "time")

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

