#' Computes feature selection on n-gram matrices
#' Result is an input to the benchmark.
#' @param paths
#' @param setup
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

      m <- ngram_matrices[[1]]

      for (i in 2:length(ngram_matrices)) {
        m <- rbind_ngram_matrices(m, ngram_matrices[[i]])
      }

      filter_ngrams(m,setup[["method"]])

    })





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

