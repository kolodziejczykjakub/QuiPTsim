#' Function returns ector indicating n-grams considered positive
#' @description Positive n-gram is either n-gram containing motif (not longer than n + (n+1) // 4)
#' or n-gram that is a part of motif
#' @param motifs
#' @param ngram_matrix
#' @export

positive_ngrams <- function(ngram_matrix) {

  if (!(class(ngram_matrix) == "simple_triplet_matrix")) {
    stop("Input is not a simple_triplet_matrix class instance!")
  }

  if (!("motifs" %in% names(attributes(ngram_matrix)))) {
    stop("Input does not have motifs attribute!")
  }

  motifs <- unique(unlist(attr(ngram_matrix, "motifs"), recursive = FALSE))
  ngrams <- colnames(ngram_matrix)

  pos <- lapply(motifs, function(motif) {

    pos <- list()

    # n-grams containing full motif
    noi <- grepl(paste0(motif, collapse = ""), decode_ngrams((ngrams)))

    # n-gram length limit
    ngram_lengths <- unlist(lapply(ngrams, function(x) nchar(decode_ngrams(x))))
    ngrams_containing_motif <- (noi & (ngram_lengths < (length(motif)+1) %/% 4 + length(motif)))

    # n-grams being a part of a motif
    ngrams_motif_part <- sapply(ngrams, function(single_ngram)
      grepl(paste0(decode_ngrams(single_ngram), collapse = ""), paste0(motif, collapse = "")))

    # exact motif
    exact_motif <- (ngrams == code_ngrams(paste0(motif, collapse = "")))

    pos[["motif"]] <- exact_motif
    pos[["positive.ngram"]] <- (ngrams_containing_motif | ngrams_motif_part | exact_motif)

    pos
    })

  logicalSumOfOccurences <- pos[[1]]
  if (length(pos) > 1) {
    for (i in 2:length(pos)) {
      logicalSumOfOccurences <- Map("|", logicalSumOfOccurences, pos[[i]])
    }
  }

  logicalSumOfOccurences

}

#' Wrapper for variuous feature selection method evaluated in a benchmark
#' @param ngram_matrix matrix of n-gram occurences
#' @param feature_selection_method feature selection method name (QuiPT, ...)

filter_ngrams <- function(ngram_matrix, feature_selection_method) {

  if (!(feature_selection_method %in% c("QuiPT", "Boruta"))) {
    stop("Unkown feature selection method!")
  }

  # QuiPT feature selection
  if (feature_selection_method == "QuiPT") {
    out <- data.frame(test_features(target = attr(ngram_matrix, "target"),
                                             features = ngram_matrix,
                                    threshold = 0))
    out[["score"]] <- out[["p.value"]]
  }

  # Boruta algorithm
  if (feature_selection_method == "Boruta") {

    browser()
    Boruta(as.matrix(m), attr(m, "target"))
  }



  # add a column indicating if given ngram is positive
  posNgrams <- positive_ngrams(ngram_matrix)
  out[names(posNgrams)] <- posNgrams

  out
}

#'
#'

calculate_score <- function(scores, setup) {

  method <- setup[["method"]]

  if (!(method %in% c("QuiPT", "Boruta"))) {
    stop("Unkown feature selection method!")
  }

  # QuiPT setup
  if (method == "QuiPT") {

    if (!("pval_thresholds" %in% names(setup))) {
      stop("P-value thresholds have not been set!")
    }

    if (!("pval_adjustments" %in% names(setup))) {
      stop("P-value adjustments have not been set!")
    }

    results <- expand.grid(pval_thresholds = setup[["pval_thresholds"]],
                          pval_adjustments = setup[["pval_adjustments"]])

    for (pval_th in setup[["pval_thresholds"]]) {
      for (pval_adj in setup[["pval_adjustments"]]) {

        out <- list()

        for (i in 1:length(scores)) {

          pvals <- scores[[i]][["score"]]
          y_true <- scores[[i]][["positive.ngram"]]

          if (!(pval_adj == "")) {
            y_pred <- p.adjust(pvals, method = pval_adj)
          } else {
            y_pred <- pvals
          }

          y_pred <- (y_pred < pval_th)

          out[[i]] <- compute_metrics(y_true, y_pred)

        }

        metrics_names <- names(out[[1]])

        aggregated_metrics <- lapply(lapply(metrics_names, function(metric_name)
          lapply(out, function(r)
            r[[metric_name]])), unlist)

        aggregated_out <- c(setNames(sapply(aggregated_metrics, mean),
                                     paste0(metrics_names, "_mean")),
                            setNames(sapply(aggregated_metrics, sd),
                                     paste0(metrics_names, "_std")))

        results[results[["pval_thresholds"]] == pval_th & results[["pval_adjustments"]] == pval_adj,
               names(aggregated_out)] = aggregated_out

      }
    }
  }

  # Boruta setup
  if (method == "Boruta") {
    stop("Not done yet")
  }

  results

}









