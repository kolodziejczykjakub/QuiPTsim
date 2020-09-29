#' Function returns ector indicating n-grams considered positive
#' @description
#' Positive n-gram is either n-gram containing motif (not longer than n + (n+1) // 4)
#' or n-gram that is a part of motif
#' @param ngram_matrix matrix of n-gram occurences
#' @importFrom biogram decode_ngrams code_ngrams
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

#' Wrapper for feature selection methods
#' @description
#' Wrapper for variuous feature selection method evaluated in a QuiPTsim benchmark
#' @param ngram_matrix matrix of n-gram occurences
#' @param feature_selection_method feature selection method name (QuiPT, ...)
#' @importFrom biogram test_features
#' @importFrom stats setNames
#' @importFrom stats sd
#' @importFrom FCBF fcbf
#' @importFrom FSelectorRcpp information_gain
#' @importFrom praznik MIM MRMR JMI JMIM DISR NJMIM CMIM
#' @export

filter_ngrams <- function(ngram_matrix, feature_selection_method) {

  # Feature selection methods from praznik package
  praznik_filters <- c("MIM", "MRMR", "JMI", "JMIM", "DISR", "NJMIM", "CMIM", "CMI")
  FSelectorRcpp_measures <- c("infogain", "gainratio", "symuncert")

  if (!(feature_selection_method %in% c("QuiPT",
                                        "FCBF",
                                        "Chi-squared",
                                        FSelectorRcpp_measures,
                                        praznik_filters))) {
    stop("Unkown feature selection method!")
  }

  if (feature_selection_method %in% praznik_filters) {

    x <- data.frame(as.matrix(ngram_matrix))
    y <- attr(ngram_matrix, "target")

    res <- get(feature_selection_method)(X = x,
                                         Y = y,
                                         k = round(ncol(x) * 0.05, digits = 0))

    filtered_ngrams <- rep(length(colnames(x)), length(colnames(x)))
    filtered_ngrams[res$selection] <- 1:(round(ncol(x) * 0.05, digits = 0))
    out <- data.frame(rank = filtered_ngrams)
    out[["score"]] <- out[["rank"]] > 0
  }

  if (feature_selection_method %in% FSelectorRcpp_measures) {

    x <- data.frame(as.matrix(ngram_matrix))
    y <- attr(ngram_matrix, "target")

    score <- information_gain(x = x, y = y, discIntegers = FALSE,
                           type = feature_selection_method)

    if (feature_selection_method == "gainratio") {
      score[["importance"]][is.na(score[["importance"]])] <- 0
    }

    out <- data.frame(ngram = score$attributes,
                      score = score$importance,
                      rank = order(score$importance, decreasing = TRUE))

  }

  # QuiPT feature selection
  if (feature_selection_method == "QuiPT") {
    out <- data.frame(test_features(target = attr(ngram_matrix, "target"),
                                             features = ngram_matrix,
                                    threshold = -1))
    out[["score"]] <- out[["p.value"]]
    out[["rank"]] <- order(out[["p.value"]])
  }

  # Chi-squared test
  if (feature_selection_method == "Chi-squared") {

    y <- attr(ngram_matrix, "target")

    pvalues <- lapply(1:ncol(ngram_matrix), function(i) {

      x <- as.vector(ngram_matrix[, i])

      if (length(unique(x)) == 1) {
        pval = 1
      } else {
        pval <- chisq.test(x, y, B = 200)$p.value
      }
      pval
      })

    out <- data.frame(score = unlist(pvalues),
                      rank = order(unlist(pvalues)))
  }

  if (feature_selection_method == "FCBF") {

    features <- apply(ngram_matrix, 1, as.factor)
    ngram_names <- rownames(features)
    n_ngrams <- length(ngram_names)

    fcbf_results <- fcbf(features, attr(ngram_matrix, "target"), verbose = FALSE)
    fcbf_id <- fcbf_results[["index"]]

    pred <- logical(n_ngrams)
    pred[fcbf_id] <- TRUE

    SU <- logical(n_ngrams)
    SU[fcbf_id] <- fcbf_results[["SU"]]

    n_features <- length(fcbf_id)
    rank_ <- rep(-1, n_ngrams)
    rank_[fcbf_id] <- 1:n_features
    rank_[-fcbf_id] <- (n_features + 1):n_ngrams

    out <- data.frame(ngram = ngram_names,
                      score = pred,
                      rank = rank_,
                      SU = SU)

  }

  # add a column indicating if given ngram is positive
  posNgrams <- positive_ngrams(ngram_matrix)
  out[names(posNgrams)] <- posNgrams

  out
}

#' Aggregates metrics
#' @export
#'
aggregate_metrics <- function(metrics) {

  metrics_names <- names(metrics[[1]])

  aggregated_metrics <- lapply(lapply(metrics_names, function(metric_name)
    lapply(metrics, function(r)
      r[[metric_name]])), unlist)

  aggregated_metrics <- c(setNames(sapply(aggregated_metrics, mean),
                                   paste0(metrics_names, "_mean")),
                          setNames(sapply(aggregated_metrics, sd),
                                   paste0(metrics_names, "_std")))

  data.frame(as.list(aggregated_metrics))

}

#' Calculates metrics for a given feature selection method
#' @description
#' Function returns metrics average based on scores computed on n-gram matrices
#' @param scores list of filter scores
#' @param setup evaluation details
#' @export

calculate_score <- function(scores, setup) {

  method <- setup[["method"]]

  # Feature selection methods from praznik package
  praznik_filters <- c("MIM", "MRMR", "JMI", "JMIM", "DISR", "NJMIM", "CMIM", "CMI")
  FSelectorRcpp_measures <- c("infogain", "gainratio", "symuncert")


  if (!(method %in% c("QuiPT",
                      "FCBF",
                      "Chi-squared",
                      FSelectorRcpp_measures,
                      praznik_filters))) {
    stop("Unkown feature selection method!")
  }

  if (method %in% FSelectorRcpp_measures) {

    numFeatures <- sapply(setup[["fraction"]], function(frac)
      round(nrow(scores[[1]]) * frac, digits = 0))

    res <- lapply(numFeatures, function(nFeat) {
      metrics <- lapply(scores, function(res) {

        if (nFeat == 0 & 0 %in% setup[["fraction"]]) {

          sorted_crit <- sort(res[[metric]], decreasing = TRUE)
          id <- cpt.mean(sorted_crit, penalty="SIC", method="AMOC",class=FALSE)[1]
          y_pred <- res[["score"]] > sorted_crit[id]

        } else {

          y_pred <- res[["score"]] > sort(res[["score"]],
                                          decreasing = TRUE)[nFeat]
        }

        compute_metrics(res[["positive.ngram"]], y_pred)
      })

      aggregate_metrics(metrics)
    })

    results <- data.frame(num_features = numFeatures,
                          top_fraction = setup[["fraction"]],
                          do.call(rbind, res))
  }

  # FCBF setup
  if (method %in% c("FCBF", praznik_filters)) {
    out <- lapply(scores, function(res)
      compute_metrics(res[["positive.ngram"]], res[["score"]]))

    results <- aggregate_metrics(out)
  }

  # QuiPT setup
  if (method %in% c("QuiPT", "Chi-squared")) {

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

        aggregated_out <- aggregate_metrics(out)

        results[results[["pval_thresholds"]] == pval_th & results[["pval_adjustments"]] == pval_adj,
               names(aggregated_out)] = aggregated_out

      }
    }
  }

  results
}









