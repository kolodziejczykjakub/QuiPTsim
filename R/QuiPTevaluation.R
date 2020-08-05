#' function reads results of `create_simulation_data` function
#' attribute contains simulation details
#' @param path directory containing result files
#' @param title simulation title
#' @importFrom utils read.csv
#' @export
read_simulation_data <- function(path, title) {

  results <- read.csv(paste0(path, title, ".csv"))

  attr(results, "details") <- lapply(results[c("replication", "n_seq", "l_seq", "n_motifs",
                                               "n", "d","seqProbs", "motifProbs")], unique)

  results
}

#' subsets simple triplet matrix
#' @param filePath path of RDS file
#' @param n number of sequences to be subsetted
#' @param fraction fraction of positive sequences
#' @export
subset_matrix <- function(filePath, n, fraction) {
  matrix <- readRDS(filePath)
  target <- attr(matrix, "target")

  n_pos <- round(n * fraction)
  n_neg <- n - n_pos

  pos_indices <- sample(which(target), size = n_pos)
  neg_indices <- sample(which(!target), size = n_neg)

  sampled_matrix <- matrix[c(pos_indices, neg_indices), ]
  attr(sampled_matrix, "sequences") <- attr(matrix, "sequences")[c(pos_indices, neg_indices), ]
  attr(sampled_matrix, "motifs") <- attr(matrix, "motifs")[pos_indices]
  attr(sampled_matrix, "masks") <- attr(matrix, "masks")[pos_indices]
  attr(sampled_matrix, "target") <- c(rep(TRUE, n_pos), rep(FALSE, n_neg))

  sampled_matrix
}

#' function reads output files of `create_simulation_data`
#' @seealso [create_simulation_data()] [subset_matrix()]
#' @param filePath path of RDS file
#' @param n number of sequences to be subsetted (optional)
#' @param fraction fraction of positive sequences (optional)
#' @export

read_ngram_matrix <- function(filePath, n = NULL, fraction = 0.5) {
  if (is.null(n)) {
    matrix <- readRDS(filePath)
  } else {
    matrix <- subset_matrix(filePath, n, fraction)
  }
  matrix
}

#' test QuiPT
#' @param ngram_matrix matrix of n-gram occurences (computed in [create_simulation_data()])
#' @param adjustments p-value adjustments for multiple testing
#' @param thresholds p-value thresholds
#' @return data frame with results of QuiPT testing
#' @importFrom stats p.adjust
#' @importFrom biogram decode_ngrams
#' @importFrom biogram code_ngrams
#' @export

QuiPT_summary <- function(ngram_matrix,
                          adjustments = c("holm", "BH"),
                          thresholds = c(0.05, 0.01)) {

  # create list of unique motifs
  unique_motifs <- unique(unlist(attr(ngram_matrix, "motifs"), recursive = FALSE))

  res <- data.frame(biogram::test_features(target = attr(ngram_matrix, "target"),
                                              features = ngram_matrix))

  res[adjustments] <- lapply(adjustments, function(p.adj) p.adjust(res[["p.value"]], p.adj))
  res[paste0("pval", thresholds)] <- lapply(thresholds, function(th) res[["p.value"]] < th)

  # iterate over motifs to create vector of occurences
  motifOcc <- lapply(unique_motifs, function(motif) {

    motifOcc <- list()

    # ngram contains motif
    noi <- grepl(gsub("_", ".", paste0(motif, collapse = "")), decode_ngrams(res[["ngram"]]))
    motifOcc[["contains.motif"]] <- noi

    # ngram is a part of a motif
    noi2 <- sapply(res[["ngram"]], function(single_ngram)
      grepl(gsub("_", ".", paste0(decode_ngrams(single_ngram), collapse = "")), paste0(motif, collapse = "")))
    motifOcc[["motif.part"]] <- noi2

    # exact motif
    noi3 <- (res[["ngram"]] == code_ngrams(paste0(motif, collapse = "")))
    motifOcc[["motif"]] <- noi3

    motifOcc
  })

  logicalSumOfOccurences <- motifOcc[[1]]
  if (length(motifOcc) > 1) {
    for (i in 2:length(motifOcc)) {
      logicalSumOfOccurences <- Map("|", logicalSumOfOccurences, motifOcc[[i]])
    }
  }

  res[names(logicalSumOfOccurences)] <- logicalSumOfOccurences

  res
}



