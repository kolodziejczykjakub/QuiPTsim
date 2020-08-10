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
#' @param mc markovchain object to calculate n-gram probability
#' @param adjustments p-value adjustments for multiple testing
#' @param thresholds p-value thresholds
#' @return data frame with results of QuiPT testing
#' @importFrom stats p.adjust
#' @importFrom biogram decode_ngrams
#' @importFrom biogram code_ngrams
#' @export

QuiPT_summary <- function(ngram_matrix,
                          mc = NULL,
                          adjustments = c("holm", "BH"),
                          thresholds = c(0.05, 0.01)) {

  # create list of unique motifs
  unique_motifs <- unique(unlist(attr(ngram_matrix, "motifs"), recursive = FALSE))

  res <- data.frame(biogram::test_features(target = attr(ngram_matrix, "target"),
                                              features = ngram_matrix))

  res[adjustments] <- lapply(adjustments, function(p.adj) p.adjust(res[["p.value"]], p.adj))
  res[paste0("pval", thresholds)] <- lapply(thresholds, function(th) res[["p.value"]] < th)

  if (!is.null(mc)) {
    res[["prob"]] <- pblapply(strsplit(decode_ngrams(res[["ngram"]]), ""),
                              function(ngram) calculate_ngram_prob(mc, ngram))
  }

  # iterate over motifs to create vector of occurences
  motifOcc <- lapply(unique_motifs, function(motif) {

    motifOcc <- list()

    # ngram contains motif
    noi <- grepl(paste0(motif, collapse = ""), decode_ngrams(res[["ngram"]]))
    motifOcc[["contains.motif"]] <- noi
    ngram_lengths <- unlist(lapply(res[["ngram"]], function(x) nchar(decode_ngrams(x))))
    maxLen <- (length(motif)+1) %/% 4 + length(motif)
    motifOcc[["positive.ngram"]] <- ((ngram_lengths < maxLen) & noi)
    # ngram is a part of a motif
    noi2 <- sapply(res[["ngram"]], function(single_ngram)
      grepl(paste0(decode_ngrams(single_ngram), collapse = ""), paste0(motif, collapse = "")))
    motifOcc[["motif.part"]] <- noi2

    # exact motif
    noi3 <- (res[["ngram"]] == code_ngrams(paste0(motif, collapse = "")))
    motifOcc[["motif"]] <- noi3


    positiveMotifs <- lapply(positive_motifs(motif), function(x) paste0(x, collapse=""))
    noi4 <- sapply(res[["ngram"]], function(single_ngram)
      decode_ngrams(single_ngram) %in% positiveMotifs)
    motifOcc[["positive"]] <- noi4

    #TODO: define final positive ngram group

    motifOcc
  })

  #TODO: more sensible implementation
  logicalSumOfOccurences <- motifOcc[[1]]
  if (length(motifOcc) > 1) {
    for (i in 2:length(motifOcc)) {
      logicalSumOfOccurences <- Map("|", logicalSumOfOccurences, motifOcc[[i]])
    }
  }

  res[names(logicalSumOfOccurences)] <- logicalSumOfOccurences

  res
}

#' function combines two ngram matrices
#' @param m1
#' @param m2
#' @return combined ngram matrix
#' @export

rbind_ngram_matrices <- function(m1, m2) {

  m1_unique_colnames <- setdiff(colnames(m1), colnames(m2))
  m2_unique_colnames <- setdiff(colnames(m2), colnames(m1))

  m1_extended <- cbind(m1, matrix(0, nrow = nrow(m1), ncol = length(m2_unique_colnames)))
  m2_extended <- cbind(m2, matrix(0, nrow = nrow(m2), ncol = length(m1_unique_colnames)))


  colnames(m1_extended) <- c(colnames(m1), m2_unique_colnames)
  colnames(m2_extended) <- c(colnames(m2), m1_unique_colnames)

  m_extended <- rbind(m1_extended, m2_extended[, colnames(m1_extended)])

  attr(m_extended, "sequences") <- rbind(attr(m1, "sequences"), attr(m2, "sequences"))
  attr(m_extended, "motifs") <- c(attr(m1, "motifs"), attr(m2, "motifs"))
  attr(m_extended, "masks") <- c(attr(m1, "masks"), attr(m2, "masks"))
  attr(m_extended, "target") <- c(attr(m1, "target"), attr(m2, "target"))

  m_extended
}

