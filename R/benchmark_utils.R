#' Function returns ector indicating n-grams considered positive
#' @description Positive n-gram is either n-gram containing motif (not longer than n + (n+1) // 4)
#' or n-gram that is a part of motif
#' @param motifs
#' @param ngram_matrix
#' @export

positive_ngrams <- function(motifs, ngram_matrix) {

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


