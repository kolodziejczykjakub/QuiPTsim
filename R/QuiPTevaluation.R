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

# wrapper for subset_matrix
read_ngram_matrix <- function(filePath, n = NULL, fraction = 0.5) {
  if (is.null(n)) matrix <- readRDS(filePath)
  else matrix <- subset_matrix(filePath, n, fraction)
  matrix
}

#' test QuiPT
#' @param simulated_seqs sequences that will be tested
#' @param n_seq number of sequences in each category
#' @param criterion
#' @param motifs
#' @param adjustment
#'
#'
#' @param motifs list of injected motifs
#' @return data frame with results of QuiPT testing
#' @examples
#' n_seq <- 20
#' len <- 10
#' alphabet <- 1L:4
#' motifs <- generate_motifs(alphabet, 3)
#' simulated_sequences <- generate_sequences(n_seq, len, alph, motifs, 1)
#' test_quipt(simulated_sequences, n_seq, criterion = "ig", motifs = motifs)

test_quipt <- function(simulated_seqs, n_seq, criterion, motifs, adjustment = "BH") {

  res_df <- data.frame(biogram::test_features(target = c(rep(1, n_seq),
                                              rep(0, n_seq)),
                                   features = simulated_seqs))


  data.frame(#n_seq = n_seq,
    #l_seq = l_seq,
    res_df,
    motif = res_df[["ngram"]] %in% biogram::code_ngrams(sapply(motifs,
                                                      paste0,
                                                      collapse = "")),
    p.value.adj = p.adjust(res_df[["p.value"]], adjustment))
}

QuiPT_summary <- function(ngram_matrix) {

  # create list of unique motifs
  unique_motifs <- unique(unlist(attr(dat, "motifs"), recursive = FALSE))

  # iterate over motifs
  lapply(unique_motifs, function(motif) {

    res_df <- data.frame(biogram::test_features(target = attr(ngram_matrix, "target"),
                                                features = ngram_matrix))

    # TODO: all adjustments at once?
    adjustment <- "BH"

    data.frame(
      res_df,
      motif = res_df[["ngram"]] %in% biogram::code_ngrams(sapply(motifs,
                                                                 paste0,
                                                                 collapse = "")),
      p.value.adj = p.adjust(res_df[["p.value"]], adjustment)) -> res

    # ngram contains motif
    noi <- grepl(decode_ngrams(single_motif), decode_ngrams(res[["ngram"]]))
    res[["contains.motif"]] <- noi

    # ngram is a part of a motif
    sapply(res[["ngram"]], function(single_ngram)
      grepl(decode_ngrams(single_ngram), decode_ngrams(single_motif))) -> noi2
    res[["motif.part"]] <- noi2

    res
  })
}



