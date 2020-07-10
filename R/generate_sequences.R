library(biogram)
library(dplyr)
library(pbapply)

#' function generates motif from an alphabet
#' TODO: motif length control
#' @param alphabet elements used to build motif
#' @return motif based on given alphabet
#' @export
#' @examples
#' generate_single_motif(1:4)

generate_single_motif <- function(alphabet) {
  ns <- c(rep(2, 6), rep(3, 9))
  ds <- c(as.list(0L:5),
          split(expand.grid(0L:2, 0L:2), 1L:nrow(expand.grid(0L:2, 0L:2))))

  ngram_id <- sample(1L:length(ns), 1)

  n <- ns[ngram_id]
  d <- ds[[ngram_id]]

  c(unlist(lapply(1L:length(d), function(d_id) {
    c(sample(alphabet, 1), rep("_", d[d_id]))
  })), sample(alphabet, 1))
}

#' generate multiple motifs from alphabet
#' @param alphabet elements used to build motif
#' @param n_motif number of motifs to be generated
#' @return list of generated motifs
#' @export
#' @examples
#' generate_motifs(1:4, 5)
generate_motifs <- function(alphabet, n_motif) {
  lapply(1L:n_motif, function(dummy) generate_single_motif(alphabet))
}
