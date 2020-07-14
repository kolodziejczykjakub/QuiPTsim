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

#' generates sequence of elements from alphabet with replacement
#' @param alphabet elements used to build sequence
#' @param len length of generated sample sequence
#' @return randomly generated sequence
#' @export
#' @examples
#' simulate_single_sequence(5, 1L:4)
#' simulate_single_sequence(10, c("a", "b", "c"))

simulate_single_sequence <- function(len, alphabet){
  sample(alphabet, size = len, replace = TRUE)
}

#' injects motifs to a sequence
#' TODO: examples
#' @param motifs list of motifs to be injected
#' @param sequence vector of alphabet elements
#' @return list(sequence, motifs, masks)
#' @export
#' @examples
#' # simple injection
#' add_motifs(list(c(1, "_", 1), c(1, 1)), c(2, 2, 3, 4))
#' # little bit more interesting
#' alph <- as.character(1L:4)
#' motifs <- generate_motifs(alph, 2)
#' example_sequence <- simulate_single_sequence(10, alph)
#' add_motifs(motifs, example_sequence)

add_motifs <- function(motifs, sequence) {
  sequence_len <- length(sequence)

  #create grid of possible motifs' positions
  maximum_motifs_positions <- lapply(motifs, function(x)
    seq(sequence_len - length(x) + 1))
  motifs_grid <- expand.grid(maximum_motifs_positions)
  motifs_grid <- motifs_grid[sample(1:nrow(motifs_grid)), , drop = FALSE]

  for (i in 1:nrow(motifs_grid)) {

    list_of_masks <- list()
    injected_sequence <- sequence
    injected_positions <- logical(length(sequence))

    for (j in 1:ncol(motifs_grid)) {
      mask <- rep(FALSE, sequence_len)
      new_injected_sequence <- injected_sequence
      motif <- motifs[[j]]
      ids <- 0:(length(motif) - 1)
      ids <- ids[motif != "_"] + motifs_grid[i, j]
      mask[ids] <- TRUE
      new_injected_sequence[ids] <- motif[motif != "_"]

      if (j == 1) {
        injected_sequence <- new_injected_sequence
        injected_positions <- mask
      } else {
        if (all(injected_sequence[injected_positions] == new_injected_sequence[injected_positions])){
          injected_sequence <- new_injected_sequence
          injected_positions <- (injected_positions | mask)
        } else {
          break
        }
      }

      list_of_masks[[j]] <- mask

      if (j == ncol(motifs_grid)){
        attr(injected_sequence, "motifs") <- motifs
        attr(injected_sequence, "masks") <- list_of_masks
        return(injected_sequence)
      }
    }
  }
  stop("Given motifs cannot be injected to a sequence!")
}

