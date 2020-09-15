#' function generates motif from an alphabet
#' TODO: motif length control
#' @param alphabet elements used to build motif
#' @param n number of elements from alphabet
#' @param d total sum of gaps
#' @param motifProbs alphabet probabilites for motifs
#' @return motif based on given alphabet
#' @export
#' @examples
#' generate_single_motif(1:4)
#' generate_single_motif(c("a", "b", "c"))
#' generate_single_motif(1:4, n = 6, d = 2, motifProbs = c(0.7, 0.1, 0.1, 0.1))

generate_single_motif <- function(alphabet, n = 4, d = 6, motifProbs = NULL) {

  motif <- sample(alphabet, sample(2:n, 1), replace = TRUE, prob = motifProbs)
  gaps <- expand.grid(list(0:d)[rep(1, length(motif) - 1)])
  possibleGaps <- gaps[apply(gaps, 1, sum) <= d, , drop = FALSE]
  gap <- possibleGaps[sample(1:nrow(possibleGaps), 1), , drop = FALSE]

  newSeq <- vector()
  for (i in 1:(length(motif) -1)) {
    newSeq <- c(newSeq, motif[i], rep("_", gap[1, i]))
  }
  newSeq <- c(newSeq, motif[length(motif)])

  newSeq
}

#' generate multiple motifs from alphabet
#' @param alphabet elements used to build motif
#' @param n_motif number of motifs to be generated
#' @param motifProbs alphabet probabilites for motifs
#' @param n maximum number of alphabet elements in motif
#' @param d maximum number of gaps in motif
#' @return list of generated motifs
#' @export
#' @examples
#' generate_motifs(1:4, 5)
#' generate_motifs(1:4, 5, n = 6, d = 2, motifProbs = c(0.7, 0.1, 0.1, 0.1))
generate_motifs <- function(alphabet, n_motif, n = 4, d = 6, motifProbs = NULL) {
  # check if generated motifs can be injected to sequence of length 10
  validated <- FALSE
  while (!validated) {
    motifs <- lapply(1L:n_motif, function(dummy) generate_single_motif(alphabet, n, d, motifProbs))
    validated <- validate_motifs(motifs, 10)
  }
  motifs
}

#' generates sequence of elements from alphabet with replacement
#' @param alphabet elements used to build sequence
#' @param len length of generated sample sequence
#' @param seqProbs alphabet probabilites for sequences
#' @return randomly generated sequence
#' @export
#' @examples
#' simulate_single_sequence(5, 1L:4)
#' simulate_single_sequence(10, c("a", "b", "c"))
#' simulate_single_sequence(10, c("a", "b", "c"), c(0.6, 0.2, 0.2))

simulate_single_sequence <- function(len, alphabet, seqProbs = NULL){
  sample(alphabet, size = len, replace = TRUE, prob = seqProbs)
}

#' injects motifs to a sequence
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
#' example_sequence <- simulate_single_sequence(12, alph)
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

#' function generates sequences (both positive & negative)
#' @param n_seq number of sequences to be generated
#' @param len sequence length
#' @param alphabet elements used to build sequence
#' @param motifs_list list of injected motifs
#' @param n_motifs number of motifs injected to each positive sequence
#' @param fraction of positive sequences
#' @param seqProbs alphabet probabilites for sequences
#' @return generated sequences
#' @export
#' @examples
#' n_seq <- 20
#' len <- 10
#' alph <- 1L:4
#' motifs <- generate_motifs(alph, 3)
#' simulate_sequences(n_seq, len, alph, motifs, 1)
#' simulate_sequences(n_seq, len, alph, motifs, 1, fraction = 0.8)
#' simulate_sequences(n_seq, len, alph, motifs, 1, seqProbs = c(0.7, 0.1, 0.1, 0.1))

simulate_sequences <- function(n_seq,
                               len,
                               alphabet,
                               motifs_list,
                               n_motifs,
                               fraction = 0.5,
                               seqProbs = NULL) {
  n_pos <- round(fraction*n_seq, 0)

  list_of_motifs <- list()
  list_of_masks <- list()

  target <- logical(n_seq)
  target[1:n_pos] <- TRUE
  sequences <- matrix(nrow = n_seq, ncol = len)

  for (i in 1:n_pos) {
    motifs <- motifs_list[sample(1:length(motifs_list), n_motifs)]
    new_seq <- add_motifs(motifs, simulate_single_sequence(len, alphabet))
    list_of_motifs[[i]] <- attr(new_seq, "motifs")
    list_of_masks[[i]] <- attr(new_seq, "masks")
    sequences[i, ] <- new_seq
  }
  for (i in 1:(n_seq - n_pos)) {
    sequences[n_pos + i, ] <- simulate_single_sequence(len, alphabet, seqProbs)
  }
  attr(sequences, "motifs") <- list_of_motifs
  attr(sequences, "masks") <- list_of_masks
  attr(sequences, "target") <- target
  sequences
}

#' function counts n-grams in given sequences
#' @param n_seq number of sequences to be generated
#' @param l_seq sequence length
#' @param alphabet elements used to build sequence
#' @param motifs_list list of injected motifs
#' @param n_motifs number of motifs injected to each positive sequence
#' @param fraction TODO: add fraction: of positive sequences / change approach
#' @param seqProbs alphabet probabilites for sequences
#' @param n maximum number of alphabet elements in n-gram
#' @param d maximum number of gaps in n-gram
#' @return generated sequences
#' @export
#' @importFrom seqR count_multimers
#' @examples
#' n_seq <- 20
#' len <- 1200
#' alph <- letters[1:4]
#' motifs <- generate_motifs(alph, 2)
#' results <- generate_sequences(n_seq, len, alph, motifs, 1)
#' results <- generate_sequences(n_seq, len, alph, motifs, 1, seqProbs = c(0.7, 0.1, 0.1, 0.1))
#' results
#' attributes(results)
#'
generate_sequences <- function(n_seq,
                               l_seq,
                               alphabet,
                               motifs_list,
                               n_motifs,
                               fraction = 0.5,
                               seqProbs = NULL,
                               n = 4,
                               d = 6) {
  # generate sequence data
  test_dat <- simulate_sequences(n_seq, l_seq, alphabet, motifs_list, n_motifs, fraction, seqProbs)


  # element & gaps' positions for count_multigrams
  ns <- c()
  ds <- c()
  for (i in 1:(n-1)) {
    ds_ <- expand.grid(list(0:d)[rep(1, i)])
    ds_ <- ds_[apply(ds_, 1, sum) <= d, , drop = FALSE]
    ns <- c(ns, rep(i+1, nrow(ds_)))
    ds <- c(ds, split(ds_, 1:nrow(ds_)))
  }
  ds <- lapply(ds, unlist)

  test_res <- count_multimers(test_dat,
                              c(1, ns),
                              alphabet,
                              kmer_gaps_list = c(list(c()), ds),
                              with_kmer_counts = FALSE)



  attr(test_res, "sequences") <- matrix(test_dat, nrow = nrow(test_dat), ncol = ncol(test_dat))
  attr(test_res, "motifs") <- attr(test_dat, "motifs")
  attr(test_res, "masks") <- attr(test_dat, "masks")
  attr(test_res, "target") <- attr(test_dat, "target")
  test_res
}

#' function validates if given motifs can be injected to a sequence of given length
#' @param motifs list of motifs we are checking
#' @param sequence_length length of sequence we want to inject
#' @return logical value if such injection is possible
#' @export
#' @examples
#' set.seed(42)
#' motifs <- generate_motifs(1:4, 3, d = 3)
#' validate_motifs(motifs, 7)
#' validate_motifs(motifs, 9)
validate_motifs <- function(motifs, sequence_length) {
  result <- tryCatch(add_motifs(motifs, rep("*", sequence_length)),
                     error = function(dummy) FALSE)
  ifelse(class(result) == "character", TRUE, FALSE)
}
