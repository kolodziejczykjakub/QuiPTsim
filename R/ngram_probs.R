
#' function computes transition matrix
#' @param sequences matrix of sequences (each row represents single sequence)
#' @param alphabet elements used to build sequences
#' @return combined ngram matrix
#' @export

sequenceTransitionMatrix <- function(sequences, alphabet) {
  seqList <- split(sequences, rep(1:nrow(sequences), each = ncol(sequences)))
  transitionMatrix <- markovchainFit(seqList, possibleStates = alphabet)

  new("markovchain", states = alphabet,
      transitionMatrix = attr(transitionMatrix$estimate, "transitionMatrix"),
      name = "Markov Chain object")
}

calculate_ngram_prob <- function(mc, seq) {

  alphabet <- attr(mc, "states")
  #browser()

  #asserts start and end as letter
  # assert len > 1

  seqlist <- as.list(seq)
  seqlist_gaps <- lapply(seq, function(x) {
    if (x == "_") {
      alphabet
    } else {
      x
    }
  })

  seqGrid <- expand.grid(seqlist_gaps)
  sum(apply(seqGrid, 1, function(seq) calculate_seq_prob(mc, seq)))

}

calculate_seq_prob <- function(mc, seq) {

  prob <- 1
  for (i in 2:length(seq)) {
    prob <- prob * transitionProbability(mc, seq[i-1], seq[i])
  }

  prob
}
