#' function computes transition matrix
#' @param sequences matrix of sequences (each row represents single sequence)
#' @param alphabet elements used to build sequences
#' @return combined ngram matrix
#' @importFrom markovchain markovchainFit
#' @importFrom methods new
#' @export
#' @examples
#' alphabet <- letters[1:4]
#' sequences <- matrix(sample(alphabet, size = 50, replace=TRUE), nrow = 5, ncol = 10)
#' sequenceTransitionMatrix(sequences, alphabet)
#'
sequenceTransitionMatrix <- function(sequences, alphabet) {
  seqList <- split(sequences, rep(1:nrow(sequences), each = ncol(sequences)))
  transitionMatrix <- markovchainFit(seqList, possibleStates = alphabet)

  new("markovchain", states = alphabet,
      transitionMatrix = attr(transitionMatrix$estimate, "transitionMatrix"),
      name = "Markov Chain object")
}

#' function computes probability of given n-gram
#' @param mc markovchain object containing transition matrix
#' @param seq sequence of alphabet elements or gaps("_")
#' @export
#' @examples
#' alphabet <- letters[1:4]
#' sequences <- matrix(sample(alphabet, size = 50, replace=TRUE), nrow = 5, ncol = 10)
#' mc <- sequenceTransitionMatrix(sequences, alphabet)
#' example_ngram <- sample(c(alphabet, "_"), size = 5, replace = TRUE)
#' print(example_ngram)
#' calculate_ngram_prob(mc, example_ngram)

calculate_ngram_prob <- function(mc, seq) {

  if (length(seq) < 2)
    stop("Error: Sequence is too short!")

  if (!(class(mc) == "markovchain"))
    stop("Error: mc parameter is not a markovchain object instance")

  alphabet <- attr(mc, "states")

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

#' function computes probability of given sequence
#' @param mc markovchain object containing transition matrix
#' @param seq sequence of alphabet elements
#' @importFrom markovchain transitionProbability
#' @export
#' @examples
#' alphabet <- letters[1:4]
#' sequences <- matrix(sample(alphabet, size = 50, replace=TRUE), nrow = 5, ncol = 10)
#' mc <- sequenceTransitionMatrix(sequences, alphabet)
#' example_seq <- sample(alphabet, size = 5, replace = TRUE)
#' print(example_seq)
#' calculate_seq_prob(mc, example_seq)

calculate_seq_prob <- function(mc, seq) {

  prob <- 1
  for (i in 2:length(seq)) {
    prob <- prob * transitionProbability(mc, seq[i-1], seq[i])
  }

  prob
}
