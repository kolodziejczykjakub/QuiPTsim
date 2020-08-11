context("Probabilites of n-grams")

test_that("Transition matrix generation", {

  set.seed(1)
  alphabet <- letters[1:4]
  sequences <- matrix(sample(alphabet, size = 50, replace=TRUE), nrow = 5, ncol = 10)
  mc <- sequenceTransitionMatrix(sequences, alphabet)
  expect_equal(mc, new("markovchain", states = c("a", "b", "c", "d"), byrow = TRUE,
                       transitionMatrix = structure(c(0.375, 0.230769230769231,
                                                      0.444444444444444, 0.142857142857143, 0.25, 0.538461538461538,
                                                      0.222222222222222, 0.142857142857143, 0.1875, 0.153846153846154,
                                                      0.222222222222222, 0.428571428571429, 0.1875, 0.0769230769230769,
                                                      0.111111111111111, 0.285714285714286), .Dim = c(4L, 4L), .Dimnames = list(
                                                        c("a", "b", "c", "d"), c("a", "b", "c", "d"))), name = "Markov Chain object"))
})

test_that("Calculating ngram probabilites", {

  set.seed(1)
  alphabet <- letters[1:4]
  sequences <- matrix(sample(alphabet, size = 50, replace=TRUE), nrow = 5, ncol = 10)
  mc <- sequenceTransitionMatrix(sequences, alphabet)
  example_ngram <- sample(c(alphabet, "_"), size = 5, replace = TRUE)
  prob <- calculate_ngram_prob(mc, example_ngram)
  expect_equal(prob, 0.00103021978021978)

})

test_that("Calculating ngram probabilites", {

  set.seed(1)
  alphabet <- letters[1:4]
  sequences <- matrix(sample(alphabet, size = 50, replace=TRUE), nrow = 5, ncol = 10)
  mc <- sequenceTransitionMatrix(sequences, alphabet)
  example_seq <- sample(alphabet, size = 5, replace = TRUE)
  prob <- calculate_seq_prob(mc, example_seq)
  expect_equal(prob, 0.00176609105180534)
})

