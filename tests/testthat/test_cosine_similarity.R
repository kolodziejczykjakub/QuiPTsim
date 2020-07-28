context("Cosine similarity")

test_that("Euclidean norm", {
  expect_true( euclidean_norm(1:10) - 19.62142 < 10e-5 )
})

test_that("Cosine similarity calculation", {
  expect_error(cosine_similarity(1, 1:4))
  expect_true( cosine_similarity(1:4, 5:8) - 0.9688639 < 10e-5 )
})

test_that("Probabilities generation", {
  set.seed(42)
  probs <- generate_probs(10, 0.625)
  expect_true(  cosine_similarity(probs[1,], probs[2,]) - 0.625 < 10e-5)
  expect_equal(probs, structure(c(0.096587795119055, 0.0124073614328669, 0.106567152951375,
                                  0.0739828281759341, 0.066253659102148, 0.197612868204078, 0.141957980992371,
                                  0.0103716555545052, 0.105632541137834, 0.179581976267067, 0.0990077169001684,
                                  0.250027774185551, 0.102606560647005, 0.0417350832596108, 0.136037776641259,
                                  0.0122294763370074, 0.0781929560445212, 0.00513954989484627,
                                  0.0671558604642629, 0.216911426688534), .Dim = c(2L, 10L), .Dimnames = list(
                                    c("u", "w"), NULL)))
})
