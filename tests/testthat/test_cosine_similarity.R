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
  expect_equal(probs, structure(c(0.143778152017506, 0.0241113304194032, 0.147278182280115,
                                  0.101360968594866, 0.0449719520582651, 0.250034011037696, 0.130519716033949,
                                  -0.0271629944044605, 0.100861800627489, 0.0601686657281429, 0.081585224339416,
                                  0.222098633054888, 0.11576804441742, 0.20619663337668, 0.0211652673558396,
                                  0.0188998941585601, 0.103258103808209, 0.062125428420114, 0.110813557061791,
                                  0.0821674296141103), .Dim = c(2L, 10L), .Dimnames = list(c("u",
                                                                                             "w"), NULL)))
})
