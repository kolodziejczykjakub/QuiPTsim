library(QuiPTsim)
context("Ranking functions")

test_that("kmers_for_nonranking_methods", {

  set.seed(42)
  example_filtering_results <- data.frame(score = runif(100))
  res_QuiPT <- kmers_for_nonranking_methods(example_filtering_results, "QuiPT", c(0.01, 0.05))
  res_Chi <- kmers_for_nonranking_methods(example_filtering_results,
                                          "Chi-squared", c(0.01, 0.05))

  res_gainratio <- kmers_for_nonranking_methods(example_filtering_results, "gainratio")
  res_infogain <- kmers_for_nonranking_methods(example_filtering_results, "infogain")
  res_symuncert <- kmers_for_nonranking_methods(example_filtering_results, "symuncert")

  FCBF_filtering_results = data.frame(score = c(rep(0, 50), rep(1, 50)))
  res_FCBF <- kmers_for_nonranking_methods(FCBF_filtering_results, "FCBF")

  expect_equal(c(5L, 8L), res_QuiPT)
  expect_equal(c(5L, 8L), res_Chi)
  expect_equal(c(cpt = 82), res_gainratio)
  expect_equal(c(cpt = 82), res_infogain)
  expect_equal(c(cpt = 82), res_symuncert)
  expect_equal(50, res_FCBF)
})

test_that("kmers_for_nonranking_methods", {

  set.seed(1)
  X_train <- data.frame(matrix(sample(0:1, 100, replace = T), ncol = 10))
  y_train <- c(rep(1, 5), rep(0, 5))
  X_test <- data.frame(matrix(sample(0:1, 100, replace = T), ncol = 10))
  y_test <- c(rep(1, 5), rep(0, 5))


  model_knn <- build_model(X_train, y_train, X_test, y_test, c(1, 2), "knn")
  model_lm <- build_model(X_train, y_train, X_test, y_test, 0, "lm")
  model_rf <- build_model(X_train, y_train, X_test, y_test, 100:101, "rf")
  model_nb <- build_model(X_train, y_train, X_test, y_test, 0, "naive bayes")


  expect_equal(model_knn,
               structure(c(0.666666666666667, 1, 0.5, 0.6, 1, 0, 0.333333333333333,
                           0.333333333333333, 0.5, 0, 0.666666666666667, 1, 0.5, 0.6, 0.6,
                           0.5, 0.333333333333333, 0.333333333333333, 0.5, 0.5), .Dim = c(10L,
                                                                                          2L)))
  expect_equal(model_lm,
               structure(c(1.24119544858363e-06, 5.301381376178e-06, 0.999999988948047,
                           1.55422722196123e-05, 2.35114517911265e-20, 8.65369394152444e-05,
                           7.85464234716512e-08, 0.0227888622381538, 8.267762676018e-23,
                           4.29606759781234e-08), .Dim = c(10L, 1L), .Dimnames = list(NULL,
                                                                                      "s0"), lambda = 0))
  expect_equal(model_rf,
               structure(c(0.522, 0.522, 0.522, 0.522, 0.522, 0.522, 0.522,
                           0.522, 0.522, 0.522, 0.51980198019802, 0.51980198019802, 0.51980198019802,
                           0.51980198019802, 0.51980198019802, 0.51980198019802, 0.51980198019802,
                           0.51980198019802, 0.51980198019802, 0.51980198019802), .Dim = c(10L,
                                                                                           2L)))
  expect_equal(model_nb,
               structure(c(0.889105085232103, 0.829435580770777, 0.999998598834015,
                           0.714015228751967, 2.71982681933708e-06, 0.561757797943185, 0.396907433961698,
                           0.889105085232103, 7.16940001698361e-07, 0.437404884724875), .Dim = c(10L,
                                                                                                 1L), .Dimnames = list(NULL, "1")))
})
# test_that("kmers_for_nonranking_methods", {
#
#
# })
