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
test_that("kmers_for_nonranking_methods", {

  set.seed(1)
  models_details = list(
    list(model = "lm",
         param_name = "lambda",
         param_value = 0),
    list(model = "knn",
         param_name = "neighbors",
         param_value = 2),
    list(model = "rf",
         param_name = "num.trees",
         param_value = 50),
    list(model = "naive bayes",
         param_name = "laplace",
         param_value = 0))

  validation_scheme = list(type = "cv",
                           folds = 2,
                           n_kmers = 2,
                           cv_reps = 1,
                           models_details = models_details)

  df <- data.frame(matrix(sample(0:1, 1000, replace = T), ncol = 10), y = c(rep(1, 50), rep(0, 50)))

  evaluated_kmers <- evaluate_selected_kmers(df, validation_scheme)

  expect_equal(evaluated_kmers,
               structure(list(
                 model = c("lm", "knn", "rf", "naive bayes", "lm",
                           "knn", "rf", "naive bayes"),
                 param = c("lambda", "neighbors", "num.trees", "laplace",
                           "lambda", "neighbors", "num.trees", "laplace")
                 , value = c(0, 2, 50, 0, 0, 2, 50, 0),
                 accuracy = structure(list(0.38, 0.84, 0.44, 0.46, 0.54, 0.92, 0.44, 0.44),
                                      .Names = c(NA, "accuracy", "accuracy",
                                                 "", NA, "accuracy", "accuracy", "")),
                 sensitivity = structure(list(0.36, 0.68, 0.44, 0.4, 0.56,
                                              0.84, 0.6, 0.52),
                                         .Names = c(NA, "sensitivity", "sensitivity",
                                                    "", NA, "sensitivity", "sensitivity", "")),
                 specificity = structure(list(0.4, 1, 0.44, 0.52, 0.52, 1, 0.28, 0.36),
                                         .Names = c(NA, "specificity", "specificity",
                                                    "", NA, "specificity", "specificity", "")),
                 F1score = structure(list(0.36734693877551, 0.80952380952381, 0.44,
                                          0.425531914893617, 0.549019607843137,
                                          0.91304347826087, 0.517241379310345,
                                          0.481481481481481),
                                     .Names = c(NA, "F1score", "F1score", "",
                                                NA, "F1score", "F1score", "")),
                 precision = structure(list(0.375, 1, 0.44, 0.454545454545455,
                                            0.538461538461538, 1, 0.454545454545455,
                                            0.448275862068966),
                                       .Names = c(NA, "precision", "precision", "",
                                                  NA, "precision", "precision", "")),
                 recall = structure(list(0.36, 0.68, 0.44, 0.4, 0.56, 0.84, 0.6, 0.52),
                                    .Names = c(NA, "recall", "recall", "", NA, "recall",
                                               "recall", "")),
                 auc = structure(list(0.5792, 0.9424, 0.5584, 0.5224, 0.5152, 0.984,
                                      0.5728, 0.528), .Names = c(NA, "auc", "auc", "",
                                                                 NA, "auc", "auc", "")),
                 fold = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L)),
                 class = "data.frame", row.names = c(NA, -8L)))
})
