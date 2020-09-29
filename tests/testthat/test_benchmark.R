context("Benchmark tests")

test_that("QuiPT benchmark", {

  path <- c("./example_QuiPT_simulation_data_result.Rds")
  bm_QuiPT <- create_benchmark_data(path, list(method = "QuiPT"))

  results_QuiPT <- benchmark_summary(bm_QuiPT, list(method = "QuiPT",
                                                    pval_thresholds = c(0.05, 0.01),
                                                    pval_adjustments = c("", "BH")))

  expect_equal(results_QuiPT,
               structure(list(pval_thresholds = c(0.05, 0.01, 0.05, 0.01),
                              pval_adjustments = structure(c(1L, 1L, 2L, 2L),
                                                           .Label = c("", "BH"),
                                                           class = "factor"),
                              sensitivity_mean = c(0.777777777777778,0.722222222222222, 0.722222222222222, 0.722222222222222),
                              specificity_mean = c(0.82821908569777, 0.890794354255985, 0.90158485476179, 0.923214027650658),
                              F1score_mean = c(0.00778210116731518, 0.0113141862489121, 0.0125361620057859, 0.016),
                              precision_mean = c(0.00391061452513967, 0.00570175438596491, 0.00632295719844358, 0.00808960796515246),
                              recall_mean = c(0.777777777777778, 0.722222222222222, 0.722222222222222, 0.722222222222222),
                              sensitivity_std = c(NA_real_, NA_real_, NA_real_, NA_real_),
                              specificity_std = c(NA_real_, NA_real_, NA_real_, NA_real_),
                              F1score_std = c(NA_real_, NA_real_, NA_real_, NA_real_),
                              precision_std = c(NA_real_, NA_real_, NA_real_, NA_real_),
                              recall_std = c(NA_real_, NA_real_, NA_real_, NA_real_)), out.attrs = list(
                                dim = c(pval_thresholds = 2L, pval_adjustments = 2L), dimnames = list(
                                  pval_thresholds = c("pval_thresholds=0.05", "pval_thresholds=0.01"
                                  ), pval_adjustments = c("pval_adjustments=", "pval_adjustments=BH"
                                  ))), row.names = c(NA, -4L), class = "data.frame", setup = list(
                                    method = "QuiPT", pval_thresholds = c(0.05, 0.01), pval_adjustments = c("", "BH"))))
})

test_that("FCBF benchmark", {

  path <- c("./example_QuiPT_simulation_data_result.Rds")
  bm_FCBF <- create_benchmark_data(path, list(method = "FCBF"))

  results_FCBF <- benchmark_summary(bm_FCBF, list(method = "FCBF"))

  expect_equal(results_FCBF,
               structure(list(sensitivity_mean = 0.111111111111111, specificity_mean = 0.999181078086613,
                              F1score_mean = 0.108108108108108, precision_mean = 0.105263157894737,
                              recall_mean = 0.111111111111111, sensitivity_std = NA_real_,
                              specificity_std = NA_real_, F1score_std = NA_real_, precision_std = NA_real_,
                              recall_std = NA_real_), class = "data.frame", row.names = c(NA,
                                                                                          -1L), setup = list(method = "FCBF")))
})

test_that("Chi-squared benchmark", {

  path <- c("./example_QuiPT_simulation_data_result.Rds")
  bm_chi <- create_benchmark_data(path, list(method = "Chi-squared"))

  results_chi <- benchmark_summary(bm_chi, list(method = "Chi-squared",
                                                pval_thresholds = 0.1,
                                                pval_adjustments = "BH"))

  expect_equal(results_chi,
               structure(list(pval_thresholds = 0.1, pval_adjustments = structure(1L, .Label = "BH", class = "factor"),
                              sensitivity_mean = 0.777777777777778, specificity_mean = 0.823883616744545,
                              F1score_mean = 0.00759219088937093, precision_mean = 0.00381471389645777,
                              recall_mean = 0.777777777777778, sensitivity_std = NA_real_,
                              specificity_std = NA_real_,
                              F1score_std = NA_real_,
                              precision_std = NA_real_,
                              recall_std = NA_real_),
                         out.attrs = list(dim = c(pval_thresholds = 1L,
                                                  pval_adjustments = 1L),
                                          dimnames = list(pval_thresholds = "pval_thresholds=0.1",
                                                          pval_adjustments = "pval_adjustments=BH")),
                         row.names = c(NA, -1L),
                         class = "data.frame",
                         setup = list(method = "Chi-squared",
                                      pval_thresholds = 0.1, pval_adjustments = "BH")))
})

test_that("FSelectorRcpp benchmark", {

  path <- c("./example_QuiPT_simulation_data_result.Rds")
  #path <- c("./tests/testthat/example_QuiPT_simulation_data_result.Rds")
  bm_fselector <- create_benchmark_data(path, list(method = "gainratio"))

  results_fselector <- benchmark_summary(bm_fselector,
                                         list(method = "gainratio",
                                              fraction = 0.025))

  expect_equal(results_fselector,
               structure(list(num_features = 519, top_fraction = 0.025, sensitivity_mean = 0.722222222222222,
                              specificity_mean = 0.976299436389036, F1score_mean = 0.0497131931166348,
                              precision_mean = 0.0257425742574257, recall_mean = 0.722222222222222,
                              sensitivity_std = NA_real_, specificity_std = NA_real_, F1score_std = NA_real_,
                              precision_std = NA_real_, recall_std = NA_real_), class = "data.frame", row.names = c(NA,
                                                                                                                    -1L), setup = list(method = "gainratio", fraction = 0.025)))
})

test_that("MRMR benchmark", {

  path <- c("./example_QuiPT_simulation_data_result.Rds")

  # QuiPT
  bm_MRMR <- create_benchmark_data(path, list(method = "MRMR"))

  results_MRMR <- benchmark_summary(bm_MRMR, list(method = "MRMR"))

  expect_equal(results_MRMR,
               structure(list(sensitivity_mean = 1, specificity_mean = 0, F1score_mean = 0.00173118538110123,
                              precision_mean = 0.000866342590364345, recall_mean = 1, sensitivity_std = NA_real_,
                              specificity_std = NA_real_, F1score_std = NA_real_, precision_std = NA_real_,
                              recall_std = NA_real_), class = "data.frame", row.names = c(NA,
                                                                                          -1L), setup = list(method = "MRMR")))
})
