library(QuiPTsim)
library(drake)

plan <- drake_plan(

  ###### data frame containing datasets' details and their paths to RDS files
  df = read.csv("~/projects/QuiPTsim-data/reduced_alph_enc_amylogram_encoding_unigram/amylogram_encoding.csv"),

  ###### selected paths
  # workaround after adding unigrams
  paths1 = paste0("~/projects/QuiPTsim-data/reduced_alph_enc_amylogram_encoding_unigram/",
		sapply(strsplit(x = df[["path"]], split = "/"), function(x) x[[3]])),
  paths = paths1[1:2],
  ###### details of models used in ranking comparison
  models_details = list(
    list(model = "lm",
         param_name = "lambda",
         param_value = NULL),
    list(model = "knn",
         param_name = "neighbors",
         param_value = 2^(0:4)),
    list(model = "rf",
         param_name = "num.trees",
         param_value = c(500, 1000)),
    list(model = "naive bayes",
         param_name = "laplace",
         param_value = 0)
  ),

  ###### Validation scheme
  validation_scheme = list(type = "cv",
                            folds = 5,
                            n_kmers = c(2:128, 2^(8:13)),
                            cv_reps = 3,
                            models_details = models_details),

  validation_scheme_nonranking = list(type = "cv",
                           folds = 5,
                           cv_reps = 3,
                           models_details = models_details),


  ###### Validation scheme
  filter_names = c("QuiPT",
                    "Chi-squared",
                    "FCBF",
                    "infogain", "gainratio", "symuncert",
                    "MRMR", "JMI", "JMIM", "DISR", "CMIM", "NJMIM"),
  ###### Filtering results
  # QuiPT
  filtering_results_QuiPT = lapply(paths, function(x)
    filter_ngrams(read_ngram_matrix(x, n = 100, fraction = 0.5), feature_selection_method = "QuiPT")),

  # Chi-squared test
  filtering_results_Chi = lapply(paths, function(x)
    filter_ngrams(read_ngram_matrix(x, n = 100, fraction = 0.5), feature_selection_method = "Chi-squared")),

  #  FCBF
  filtering_results_FCBF = lapply(paths, function(x)
    filter_ngrams(read_ngram_matrix(x, n = 100, fraction = 0.5), feature_selection_method = "FCBF")),

  # "infogain", "gainratio", "symuncert"
  filtering_results_infogain = lapply(paths, function(x)
    filter_ngrams(read_ngram_matrix(x, n = 100, fraction = 0.5), feature_selection_method = "infogain")),

  filtering_results_gainratio = lapply(paths, function(x)
    filter_ngrams(read_ngram_matrix(x, n = 100, fraction = 0.5), feature_selection_method = "gainratio")),

  filtering_results_symuncert = lapply(paths, function(x)
    filter_ngrams(read_ngram_matrix(x, n = 100, fraction = 0.5), feature_selection_method = "symuncert")),

  # MRMR
  filtering_results_MRMR = lapply(paths, function(x)
    filter_ngrams(read_ngram_matrix(x, n = 100, fraction = 0.5), feature_selection_method = "MRMR")),

  filtering_results_JMI = lapply(paths, function(x)
    filter_ngrams(read_ngram_matrix(x, n = 100, fraction = 0.5), feature_selection_method = "JMI")),

  filtering_results_DISR = lapply(paths, function(x)
    filter_ngrams(read_ngram_matrix(x, n = 100, fraction = 0.5), feature_selection_method = "DISR")),

  filtering_results_JMIM = lapply(paths, function(x)
    filter_ngrams(read_ngram_matrix(x, n = 100, fraction = 0.5), feature_selection_method = "JMIM")),

  filtering_results_CMIM = lapply(paths, function(x)
    filter_ngrams(read_ngram_matrix(x, n = 100, fraction = 0.5), feature_selection_method = "CMIM")),


  filtering_results_NJMIM = lapply(paths, function(x)
    filter_ngrams(read_ngram_matrix(x, n = 100, fraction = 0.5), feature_selection_method = "NJMIM")),

  ###### Classifiers on selected k-mers

  results_QuiPT = lapply(1:length(paths), function(i) {
    evaluate_filtering_results(read_ngram_matrix(paths[[i]]),
                               filtering_results_QuiPT[[i]],
                               list(method="QuiPT"),
                               validation_scheme)
  }),

  results_Chi = lapply(1:length(paths), function(i) {
    evaluate_filtering_results(read_ngram_matrix(paths[[i]]),
                               filtering_results_Chi[[i]],
                               list(method="Chi-squared"),
                               validation_scheme)
  }),

  results_infogain = lapply(1:length(paths), function(i) {
    evaluate_filtering_results(read_ngram_matrix(paths[[i]]),
                               filtering_results_infogain[[i]],
                               list(method="infogain"),
                               validation_scheme)
  }),

  results_gainratio = lapply(1:length(paths), function(i) {
    evaluate_filtering_results(read_ngram_matrix(paths[[i]]),
                               filtering_results_gainratio[[i]],
                               list(method="gainratio"),
                               validation_scheme)
  }),

  results_symuncert = lapply(1:length(paths), function(i) {
    evaluate_filtering_results(read_ngram_matrix(paths[[i]]),
                               filtering_results_symuncert[[i]],
                               list(method="symuncert"),
                               validation_scheme)
  }),

  results_MRMR = lapply(1:length(paths), function(i) {
    evaluate_filtering_results(read_ngram_matrix(paths[[i]]),
                               filtering_results_MRMR[[i]],
                               list(method="MRMR"),
                               validation_scheme)
  }),

  results_JMI = lapply(1:length(paths), function(i) {
    evaluate_filtering_results(read_ngram_matrix(paths[[i]]),
                               filtering_results_JMI[[i]],
                               list(method="JMI"),
                               validation_scheme)
  }),

  results_JMIM = lapply(1:length(paths), function(i) {
    evaluate_filtering_results(read_ngram_matrix(paths[[i]]),
                               filtering_results_JMIM[[i]],
                               list(method="JMIM"),
                               validation_scheme)
  }),

  results_DISR = lapply(1:length(paths), function(i) {
    evaluate_filtering_results(read_ngram_matrix(paths[[i]]),
                               filtering_results_DISR[[i]],
                               list(method="DISR"),
                               validation_scheme)
  }),

  results_CMIM = lapply(1:length(paths), function(i) {
    evaluate_filtering_results(read_ngram_matrix(paths[[i]]),
                               filtering_results_CMIM[[i]],
                               list(method="CMIM"),
                               validation_scheme)
  }),

  results_NJMIM = lapply(1:length(paths), function(i) {
    evaluate_filtering_results(read_ngram_matrix(paths[[i]]),
                               filtering_results_NJMIM[[i]],
                               list(method="NJMIM"),
                               validation_scheme)
  }),

    # results for non-ranking methods
  thresholds = c(0.01, 0.05),

  results_QuiPT_nonranking = lapply(1:length(filtering_results_QuiPT), function(i) {

    fr = filtering_results_QuiPT[[i]]
    validation_scheme_nonranking_tmp = validation_scheme_nonranking
    validation_scheme_nonranking_tmp[["n_kmers"]] =
      kmers_for_nonranking_methods(fr, "QuiPT", thresholds)

    m <- read_ngram_matrix(paths[[i]])
    evaluate_filtering_results(m,
                               fr,
                               list(method="QuiPT"),
                               validation_scheme_nonranking_tmp)

  }),

  results_Chi_nonranking = lapply(1:length(filtering_results_Chi), function(i) {

    fr = filtering_results_Chi[[i]]
    validation_scheme_nonranking_tmp = validation_scheme_nonranking
    validation_scheme_nonranking_tmp[["n_kmers"]] =
      kmers_for_nonranking_methods(fr, "Chi-squared", thresholds)

    m <- read_ngram_matrix(paths[[i]])
    evaluate_filtering_results(m,
                               fr,
                               list(method="Chi-squared"),
                               validation_scheme_nonranking_tmp)

  }),

  results_gainratio_nonranking = lapply(1:length(filtering_results_gainratio), function(i) {

    fr = filtering_results_gainratio[[i]]
    validation_scheme_nonranking_tmp = validation_scheme_nonranking
    validation_scheme_nonranking_tmp[["n_kmers"]] =
      kmers_for_nonranking_methods(fr, "gainratio")

    m <- read_ngram_matrix(paths[[i]])
    evaluate_filtering_results(m,
                               fr,
                               list(method="gainratio"),
                               validation_scheme_nonranking_tmp)

  }),

  results_infogain_nonranking = lapply(1:length(filtering_results_infogain), function(i) {

    fr = filtering_results_infogain[[i]]
    validation_scheme_nonranking_tmp = validation_scheme_nonranking
    validation_scheme_nonranking_tmp[["n_kmers"]] =
      kmers_for_nonranking_methods(fr, "infogain")

    m <- read_ngram_matrix(paths[[i]])
    evaluate_filtering_results(m,
                               fr,
                               list(method="infogain"),
                               validation_scheme_nonranking_tmp)

  }),

  results_symuncert_nonranking = lapply(1:length(filtering_results_symuncert), function(i) {

    fr = filtering_results_symuncert[[i]]
    validation_scheme_nonranking_tmp = validation_scheme_nonranking
    validation_scheme_nonranking_tmp[["n_kmers"]] =
      kmers_for_nonranking_methods(fr, "symuncert")

    m <- read_ngram_matrix(paths[[i]])
    evaluate_filtering_results(m,
                               fr,
                               list(method="symuncert"),
                               validation_scheme_nonranking_tmp)

  }),

  results_FCBF_nonranking = lapply(1:length(filtering_results_FCBF), function(i) {

    fr = filtering_results_FCBF[[i]]
    validation_scheme_nonranking_tmp = validation_scheme_nonranking
    validation_scheme_nonranking_tmp[["n_kmers"]] =
      kmers_for_nonranking_methods(fr, "FCBF")

    m <- read_ngram_matrix(paths[[i]])
    evaluate_filtering_results(m,
                               fr,
                               list(method="FCBF"),
                               validation_scheme_nonranking_tmp)

  })
)


# supressing warnings from Chi2 and AUC
oldw <- getOption("warn")
options(warn = -1)

make(
  plan,
  parallelism = "future",
  jobs = 8,
  log_make = "drake.log"
)
options(warn = oldw)
