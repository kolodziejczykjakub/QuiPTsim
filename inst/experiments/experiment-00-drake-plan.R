library(QuiPTsim)
library(drake)

# TODO:
# - parallelization

plan <- drake_plan(

  ###### data frame containing datasets' details and their paths to RDS files
  df = read.csv("./reduced_alph_enc_alph4_const/alph4_const.csv"),

  ###### selected paths
  paths = df[df$l_seq == 10 & df$n_motifs == 2, "path"],

  ###### details of models used in ranking comparison
  models_details = list(
    list(model = "lm",
         param_name = "lambda",
         param_value = NULL),
    list(model = "knn",
         param_name = "neighbors",
         param_value = 2^(0:4)),
    list(model = "rf",
         param_name = "num.trees", # 500, 1000
         param_value = c(500, 100)),
    list(model = "naive bayes",
         param_name = "laplace",
         param_value = 0)
  ),

  ###### Validation scheme
  validation_scheme = list(type = "cv",
                            folds = 5,
                            n_kmers = c(1:128, 2^(8:13)),
                            cv_reps = 2,
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
    filter_ngrams(read_ngram_matrix(x), feature_selection_method = "QuiPT")),

  # Chi-squared test
  filtering_results_Chi = lapply(paths, function(x)
    filter_ngrams(read_ngram_matrix(x), feature_selection_method = "Chi-squared")),

  #  FCBF
  filtering_results_FCBF = lapply(paths, function(x)
    filter_ngrams(read_ngram_matrix(x), feature_selection_method = "FCBF")),

  # "infogain", "gainratio", "symuncert"
  filtering_results_infogain = lapply(paths, function(x)
    filter_ngrams(read_ngram_matrix(x), feature_selection_method = "infogain")),

  filtering_results_gainratio = lapply(paths, function(x)
    filter_ngrams(read_ngram_matrix(x), feature_selection_method = "gainratio")),

  filtering_results_symuncert = lapply(paths, function(x)
    filter_ngrams(read_ngram_matrix(x), feature_selection_method = "symuncert")),

  # MRMR
  filtering_results_MRMR = lapply(paths, function(x)
    filter_ngrams(read_ngram_matrix(x), feature_selection_method = "MRMR")),

  filtering_results_JMI = lapply(paths, function(x)
    filter_ngrams(read_ngram_matrix(x), feature_selection_method = "JMI")),

  filtering_results_DISR = lapply(paths, function(x)
    filter_ngrams(read_ngram_matrix(x), feature_selection_method = "DISR")),

  filtering_results_JMIM = lapply(paths, function(x)
    filter_ngrams(read_ngram_matrix(x), feature_selection_method = "JMIM")),

  filtering_results_CMIM = lapply(paths, function(x)
    filter_ngrams(read_ngram_matrix(x), feature_selection_method = "CMIM")),


  filtering_results_NJMIM = lapply(paths, function(x)
    filter_ngrams(read_ngram_matrix(x), feature_selection_method = "NJMIM")),

  ###### Classifiers on selected k-mers

  results_QuiPT = lapply(filtering_results_QuiPT, function(fr) {
    lapply(paths, function(path) {

      m <- read_ngram_matrix(path)
      evaluate_filtering_results(m, fr, list(method="QuiPT"), validation_scheme)

    })
  }),

  results_Chi = lapply(filtering_results_Chi, function(fr) {
    lapply(paths, function(path) {

      m <- read_ngram_matrix(path)
      evaluate_filtering_results(m, fr, list(method="Chi"), validation_scheme)

    })
  }),

  results_FCBF = lapply(filtering_results_FCBF, function(fr) {
    lapply(paths, function(path) {

      m <- read_ngram_matrix(path)
      evaluate_filtering_results(m, fr, list(method="FCBF"), validation_scheme)

    })
  }),

  results_infogain = lapply(filtering_results_infogain, function(fr) {
    lapply(paths, function(path) {

      m <- read_ngram_matrix(path)
      evaluate_filtering_results(m, fr, list(method="infogain"), validation_scheme)

    })
  }),

  results_gainratio = lapply(filtering_results_gainratio, function(fr) {
    lapply(paths, function(path) {

      m <- read_ngram_matrix(path)
      evaluate_filtering_results(m, fr, list(method="gainratio"), validation_scheme)

    })
  }),

  results_symuncert = lapply(filtering_results_symuncert, function(fr) {
    lapply(paths, function(path) {

      m <- read_ngram_matrix(path)
      evaluate_filtering_results(m, fr, list(method="symuncert"), validation_scheme)

    })
  }),

  results_MRMR = lapply(filtering_results_MRMR, function(fr) {
    lapply(paths, function(path) {

      m <- read_ngram_matrix(path)
      evaluate_filtering_results(m, fr, list(method="MRMR"), validation_scheme)

    })
  }),

  results_JMI = lapply(filtering_results_JMI, function(fr) {
    lapply(paths, function(path) {

      m <- read_ngram_matrix(path)
      evaluate_filtering_results(m, fr, list(method="JMI"), validation_scheme)

    })
  }),

  results_JMIM = lapply(filtering_results_JMIM, function(fr) {
    lapply(paths, function(path) {

      m <- read_ngram_matrix(path)
      evaluate_filtering_results(m, fr, list(method="JMIM"), validation_scheme)

    })
  }),

  results_DISR = lapply(filtering_results_DISR, function(fr) {
    lapply(paths, function(path) {

      m <- read_ngram_matrix(path)
      evaluate_filtering_results(m, fr, list(method="DISR"), validation_scheme)

    })
  }),

  results_CMIM = lapply(filtering_results_CMIM, function(fr) {
    lapply(paths, function(path) {

      m <- read_ngram_matrix(path)
      evaluate_filtering_results(m, fr, list(method="CMIM"), validation_scheme)

    })
  }),

  results_NJMIM = lapply(filtering_results_NJMIM, function(fr) {
    lapply(paths, function(path) {

      m <- read_ngram_matrix(path)
      evaluate_filtering_results(m, fr, list(method="NJMIM"), validation_scheme)

    })
  })
)

#make(plan)
vis_drake_graph(plan)
