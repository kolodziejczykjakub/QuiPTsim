library(QuiPTsim)
library(drake)

plan <- drake_plan(
  
  ###### data frame containing datasets' details and their paths to RDS files
  #df = read.csv("~/projects/QuiPTsim-data/reduced_alph_enc_amylogram_encoding_unigram/amylogram_encoding.csv"),
  
  ###### selected paths
  #paths = df[df$seqLen==10 & df$nMotif=1 & nSeq = 300, "path"],
  paths = "./tests/testthat/example_QuiPT_simulation_data_result.Rds",
  output_prefix = "./do-usuniecia/",
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
                           n_kmers = 2:4,#c(2:128, 2^(8:12)),
                           cv_reps = 1,
                           models_details = models_details),
  
  validation_scheme_nonranking = list(type = "cv",
                                      folds = 5,
                                      cv_reps = 1,
                                      models_details = models_details),
  
  
  ###### Validation scheme
  filter_names = c("QuiPT",
                   "Chi-squared",
                   "FCBF",
                   "infogain", "gainratio", "symuncert",
                   "MRMR", "JMI", "JMIM", "DISR", "CMIM", "NJMIM"),
  
  ranking_QuiPT = filter_rankings(paths, output_prefix, "QuiPT", 300, 0.5, validation_scheme),
  ranking_Chi = filter_rankings(paths, output_prefix, "Chi-squared", 300, 0.5, validation_scheme),
  ranking_FCBF = filter_rankings(paths, output_prefix, "FCBF", 300, 0.5, validation_scheme),
  ranking_infogain = filter_rankings(paths, output_prefix, "infogain", 300, 0.5, validation_scheme),
  ranking_gainratio = filter_rankings(paths, output_prefix, "gainratio", 300, 0.5, validation_scheme),
  ranking_symuncert = filter_rankings(paths, output_prefix, "symuncert", 300, 0.5, validation_scheme),
  ranking_MRMR = filter_rankings(paths, output_prefix, "MRMR", 300, 0.5, validation_scheme),
  ranking_JMI = filter_rankings(paths, output_prefix, "JMI", 300, 0.5, validation_scheme),
  ranking_JMIM = filter_rankings(paths, output_prefix, "JMIM", 300, 0.5, validation_scheme),
  ranking_DISR = filter_rankings(paths, output_prefix, "DISR", 300, 0.5, validation_scheme),
  ranking_NJMIM = filter_rankings(paths, output_prefix, "NJMIM", 300, 0.5, validation_scheme),
  
  
  thresholds = c(0.01, 0.05),
  
  nonranking_QuiPT = filter_nonrankings(paths, output_prefix, "QuiPT", 300, 0.5, validation_scheme_nonranking,
                                         thresholds),
  nonranking_Chi = filter_nonrankings(paths, output_prefix, "Chi-squared", 300, 0.5, validation_scheme_nonranking,
                                       thresholds),
  
  nonranking_gainratio = filter_nonrankings(paths, output_prefix, "gainratio", 300, 0.5, validation_scheme_nonranking), 
  nonranking_infogain = filter_nonrankings(paths, output_prefix, "infogain", 300, 0.5, validation_scheme_nonranking), 
  nonranking_symuncert = filter_nonrankings(paths, output_prefix, "symuncert", 300, 0.5, validation_scheme_nonranking), 
  
  nonranking_FCBF = filter_nonrankings(paths, "FCBF", 300, 0.5, validation_scheme_nonranking)
)


cache <- new_cache("exp01_seqLen10_nSeq300_amylogram_encoding")
make(
  plan,
  parallelism = "future",
  jobs = 8,
  log_make = "exp01_seqLen10_nSeq300_amylogram_encoding.log",
  cache = cache,
  seed = 42
)

