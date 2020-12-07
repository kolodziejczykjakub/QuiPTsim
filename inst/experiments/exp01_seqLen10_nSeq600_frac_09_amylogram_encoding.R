library(QuiPTsim)
library(drake)

plan <- drake_plan(
  
  ###### data frame containing datasets' details and their paths to RDS files
  df = read.csv("~/projects/QuiPTsim-data/reduced_alph_enc_amylogram_encoding_unigram/amylogram_encoding.csv"),
  
  numSeq = 600,
  fraction = 0.9,
  ###### selected paths
  paths = paste0("~/projects/QuiPTsim-data/reduced_alph_enc_amylogram_encoding_unigram/",
                 sapply(strsplit(x = df[df$l_seq==10 & df$n_motifs==1,"path"], split = "/"), function(x) x[[3]])),
  
  output_prefix = "~/experiment-results/exp01-seqLen10-nSeq600-amylogram-encoding/result_",
  
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
                           n_kmers = c(2:128, 2^(8:12)),
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
  
  ranking_QuiPT = filter_rankings(paths, output_prefix, "QuiPT", numSeq, fraction, validation_scheme),
  ranking_Chi = filter_rankings(paths, output_prefix, "Chi-squared", numSeq, fraction, validation_scheme),
  ranking_FCBF = filter_rankings(paths, output_prefix, "FCBF", numSeq, fraction, validation_scheme),
  ranking_infogain = filter_rankings(paths, output_prefix, "infogain", numSeq, fraction, validation_scheme),
  ranking_gainratio = filter_rankings(paths, output_prefix, "gainratio", numSeq, fraction, validation_scheme),
  ranking_symuncert = filter_rankings(paths, output_prefix, "symuncert", numSeq, fraction, validation_scheme),
  ranking_MRMR = filter_rankings(paths, output_prefix, "MRMR", numSeq, fraction, validation_scheme),
  ranking_JMI = filter_rankings(paths, output_prefix, "JMI", numSeq, fraction, validation_scheme),
  ranking_JMIM = filter_rankings(paths, output_prefix, "JMIM", numSeq, fraction, validation_scheme),
  ranking_DISR = filter_rankings(paths, output_prefix, "DISR", numSeq, fraction, validation_scheme),
  ranking_NJMIM = filter_rankings(paths, output_prefix, "NJMIM", numSeq, fraction, validation_scheme),
  
  
  thresholds = c(0.01, 0.05),
  
  nonranking_QuiPT = filter_nonrankings(paths, output_prefix, "QuiPT", numSeq, fraction, validation_scheme_nonranking,
                                        thresholds),
  nonranking_Chi = filter_nonrankings(paths, output_prefix, "Chi-squared", numSeq, fraction, validation_scheme_nonranking,
                                      thresholds),
  
  nonranking_gainratio = filter_nonrankings(paths, output_prefix, "gainratio", numSeq, fraction, validation_scheme_nonranking), 
  nonranking_infogain = filter_nonrankings(paths, output_prefix, "infogain", numSeq, fraction, validation_scheme_nonranking), 
  nonranking_symuncert = filter_nonrankings(paths, output_prefix, "symuncert", numSeq, fraction, validation_scheme_nonranking), 
  
  nonranking_FCBF = filter_nonrankings(paths, output_prefix, "FCBF", numSeq, fraction, validation_scheme_nonranking)
)


cache <- new_cache("exp01_seqLen10_nSeq600_amylogram_encoding")
make(
  plan,
  parallelism = "future",
  jobs = 8,
  log_make = "exp01_seqLen10_nSeq600_amylogram_encoding.log",
  cache = cache,
  seed = 42
)

