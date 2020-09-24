library(drake)
library(QuiPTsim)

plan <- drake_plan(
  # data frame with paths of each RDS file containing n-gram counts
  df = read.csv("./reduced_alph_enc_alph4_const/alph4_const.csv"),

  # specified file paths
  paths = as.character(df[df$l_seq == 10 & df$n_motifs == 3, "path"])[1],

  # benchmark data for various methods
  data_QuiPT  = create_benchmark_data(paths, list(method="QuiPT")),
  data_Chi  = create_benchmark_data(paths, list(method="Chi-squared")),
  data_FCBF  = create_benchmark_data(paths, list(method="FCBF")),
  data_gainratio  = create_benchmark_data(paths, list(method="gainratio")),
  data_MRMR  = create_benchmark_data(paths, list(method="MRMR")),

  # results aggregating
  results_QuiPT = benchmark_summary(data_QuiPT,
                                              list(method = "QuiPT",
                                                   pval_thresholds = c(0.025, 0.01),
                                                   pval_adjustments = c("", "BH"))),
  results_FCBF = benchmark_summary(data_FCBF, list(method = "FCBF")),

  results_Chi = benchmark_summary(data_Chi,
                                            list(method = "Chi-squared",
                                                 pval_thresholds = c(0.025, 0.01),
                                                 pval_adjustments = c("", "BH"))),
  results_gainratio = benchmark_summary(data_gainratio,
                                           list(method = "gainratio",
                                                fraction = 0)),

  results = rbind(cbind(method_name = "QuiPT", results_QuiPT),
                  cbind(method_name = "Chi-squared", results_Chi),
                  cbind(method_name = "FCBF", pval_thresholds = "",
                        pval_adjustments = "", results_FCBF))
)

make(plan
     # HPC setup
     # ,parallelism = "clustermq",
    #jobs = 2,
     # console_log_file = "drake.log"
     )

vis_drake_graph(plan)



