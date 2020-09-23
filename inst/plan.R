library(drake)

plan <- drake_plan(
  # data frame with paths of each RDS file containing n-gram counts
  df = read.csv("./reduced_alph_enc_alph4_const/alph4_const.csv"),

  # specified file paths
  paths = as.character(df[df$l_seq == 10 & df$n_motifs == 3, "path"])[1],

  # benchmark data for various methods
  bm_data_QuiPT  = create_benchmark_data(paths, list(method="QuiPT")),
  bm_data_Chi  = create_benchmark_data(paths, list(method="Chi-squared")),
  bm_data_FCBF  = create_benchmark_data(paths, list(method="FCBF")),
  bm_data_gainratio  = create_benchmark_data(paths, list(method="gainratio")),
  bm_data_MRMR  = create_benchmark_data(paths, list(method="MRMR")),

  # results aggregating
  bm_results_QuiPT = benchmark_summary(bm_data_QuiPT,
                                              list(method = "QuiPT",
                                                   pval_thresholds = c(0.05, 0.01),
                                                   pval_adjustments = c("", "BH"))),
  bm_results_FCBF = benchmark_summary(bm_data_FCBF, list(method = "FCBF")),

  bm_results_Chi = benchmark_summary(bm_data_Chi,
                                            list(method = "Chi-squared",
                                                 pval_thresholds = c(0.05, 0.01),
                                                 pval_adjustments = c("", "BH"))),

  bm_results_QuiPT2 = benchmark_summary(bm_data_QuiPT,
                                              list(method = "QuiPT",
                                                   pval_thresholds = c(0.05, 0.01),
                                                   pval_adjustments = c("", "BH"))),
  results = rbind(cbind(method_name = "QuiPT", bm_results_QuiPT),
                   cbind(method_name = "Chi-squared", bm_results_Chi))

)


make(plan
     # HPC setup
     # ,parallelism = "clustermq",
     # jobs = 4,
     # console_log_file = "drake.log"
     )

vis_drake_graph(plan)
