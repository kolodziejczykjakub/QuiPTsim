library(QuiPTsim)

df <- read.csv("./reduced_alph_enc_alph4_const/alph4_const.csv")
paths <- as.character(df[df$l_seq == 10 & df$n_motifs == 3, "path"])
############
bm <- create_benchmark_data(paths[1:2], list(method = "QuiPT"))
benchmark_summary(bm, list(method = "QuiPT",
                           pval_thresholds = c(0.05, 0.01),
                           pval_adjustments = c("", "BH")))
############
bm2 <- create_benchmark_data(paths[1:2], list(method = "QuiPT",
                                        shuffle_matrices = 2,
                                        fraction = 0.5,
                                        n = 30))
benchmark_summary(bm2, list(method = "QuiPT",
                           pval_thresholds = c(0.05, 0.01),
                           pval_adjustments = c("", "BH")))
############
setup = list(method = "QuiPT",
             fraction = 0.5,
             n = 10,
             pval_thresholds = c(0.05, 0.01),
             pval_adjustments = c("", "BH"))
results <- QuiPTsimBenchmark(paths[3:4], setup)
############
bm_FCBF <- create_benchmark_data(paths[1:2], list(method="FCBF"))
benchmark_summary(bm_FCBF, list(method = "FCBF"))
############
bm_chi <- create_benchmark_data(paths[1:2], list(method = "Chi-squared",
                                             fraction = 0.5,
                                             n = 30))
benchmark_summary(bm_chi, list(method = "Chi-squared",
                                 pval_thresholds = c(0.05, 0.01),
                                 pval_adjustments = c("", "BH")))
############
bm_fselector <- create_benchmark_data(paths[1:2], list(method = "FSelectorRcpp",
                                                 fraction = 0.5,
                                                 n = 10))
benchmark_summary(bm_fselector, list(method = "FSelectorRcpp",
                                     fractions = c(0, 0.001, 0.01, 0.05)))
############
bm_praznik <- create_benchmark_data(paths[1:2], list(method = "MRMR",
                                                       fraction = 0.5,
                                                       n = 30))
benchmark_summary(bm_praznik, list(method = "MRMR"))



