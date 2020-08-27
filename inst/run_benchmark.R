library(QuiPTsim)

df <- read.csv("./reduced_alph_enc_alph4_const/alph4_const.csv")

paths <- as.character(df[df$l_seq == 10 & df$n_motifs == 3, "path"])

bm <- create_benchmark_data(paths[1:2], list(method = "QuiPT"))

benchmark_summary(bm, list(method = "QuiPT",
                           pval_thresholds = c(0.05, 0.01),
                           pval_adjustments = c("", "BH")))

#########################################################################

bm2 <- create_benchmark_data(paths[1:2], list(method = "QuiPT",
                                        shuffle_matrices = 2,
                                        fraction = 0.5,
                                        n = 300))

benchmark_summary(bm2, list(method = "QuiPT",
                           pval_thresholds = c(0.05, 0.01),
                           pval_adjustments = c("", "BH")))





paths1 <- as.character(df[df$l_seq == 10 & df$n_motifs == 1, "path"])
paths2 <- as.character(df[df$l_seq == 10 & df$n_motifs == 2, "path"])

c(paths1, paths2)

# TODO:
list(paths1, paths2)

####################################################################
setup = list(method = "QuiPT",
             fraction = 0.5,
             n = 100,
             pval_thresholds = c(0.05, 0.01),
             pval_adjustments = c("", "BH"))

results <- QuiPTsimBenchmark(paths, setup)

#########################################################################

bm_FCBF <- create_benchmark_data(paths, list(method = "FCBF",
                                              fraction = 0.5,
                                              n = 300))
bm_QuiPT <- create_benchmark_data(paths, list(method = "QuiPT",
                                             fraction = 0.5,
                                             n = 300))

benchmark_summary(bm_FCBF, list(method = "FCBF"))
benchmark_summary(bm_QuiPT, list(method = "QuiPT",
                            pval_thresholds = c(0.05, 0.01),
                            pval_adjustments = c("", "BH")))



