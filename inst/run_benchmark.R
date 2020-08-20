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


