library(QuiPTsim)

df <- read.csv("./reduced_alph_enc_alph4_const/alph4_const.csv")

paths <- as.character(df[df$l_seq == 10 & df$n_motifs == 1, "path"])

bm <- create_benchmark_data(paths[1:3], list(method = "QuiPT"))
bm

# Boruta(as.matrix(m), attr(m, "target"))

bm[[1]]$positive.ngram


benchmark_summary(bm, list(method = "QuiPT", pval_threshold = 0.05))


