library(QuiPTsim)

paths <- c("./reduced_alph_enc_alph4_const/alph4_const_rep_9_seqNum_600_seqLen_80_nMotifs_3.Rds")
pval_thresholds <- c(0.01, 0.05)
pval_adjustments <- c("", "BH", "holm")

QuiPT_simulation_summary <- function(paths, pval_thresholds, pval_adjustments) {

  for (path in paths) {
    for (pval_th in pval_thresholds) {
      for (pval_adj in pval_adjustments) {

        m <- read_ngram_matrix(path)
        results <- QuiPT_summary(m)

        y_true <- results[["positive.ngram"]]

        if (pval_adj ==  "") {
          y_pred <- (results[["p.value"]] < pval_th)
        } else {
          y_pred <- p.adjust(results[["p.value"]] < pval_th, method = pval_adj)
        }

        print(sensitivity(y_true, y_pred))
        print(specificity(y_true, y_pred))

      }
    }
  }
}
