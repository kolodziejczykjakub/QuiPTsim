results <- read.csv("./inst/example_quipt_results.csv")


m1 <- read_ngram_matrix("./reduced_alph_enc_alph4_const/alph4_const_rep_1_seqNum_600_seqLen_10_nMotifs_3.Rds")
m2 <- read_ngram_matrix("./reduced_alph_enc_alph4_const/alph4_const_rep_2_seqNum_600_seqLen_10_nMotifs_3.Rds")

m3 <- rbind_ngram_matrices(m1, m2)

metrics <- function(results) {

  # True positive
  # both should be considered as positive?
  true_positives <- sum((results$contains.motif | results$motif.part) & (results$p.value <= 0.05))

  # True negative
  true_negatives <- sum(!(results$contains.motif | results$motif.part) & !(results$p.value <= 0.05))

  # False positive
  false_positives <- sum(!(results$contains.motif | results$motif.part) & (results$p.value <= 0.05))

  # False negative
  false_negatives <- sum((results$contains.motif | results$motif.part) & !(results$p.value <= 0.05))

  cat(paste0("TPR(sens):\n", true_positives / (true_positives + false_negatives)), "\n")
  cat(paste0("TNR(spec):\n", true_negatives / (true_negatives + false_positives)), "\n")
}

metrics2 <- function(results) {

  # True positive
  # both should be considered as positive?
  true_positives <- sum((results$motif.part) & (results$p.value <= 0.05))

  # True negative
  true_negatives <- sum(!(results$motif.part) & !(results$p.value <= 0.05))

  # False positive
  false_positives <- sum(!(results$motif.part) & (results$p.value <= 0.05))

  # False negative
  false_negatives <- sum((results$motif.part) & !(results$p.value <= 0.05))

  cat(paste0("TPR(sens):\n", true_positives / (true_positives + false_negatives)), "\n")
  cat(paste0("TNR(spec):\n", true_negatives / (true_negatives + false_positives)), "\n")
}

results1 <- QuiPT_summary(m1)
results2 <- QuiPT_summary(m2)
results3 <- QuiPT_summary(m3)
metrics(results1)
metrics(results2)
metrics(results3)


metrics2(results1)
metrics2(results2)
metrics2(results3)
