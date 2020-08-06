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
##############################

m1_3 <- read_ngram_matrix("./reduced_alph_enc_alph4_const/alph4_const_rep_1_seqNum_600_seqLen_10_nMotifs_3.Rds")
m2_3 <- read_ngram_matrix("./reduced_alph_enc_alph4_const/alph4_const_rep_2_seqNum_600_seqLen_10_nMotifs_3.Rds")
m3_3 <- read_ngram_matrix("./reduced_alph_enc_alph4_const/alph4_const_rep_3_seqNum_600_seqLen_10_nMotifs_3.Rds")
m4_3 <- read_ngram_matrix("./reduced_alph_enc_alph4_const/alph4_const_rep_4_seqNum_600_seqLen_10_nMotifs_3.Rds")

m_3 <- rbind_ngram_matrices(m3_3, m4_3)

matrices_3 <- list(m1_3, m2_3, m3_3, m4_3, m_3)
results_3 <- lapply(matrices_3, QuiPT_summary)
print("N-grams containg motifs are positve")
lapply(results_3, metrics) -> tmp
print("N-grams containg motifs are negative")
lapply(results_3, metrics2) -> tmp

##############################

m1_1 <- read_ngram_matrix("./reduced_alph_enc_alph4_const/alph4_const_rep_1_seqNum_600_seqLen_10_nMotifs_1.Rds")
m2_1 <- read_ngram_matrix("./reduced_alph_enc_alph4_const/alph4_const_rep_2_seqNum_600_seqLen_10_nMotifs_1.Rds")
m3_1 <- read_ngram_matrix("./reduced_alph_enc_alph4_const/alph4_const_rep_3_seqNum_600_seqLen_10_nMotifs_1.Rds")
m4_1 <- read_ngram_matrix("./reduced_alph_enc_alph4_const/alph4_const_rep_4_seqNum_600_seqLen_10_nMotifs_1.Rds")

m_1 <- rbind_ngram_matrices(m3_1, m4_1)

matrices_1 <- list(m1_1, m2_1, m3_1, m4_1, m_1)
results_1 <- lapply(matrices_1, QuiPT_summary)

print("N-grams containg motifs are positve")
lapply(results_1, metrics) -> tmp
print("N-grams containg motifs are negative")
lapply(results_1, metrics2) -> tmp

