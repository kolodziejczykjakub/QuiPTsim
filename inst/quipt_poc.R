results <- read.csv("./inst/example_quipt_results.csv")

# True positive
# both should be considered as positive?
true_positives <- results[(results$contains.motif | results$motif.part) & (results$p.value <= 0.05), ]

# True negative
true_negatives <- results[!(results$contains.motif | results$motif.part) & !(results$p.value <= 0.05), ]

# False positive
false_positives <- results[!(results$contains.motif | results$motif.part) & (results$p.value <= 0.05), ]

# False negative
false_negatives <- results[(results$contains.motif | results$motif.part) & !(results$p.value <= 0.05), ]
