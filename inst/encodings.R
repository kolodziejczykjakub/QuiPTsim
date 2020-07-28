# probabilites from drake repo
weights <- readRDS("./inst/weights.Rds")

vis_dat_fit_best <- read.table("./inst/vis_dat_fit_best.csv", sep = ";", dec = ",", header = 1)

for (col in colnames(vis_dat_fit_best)) {
  vis_dat_fit_best[[col]] <- sub(",", ".", vis_dat_fit_best[[col]])
}
for (col in c("AUC_mean", "MCC_mean", "Sens_mean", "Spec_mean", "AUC_sd",
              "MCC_sd", "Sens_sd", "Spec_sd", "Cosine_similarity")) {
  vis_dat_fit_best[[col]] <- as.numeric(vis_dat_fit_best[[col]])
}


#add as a benchmark two encodings from the literature
aa1 = list(`1` = c("g", "a", "p", "v", "l", "i", "m"),
           `2` = c("k", "r", "h"),
           `3` = c("d", "e"),
           `4` = c("f", "w", "y", "s", "t", "c", "n", "q"))

aa2 = list(`1` = c("g", "a", "p", "v", "l", "i", "m", "f"),
           `2` = c("k", "r", "h"),
           `3` = c("d", "e"),
           `4` = c("s", "t", "c", "n", "q", "y", "w"))

library(AmyloGram)
amylogram_model_encoding <- AmyloGram_model[["enc"]]

amylogram_model_encoding <- list(`1` = "g",
                                 `2` = c("k", "p", "r"),
                                 `3` = c("i", "l", "v"),
                                 `4` = c("f", "w", "y"),
                                 `5` = c("a", "c", "h", "m"),
                                 `6` = c("d","e", "n", "q", "s", "t"))

#
PATH = "~/amylogram/AmyloGramAnalysis/"

library(seqinr)
source(paste0(PATH, "./functions/encode_amyloids.R"))

load(paste0(PATH, "./data/aa_groups.RData"))
aa_groups <- string2list(aa_groups)

raw_seqs_list <- c(read.fasta(paste0(PATH, "./data/amyloid_pos_full.fasta"),seqtype = "AA"),
                   read.fasta(paste0(PATH, "./data/amyloid_neg_full.fasta"),seqtype = "AA"))
#sequences longer than 5 aa and shorter than 26 aa
purified_seqs_id <- lengths(raw_seqs_list) > 5 & lengths(raw_seqs_list) < 26
seqs_list <- raw_seqs_list[purified_seqs_id]

seqs_m <- tolower(t(sapply(seqs_list, function(i)
  c(i, rep(NA, max(lengths(seqs_list)) - length(i))))))


# biogram::degenerate
# example


# weights
weights$positive$`[11,19]`

motifProbs <- weights$positive$`[11,19]`
seqProbs <- weights$negative$`[11,19]`

new_names_motif <-biogram::degenerate(tolower(names(weights$positive$`[11,19]`)), amylogram_model_encoding)
new_names_seq <-biogram::degenerate(tolower(names(weights$negative$`[11,19]`)), amylogram_model_encoding)

new_motifProbs <- unlist(lapply(1:6, function(x) sum(motifProbs[new_names_motif == x])))
new_seqProbs <- unlist(lapply(1:6, function(x) sum(seqProbs[new_names_seq == x])))

cosine_similarity(new_motifProbs, new_seqProbs)


#
table(c(seqs_m)) / sum(table(c(seqs_m)))
