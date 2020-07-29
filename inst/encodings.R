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

raw_seqs_positive <- read.fasta(paste0(PATH, "./data/amyloid_pos_full.fasta"),seqtype = "AA")
raw_seqs_negative <- read.fasta(paste0(PATH, "./data/amyloid_neg_full.fasta"),seqtype = "AA")

seqs_list_positive <- raw_seqs_positive[lengths(raw_seqs_positive) > 5 & lengths(raw_seqs_positive) < 26]
seqs_list_negative <- raw_seqs_negative[lengths(raw_seqs_negative) > 5 & lengths(raw_seqs_negative) < 26]

seqs_m_pos <- tolower(t(sapply(seqs_list_positive, function(i)
  c(i, rep(NA, max(lengths(seqs_list_positive)) - length(i))))))

seqs_m_neg <- tolower(t(sapply(seqs_list_negative, function(i)
  c(i, rep(NA, max(lengths(seqs_list_negative)) - length(i))))))

occ_positive <- table(seqs_m_pos)
occ_negative <- table(seqs_m_neg)

# delete "-" from occurences
occ_positive <- occ_positive[2:length(occ_positive)]

probs_positive <- occ_positive / sum(occ_positive)
probs_negative <- occ_negative / sum(occ_negative)

# biogram::degenerate
create_encoded_probabilites <- function(encoding, probs_positive, probs_negative) {

  new_names_pos <- biogram::degenerate(names(probs_positive), encoding)
  new_names_neg <- biogram::degenerate(names(probs_negative), encoding)

  new_motifProbs <- unlist(lapply(1:length(encoding),
                                  function(x) sum(probs_positive[new_names_pos == x])))
  new_seqProbs <- unlist(lapply(1:length(encoding),
                                function(x) sum(probs_negative[new_names_neg == x])))

  list(positive = new_motifProbs,
       negative = new_seqProbs)

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


encodingProbs <- list(amylogram_encoding = create_encoded_probabilites(amylogram_model_encoding,
                                                                       probs_positive,
                                                                       probs_negative),
                      aa1 = create_encoded_probabilites(aa1,
                                                        probs_positive,
                                                        probs_negative),
                      aa2 = create_encoded_probabilites(aa2,
                                                        probs_positive,
                                                        probs_negative))


library(QuiPTsim)
lapply(encodingProbs, function(x) cosine_similarity(x$positive, x$negative))

saveRDS(encodingProbs, "./inst/encodingProbs.Rds")
