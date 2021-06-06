library(QuiPTsim)


library(dplyr)
library(AmyloGram)
library(biogram)

data("AmyloGram_model")

AmyloGram_alph <- lapply(AmyloGram_model[["enc"]], toupper)

prepare_sequences <- function(path){
  raw_seqs <- read_fasta(path)
  purified_seqs <- lengths(raw_seqs) > 5 & lengths(raw_seqs) < 26
  seqs <- raw_seqs[purified_seqs]
  
  lapply(seqs, degenerate, AmyloGram_alph)
}

pos_sequences <- prepare_sequences("https://raw.githubusercontent.com/michbur/AmyloGramAnalysis/master/data/amyloid_pos_full.fasta")
neg_sequences <- prepare_sequences("https://raw.githubusercontent.com/michbur/AmyloGramAnalysis/master/data/amyloid_neg_full.fasta")

length(pos_sequences)
length(neg_sequences)

paste0("Number of positive sequences: ", length(pos_sequences))
paste0("Number of negative sequences: ", length(neg_sequences))
paste0("Number of sequences: ", length(pos_sequences) + length(neg_sequences))
fraction = length(pos_sequences) / (length(pos_sequences) + length(neg_sequences))
paste0("Fraction of positive sequences: ", round(fraction, 2))

target <- c(rep(TRUE, length(pos_sequences)), rep(FALSE, length(neg_sequences)))

sequences <- lapply(c(pos_sequences, neg_sequences), function(x) paste0(x, collapse=""))

kmers <- QuiPTsim::count_ngrams(sequences, alphabet = as.character(1:6), n=4, d=6)
attr(kmers, "target") <- target

hist(lengths(c(pos_sequences, neg_sequences), use.names=F))
saveRDS(kmers, "./inst/AmylogramExp.Rds")

