# > readRDS("../data/alph20_prob_(19,26]/alph20_prob_(19,26]_rep_1_seqNum_6000_seqLen_80_nMotifs_1.Rds")
# A 6000x12209151 simple triplet matrix.
library(pbapply)
library(QuiPTsim)
library(parallel)

add_unigrams <- function(old_dir, new_dir, RDSnames, alphabet, num.cores=1) {
  
  mclapply(RDSnames, function(n) {
    m <- readRDS(paste0(old_dir, n))
    m_with_uni <- cbind(count_ngrams(attr(m, "sequences"), alphabet = alphabet, 1), m)
    
    for (a in c("sequences", "motifs", "masks", "target")) {
      attr(m_with_uni, a) <- attr(m, a)
    }
    new_path <- paste0(new_dir, n)
    saveRDS(object = m_with_uni, file = new_path)
  }, mc.cores = num.cores)
}

old_dir <- "../data/reduced_alph_enc_alph4_const/"
new_dir <- "../data/reduced_alph_enc_alph4_const_unigram/"
paths <- dir(old_dir)

df <- read.csv(paste0(old_dir, paths[grepl("csv", paths)]))

RDSnames <- sapply(strsplit(x = df[["path"]], split = "/"), function(x) x[[3]])


add_unigrams(old_dir, new_dir, RDSnames, letters[1:4])

# reduced_alph_enc_amylogram_encoding dataset broken

