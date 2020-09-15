# > readRDS("../data/alph20_prob_(19,26]/alph20_prob_(19,26]_rep_1_seqNum_6000_seqLen_80_nMotifs_1.Rds")
# A 6000x12209151 simple triplet matrix.
library(pbapply)

add_unigrams <- function(paths, alphabet) {
  
  pblapply(paths, function(path) {
    m <- readRDS(path)
    m_with_uni <- cbind(count_ngrams(attr(m, "sequences"), alphabet = alphabet, 1), m)
    
    for (a in c("sequences", "motifs", "masks", "target")) {
      attr(m_with_uni, a) <- attr(m, a)
    }
    
    saveRDS(object = m_with_uni, file = path)
  })
}

path <- "../data/reduced_alph_enc_amylogram_encoding/"
paths <- dir("../data/reduced_alph_enc_amylogram_encoding/")
df <- read.csv(paste0(path, paths[grepl("csv", paths)]))
ngram_matrices_paths <-  paste0("../data/", df[["path"]])

add_unigrams(ngram_matrices_paths, letters[1:6])
