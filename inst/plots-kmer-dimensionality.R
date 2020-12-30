alphabetSize <- 4
kmerLengths <- 1:4
contiguous_kmers <- alphabetSize ^ kmerLengths

motifLength <- 1
d <- 6
gaps <- expand.grid(list(0:d)[rep(1, motifLength - 1)])
possibleGaps <- gaps[apply(gaps, 1, sum) <= d, , drop = FALSE]

count_gapped_kmers <- function(alphabetSize, kmerLengths, d) {
  sapply(kmerLengths, function(kmerLength) {

    if (kmerLength == 1){
      x <-alphabetSize
    } else {
      gaps <- expand.grid(list(0:d)[rep(1, kmerLength - 1)])
      possibleGaps <- gaps[apply(gaps, 1, sum) <= d, , drop = FALSE]
      x <- (alphabetSize ^ kmerLength) * nrow(possibleGaps)
    }
    x
  })
}

library(ggplot2)
library(scales)

# total gaps
kmers <- lapply(0:d, function(x) count_gapped_kmers(alphabetSize, kmerLengths, x))
kmer_matrix <- do.call(cbind, kmers)
colnames(kmer_matrix) <- 0:d
data <- reshape2::melt(kmer_matrix)
colnames(data) <- c("kmerLength", "TotalGaps", "nKmers")
data$TotalGaps <- as.factor(data$TotalGaps)
ggplot(data, aes(x=kmerLength, y = nKmers, colour = TotalGaps)) +
  geom_point(size=2.3) +
  scale_y_log10(labels=trans_format('log10',math_format(10^.x))) +
  theme_bw()

# alphabet size
alphabetSizes <- c(4,  6, 20)
kmerLengths <- 1:4
kmers <- lapply(alphabetSizes, function(x) count_gapped_kmers(x, kmerLengths, 0))
kmer_matrix <- do.call(cbind, kmers)
colnames(kmer_matrix) <- alphabetSizes
data <- reshape2::melt(kmer_matrix)
colnames(data) <- c("kmerLength", "AlphabetSize", "nKmers")
data$AlphabetSize <- as.factor(data$AlphabetSize)
ggplot(data, aes(x=kmerLength, y = nKmers, colour = AlphabetSize)) +
  geom_point(size=2.3) +
  scale_y_log10(labels=trans_format('log10',math_format(10^.x))) +
  theme_bw()
