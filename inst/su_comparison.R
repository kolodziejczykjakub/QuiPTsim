library(QuiPTsim)
library(FSelectorRcpp)

calculate_ngram_matrix_su <- function(m){
  
  x <- data.frame(as.matrix(m))
  y <- attr(m, "target")
  
  score <- information_gain(x=x, y=y, discIntegers=F, type="symuncert")
  score$importance
}

### Read exp1 data

exp1_df <- read.csv("reduced_alph_enc_amylogram_encoding/amylogram_encoding.csv")

exp1_df_paths_nMotifs1 <- exp1_df[exp1_df$n_motifs == 1, "path"]
exp1_df_paths_nMotifs2 <- exp1_df[exp1_df$n_motifs == 2, "path"]

SU_exp1_nSeq300_nMotifs1 <- lapply(exp1_df_paths_nMotifs1, function(path) {
  m <- read_ngram_matrix(path, n=300, fraction=0.5)
  calculate_ngram_matrix_su(m)
})

SU_exp1_nSeq300_nMotifs2 <- lapply(exp1_df_paths_nMotifs2, function(path) {
  m <- read_ngram_matrix(path, n=300, fraction=0.5)
  calculate_ngram_matrix_su(m)
})


SU_exp1_nSeq600_nMotifs1 <- lapply(exp1_df_paths_nMotifs1, function(path) {
  m <- read_ngram_matrix(path, n=600, fraction=0.5)
  calculate_ngram_matrix_su(m)
})

SU_exp1_nSeq600_nMotifs2 <- lapply(exp1_df_paths_nMotifs2, function(path) {
  m <- read_ngram_matrix(path, n=600, fraction=0.5)
  calculate_ngram_matrix_su(m)
})

### Read exp3 data

exp3_setMotifs5_df <- read.csv("./exp3_reduced_alph_enc_5motifs_amylogram_encoding/amylogram_encoding.csv")
exp3_setMotifs15_df <- read.csv("./exp3_reduced_alph_enc_15motifs_amylogram_encoding/amylogram_encoding.csv")

exp3_setMotifs5_df_nMotifs1_paths <- exp3_setMotifs5_df[exp3_setMotifs5_df$n_motifs==1, "path"]
exp3_setMotifs15_df_nMotifs1_paths <- exp3_setMotifs15_df[exp3_setMotifs15_df$n_motifs==1, "path"]


exp3_setMotifs5_df_nMotifs2_paths <- exp3_setMotifs5_df[exp3_setMotifs5_df$n_motifs==2, "path"]
exp3_setMotifs15_df_nMotifs2_paths <- exp3_setMotifs15_df[exp3_setMotifs15_df$n_motifs==2, "path"]

SU_exp3_nSeq300_setMotifs5_nMotifs1 <- lapply(exp3_setMotifs5_df_nMotifs1_paths, function(path) {
  m <- read_ngram_matrix(path, n=300, fraction=0.5)
  calculate_ngram_matrix_su(m)
})

SU_exp3_nSeq600_setMotifs5_nMotifs1 <- lapply(exp3_setMotifs5_df_nMotifs1_paths, function(path) {
  m <- read_ngram_matrix(path, n=600, fraction=0.5)
  calculate_ngram_matrix_su(m)
})

SU_exp3_nSeq300_setMotifs5_nMotifs2 <- lapply(exp3_setMotifs5_df_nMotifs2_paths, function(path) {
  m <- read_ngram_matrix(path, n=300, fraction=0.5)
  calculate_ngram_matrix_su(m)
})

SU_exp3_nSeq600_setMotifs5_nMotifs2 <- lapply(exp3_setMotifs5_df_nMotifs2_paths, function(path) {
  m <- read_ngram_matrix(path, n=600, fraction=0.5)
  calculate_ngram_matrix_su(m)
})

###

SU_exp3_nSeq300_setMotifs15_nMotifs1 <- lapply(exp3_setMotifs15_df_nMotifs1_paths, function(path) {
  m <- read_ngram_matrix(path, n=300, fraction=0.5)
  calculate_ngram_matrix_su(m)
})

SU_exp3_nSeq600_setMotifs15_nMotifs1 <- lapply(exp3_setMotifs15_df_nMotifs1_paths, function(path) {
  m <- read_ngram_matrix(path, n=600, fraction=0.5)
  calculate_ngram_matrix_su(m)
})

SU_exp3_nSeq300_setMotifs15_nMotifs2 <- lapply(exp3_setMotifs15_df_nMotifs2_paths, function(path) {
  m <- read_ngram_matrix(path, n=300, fraction=0.5)
  calculate_ngram_matrix_su(m)
})

SU_exp3_nSeq600_setMotifs15_nMotifs2 <- lapply(exp3_setMotifs15_df_nMotifs2_paths, function(path) {
  m <- read_ngram_matrix(path, n=600, fraction=0.5)
  calculate_ngram_matrix_su(m)
})

### Read data exp2
SU_exp2_nSeq300_nMotifs1 <- lapply(1:10, function(dummy) {
  m <- do.call(rbind_ngram_matrices, lapply(sample(exp1_df_paths_nMotifs1, 2), function(x)
    read_ngram_matrix(x, n=150, fraction=0.5)))
  calculate_ngram_matrix_su(m)
})

SU_exp2_nSeq600_nMotifs1 <- lapply(1:10, function(dummy) {
  m <- do.call(rbind_ngram_matrices, lapply(sample(exp1_df_paths_nMotifs1, 2), function(x)
    read_ngram_matrix(x, n=300, fraction=0.5)))
  calculate_ngram_matrix_su(m)
})

SU_exp2_nSeq300_nMotifs2 <- lapply(1:10, function(dummy) {
  m <- do.call(rbind_ngram_matrices, lapply(sample(exp1_df_paths_nMotifs2, 2), function(x)
    read_ngram_matrix(x, n=150, fraction=0.5)))
  calculate_ngram_matrix_su(m)
})

SU_exp2_nSeq600_nMotifs2 <- lapply(1:10, function(dummy) {
  m <- do.call(rbind_ngram_matrices, lapply(sample(exp1_df_paths_nMotifs2, 2), function(x)
    read_ngram_matrix(x, n=300, fraction=0.5)))
  calculate_ngram_matrix_su(m)
})

## aggregate functions

calculate_top_SU <- function(results, n=4096) {
  apply(do.call(cbind, lapply(results, function(x){
    sort(x, decreasing = T)[1:4096]
  })), 1, mean)
}

calculate_valuable_kmers <- function(results, threshold=0.025) {
  sapply(results, function(x) sum(x > threshold))
}

calculate_max_su <- function(results) {
  sapply(results, max)
}

#########################
SU_results <- list(
  SU_exp1_nSeq300_nMotifs1,
  SU_exp1_nSeq600_nMotifs1,
  SU_exp1_nSeq300_nMotifs2,
  SU_exp1_nSeq600_nMotifs2,
  
  SU_exp2_nSeq300_nMotifs1,
  SU_exp2_nSeq600_nMotifs1,
  SU_exp2_nSeq300_nMotifs2,
  SU_exp2_nSeq600_nMotifs2,
  
  SU_exp3_nSeq300_setMotifs5_nMotifs1, 
  SU_exp3_nSeq600_setMotifs5_nMotifs1,
  SU_exp3_nSeq300_setMotifs5_nMotifs2,
  SU_exp3_nSeq600_setMotifs5_nMotifs2,
  
  SU_exp3_nSeq300_setMotifs15_nMotifs1,
  SU_exp3_nSeq600_setMotifs15_nMotifs1,
  SU_exp3_nSeq300_setMotifs15_nMotifs2,
  SU_exp3_nSeq600_setMotifs15_nMotifs2
)

results <- do.call(rbind, lapply(SU_results, function(x) {
  data.frame(name=deparse(substitute(x)),
             n_kmers=calculate_valuable_kmers(x),
             max_su=calculate_max_su(x)
  )
}))

SU_names <- c(
  "E1-300s-1m",
  "E1-600s-1m",
  "E1-300s-2m",
  "E1-600s-2m",
  
  "E2-300s-1m",
  "E2-600s-1m",
  "E2-300s-2m",
  "E2-600s-2m",
  
  "E3-300s-5s-1m",
  "E3-600s-5s-1m",
  "E3-300s-5s-2m",
  "E3-600s-5s-2m",
  
  "E3-300s-15s-1m",
  "E3-600s-15s-1m",
  "E3-300s-15s-2m",
  "E3-600s-15s-2m"
  )

results$Name <- rep(SU_names, each=10)
results$Sequences <- as.factor(rep(rep(c(300, 600), 8), each=10))
results$Motifs <- as.factor(rep(rep(1:2, 4,each=2), each=10))
results$Experiment <- rep(rep(c("Experiment 1", "Experiment 2", "Experiment 3 - 5 motifs", "Experiment 4 - 15 motifs"), each=4),
                          each=10)
### gg plotly phase
library(ggplot2)

library(plotly)
ggplot(results, aes(x=n_kmers, y=max_su, shape=Motifs, color=Sequences)) +
  geom_point(size=2) + 
  facet_wrap(vars(Experiment)) +
  scale_shape_manual(values=c(21, 24)) +
  theme_bw() +
  scale_colour_brewer(palette = "Set1")
#ggplotly()


su_curves <- do.call(rbind, lapply(SU_results, function(x) {
  data.frame(name=deparse(substitute(x)),
             n_kmers=1:4096,
             max_su=calculate_top_SU(x, n = 4096)
  )
}))

su_curves$Name <- rep(SU_names, each=4096)
su_curves$Sequences <- as.factor(rep(rep(c(300, 600), 8), each=4096))
su_curves$Motifs <- as.factor(rep(rep(1:2, 4,each=2), each=4096))
su_curves$Experiment <- rep(rep(c("Experiment 1", "Experiment 2", "Experiment 3 - 5 motifs", "Experiment 3 - 15 motifs"), each=4),
                          each=4096)
su_curves$MotifsSequences <- paste(su_curves$Motifs, su_curves$Sequences)


ggplot(su_curves, aes(x=n_kmers, y=max_su, color=Motifs, linetype=Sequences)) +
  facet_wrap(vars(Experiment)) +
  geom_line() +
  theme_bw() +
  scale_colour_brewer(palette = "Set1")
#ggplotly()
