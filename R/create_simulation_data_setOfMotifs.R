#' function generates set of RDS files containing n-gram matrices
#' @param replications number of matrices' generation for each set of parameters
#' @param seq_nums vector of sequences' numbers
#' @param seq_lengths vector of sequences' lengths
#' @param num_motifs vector of number of motifs to be injected
#' @param motif_set_size number of motifs to be created
#' @param alphabet elements used to build both sequence and motif
#' @param path file path of RDS files
#' @param title simulation name
#' @param save_files triggers if RDS will be saved
#' @param motifProbs alphabet probabilites for motifs
#' @param seqProbs alphabet probabilites for sequences
#' @param fraction fraction of positive sequences
#' @param n maximum number of alphabet elements in motif
#' @param d maximum number of gaps in motif
#' @return data frame with RDS files' details
#' @importFrom utils write.csv
#' @importFrom progress progress_bar
#' @importFrom tools md5sum
#' @export
create_simulation_data_set_of_motifs <- function(replications,
                                   seq_nums,
                                   seq_lengths,
                                   motif_set_size,
                                   num_motifs,
                                   alphabet,
                                   path,
                                   title,
                                   save_files = FALSE,
                                   motifProbs = NULL,
                                   seqProbs = NULL,
                                   fraction = 0.5,
                                   n = 4,
                                   d = 6) {

  if (!file.exists(path))
    stop("Output directory does not exist.")

  df <- data.frame(replication = integer(),
                   n_seq = integer(),
                   l_seq = integer(),
                   n_motifs = integer(),
                   n = integer(),
                   d = integer(),
                   seqProbs = character(),
                   motifProbs = character(),
                   path = character(),
                   md5sum = character())

  totalIterations = replications * length(seq_nums) * length(seq_lengths) *  length(num_motifs)
  pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)",
                         total = totalIterations)
  pb$tick(0)


  for (n_motifs in num_motifs) {
    for (n_seq in seq_nums) {
      for (l_seq in seq_lengths) {
        for (replication in 1:replications) {
     
          pb$tick(1)
          set.seed(replication)
          
          motifs <- lapply(1L:motif_set_size, 
                           function(dummy) generate_single_motif(alphabet, n, d, motifProbs))
          
          dat <- generate_sequences(n_seq,
                                    l_seq,
                                    alphabet,
                                    motifs,
                                    n_motifs,
                                    fraction = fraction,
                                    seqProbs = seqProbs,
                                    n = n,
                                    d = d)

          filePath = paste0(path, title,
                            "_rep_", replication,
                            "_seqNum_", n_seq,
                            "_seqLen_", l_seq,
                            "_nMotifs_", n_motifs,
                            ".Rds")

          if (save_files)
            saveRDS(dat, file = filePath)

          md5 <- ifelse(save_files, md5sum(filePath), "none")

          newdf <- data.frame(replication = replication,
                     n_seq = n_seq,
                     l_seq = l_seq,
                     n_motifs = n_motifs,
                     n = n,
                     d = d,
                     seqProbs = paste(motifProbs, collapse = ";"),
                     motifProbs = paste(motifProbs, collapse = ";"),
                     path = filePath,
                     md5sum = md5)

          df <- rbind(df, newdf)

        }
      }
    }
  }

  if (save_files)
    write.csv(df, file = paste0(path, title,".csv"))
  df
}

