
#' function generates set of RDS files containing n-gram matrices
#' @param replications number of matrices' generation for each set of parameters
#' @param seq_nums vector of sequences' numbers
#' @param seq_lengths vector of sequences' lengths
#' @param n_motifs vector of number of motifs to be injected
#' @param alphabet elements used to build both sequence and motif
#' @param path file path of RDS files
#' @param title simulation name
#' @param save_files triggers if RDS will be saved
#' @return data frame with RDS files' details
#' @importFrom pbapply pblapply
#' @importFrom utils write.csv
#' @export
#' @examples
#' alph <- 1L:4
#' reps <- 1
#' n_seq <- c(25)
#' l_seq <- c(8)
#' n_motifs <- c(1, 2)
#' path = "./"
#' create_simulation_data(reps, n_seq, l_seq, n_motifs, alph, path, "SEQ",FALSE)

create_simulation_data <- function(replications,
                                   seq_nums,
                                   seq_lengths,
                                   n_motifs,
                                   alphabet,
                                   path,
                                   title,
                                   save_files = FALSE) {

  if (!file.exists(path))
    stop("Output directory does not exist.")

  sim_quipt <- pblapply(1L:replications, function(replication) {
    lapply(seq_nums, function(n_seq) {
      lapply(seq_lengths, function(l_seq) {
        lapply(n_motifs, function(n_motifs){

          set.seed(replication)
          motifs <- generate_motifs(alphabet, n_motifs)
          dat <- generate_sequences(n_seq, l_seq, alphabet, motifs, n_motifs)


          filePath = paste0(path, title,
                            "_rep_", replication,
                            "_seqNum_", n_seq,
                            "_seqLen_", l_seq,
                            "_nMotifs_", n_motifs,
                            ".Rds")

          if (save_files)
            saveRDS(dat, file = filePath)

          data.frame(replication = replication,
                     n_seq = n_seq,
                     l_seq = l_seq,
                     n_motifs = n_motifs,
                     path = filePath
          )
        }) %>%
          do.call(rbind, .)
      }) %>%
        do.call(rbind, .)
    }) %>%
      do.call(rbind, .)
  }) %>%
    do.call(rbind, .)

  if (save_files)
    write.csv(sim_quipt, file = paste0(path, title,".csv"))

  sim_quipt

}


