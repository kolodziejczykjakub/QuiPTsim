#' function generates set of RDS files containing n-gram matrices
#' @param replications number of matrices' generation for each set of parameters
#' @param seq_nums vector of sequences' numbers
#' @param seq_lengths vector of sequences' lengths
#' @param num_motifs vector of number of motifs to be injected
#' @param alphabet elements used to build both sequence and motif
#' @param path file path of RDS files
#' @param title simulation name
#' @param save_files triggers if RDS will be saved
#' @return data frame with RDS files' details
#' @importFrom utils write.csv
#' @importFrom progress progress_bar
#' @export
#' @examples
#' alph <- 1L:4
#' reps <- 10
#' n_seq <- c(25)
#' l_seq <- c(8)
#' n_motifs <- c(1, 2)
#' path = "./"
#' create_simulation_data(reps, n_seq, l_seq, n_motifs, alph, path, "SEQ",FALSE)

create_simulation_data <- function(replications,
                                   seq_nums,
                                   seq_lengths,
                                   num_motifs,
                                   alphabet,
                                   path,
                                   title,
                                   save_files = FALSE) {

  if (!file.exists(path))
    stop("Output directory does not exist.")

  df <- data.frame(replication = integer(),
                   n_seq = integer(),
                   l_seq = integer(),
                   n_motifs = integer(),
                   path = character())

  totalIterations = replications * length(seq_nums) * length(seq_lengths) *  length(num_motifs)
  pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)",
                         total = totalIterations)
  pb$tick(0)

  for (replication in 1:replications) {
    for (n_seq in seq_nums) {
      for (l_seq in seq_lengths) {
        for (n_motifs in num_motifs) {

          pb$tick(1)
          set.seed(replication)

          dat <- generate_sequences(n_seq, l_seq, alphabet, motifs, n_motifs)

          filePath = paste0(path, title,
                            "_rep_", replication,
                            "_seqNum_", n_seq,
                            "_seqLen_", l_seq,
                            "_nMotifs_", n_motifs,
                            ".Rds")

          if (save_files)
            saveRDS(dat, file = filePath)

          newdf <- data.frame(replication = replication,
                     n_seq = n_seq,
                     l_seq = l_seq,
                     n_motifs = n_motifs,
                     path = filePath)

          df <- rbind(df, newdf)

        }
      }
    }
  }

  if (save_files)
    write.csv(df, file = paste0(path, title,".csv"))

  df
}


alph <- 1L:4
reps <- 10
n_seq <- c(25)
l_seq <- c(8)
n_motifs <- c(1, 2)
path = "./results/"
create_simulation_data(reps, n_seq, l_seq, n_motifs, alph, path, "SEQ",TRUE)
