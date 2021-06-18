# Setup ###############################################################################################
experiment_name <- "Experiment 3 (1 motif injected, 15 random motifs, 600 sequences, unbalanced datasets)"
exp_prefix <- "exp3_set15_1m_300s_unbalanced"

data_paths <- paste0("../data/experiment-results/exp-03/exp3_reduced_alph_enc_15motifs_amylogram_encoding/exp03-nSeq300-nMotifs1",
                     c("-frac01/", "-frac025/", "-frac075/", "-frac09/"))
cache_paths <- paste0("../data/experiment-drake-caches/exp03-15motifs/exp03_15motifs_nSeq300_nMotifs1_",
                     c("frac01/", "frac025/", "frac075/", "frac09/"))

#######################################################################################################
ranking_methods <- c(
  "QuiPT",
  "Chi-squared",
  
  "symuncert",
  "infogain",
  "gainratio",
  
  "NJMIM",
  "MRMR",
  "JMIM",
  "JMI",
  "DISR"
)

nonranking_methods <- c(
  #"symuncert_nonranking",
  #"infogain_nonranking",
  #"gainratio_nonranking",
  "QuiPT_nonranking",
  "Chi-squared_nonranking",
  "FCBF_nonranking"
)

nice_names <- list(
  symuncert="SU", QuiPT="QuiPT", NJMIM="NJMIM", MRMR="MRMR",
  JMIM="JMIM", JMI="JMI", infogain="IG", gainratio="GR", DISR="DISR",
  "Chi-squared"="Chi-squared", symuncert_nonranking="SU (nr)",
  QuiPT_nonranking ="QuiPT (nr)", infogain_nonranking = "IG (nr)",
  gainratio_nonranking = "GR (nr)", FCBF_nonranking = "FCBF",
  "Chi-squared_nonranking" = "Chi-squared (nr)"
)

# libraries
library(QuiPTsim)
library(drake)
library(pbapply)
library(ggplot2)
library(dplyr)
library(xtable)

theme_set(theme_bw(base_size = 16))

# objects
#cache <- drake_cache(cache_path)
#paths <- readd(paths, cache = cache)

fractions <- c(
  "10% positive sequences",
  "25% positive sequences",
  "75% positive sequences",
  "90% positive sequences"
)

ranking_results <- lapply(data_paths, function(data_path) {
                            lapply(ranking_methods, function(method){
                               output_files <- dir(data_path)[grep(paste0(method, "_", "[0-9]"), dir(data_path))]
                               result_files <- paste0(data_path, output_files)
                               aggregate_results(parse_results(result_files))
  })
})

nonranking_results <- lapply(data_paths, function(data_path) {
  lapply(nonranking_methods, function(method){
    output_files <- dir(data_path)[grep(paste0(method, "_", "[0-9]"), dir(data_path))]
    result_files <- paste0(data_path, output_files)
    aggregate_results(parse_results(result_files), ranking = FALSE)
  })
})

for (i in 1:length(ranking_results)) {
  names(ranking_results[[i]]) <- ranking_methods
  names(ranking_results[[i]]) <- unlist(lapply(names(ranking_results[[i]]), function(x) nice_names[[x]]))
}

for (i in 1:length(nonranking_results)) {
  names(nonranking_results[[i]]) <- nonranking_methods
  names(nonranking_results[[i]]) <- unlist(lapply(names(nonranking_results[[i]]), function(x) nice_names[[x]]))

  nonranking_results[[i]][["QuiPT 0.01"]] <- nonranking_results[[i]]$`QuiPT (nr)`[nonranking_results[[i]]$`QuiPT (nr)`$threshold == 0.01,]
  nonranking_results[[i]][["QuiPT 0.05"]] <- nonranking_results[[i]]$`QuiPT (nr)`[nonranking_results[[i]]$`QuiPT (nr)`$threshold == 0.05,]
  nonranking_results[[i]]$`QuiPT (nr)` <- NULL
  
  nonranking_results[[i]][["Chi-squared 0.01"]] <- nonranking_results[[i]]$`Chi-squared (nr)`[nonranking_results[[i]]$`Chi-squared (nr)`$threshold == 0.01,]
  nonranking_results[[i]][["Chi-squared 0.05"]] <- nonranking_results[[i]]$`Chi-squared (nr)`[nonranking_results[[i]]$`Chi-squared (nr)`$threshold == 0.05,]
  nonranking_results[[i]]$`Chi-squared (nr)` <- NULL
}


for (i in 1:length(ranking_results)){
  for (j in 1:length(ranking_methods)){
    ranking_results[[i]][[j]]["Filter"] <- ranking_methods[j] 
  }
}

for (i in 1:length(nonranking_results)){
  for (j in 1:length(names(nonranking_results[[i]]))){
    nonranking_results[[i]][[j]]["Filter"] <- names(nonranking_results[[i]])[j] 
  }
}

ranking_results <- lapply(ranking_results, function(x) do.call(rbind, x))
nonranking_results <- lapply(nonranking_results, function(x) do.call(rbind, x))

for (i in 1:length(ranking_results)){
  ranking_results[[i]]["Fraction"] <- fractions[i] 
}

for (i in 1:length(nonranking_results)){
  nonranking_results[[i]]["Fraction"] <- fractions[i] 
}

ranking_results <- do.call(rbind, ranking_results)
nonranking_results <- do.call(rbind, nonranking_results)

ranking_results_lm <- ranking_results[ranking_results$model == "lm", ]
nonranking_results_lm <-nonranking_results[nonranking_results$model == "lm", ]
ranking_results_lm$Filter <- unlist(lapply(ranking_results_lm$Filter, function(x) nice_names[[x]]))


metrics = "AUC"
ncol=2

new_levels <- c("Chi-squared", "QuiPT", "SU",
                "IG", "GR", "NJMIM","MRMR",
                "JMIM", "JMI", "DISR")

new_levels_nonranking <- c("FCBF", "QuiPT 0.01", "QuiPT 0.05", "Chi-squared 0.01", "Chi-squared 0.05")


ranking_results_lm[["Filter"]] <- factor(ranking_results_lm[["Filter"]], levels = new_levels)
nonranking_results_lm[["Filter"]] <- factor(nonranking_results_lm[["Filter"]], levels = new_levels_nonranking)
nonranking_results_lm[['n_kmers']] <- nonranking_results_lm[['n_kmers_mean']]
