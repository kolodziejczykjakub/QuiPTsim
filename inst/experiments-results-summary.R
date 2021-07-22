# Setup ###############################################################################################
experiment_name <- "Experiment 3 (2 motifs injected, 15 random motifs, 300 sequences, 50\\% positive sequences)"
exp_prefix <- "exp3_2m_300s_set15"
#nSeq <- 600
#motif <- 2
data_path <- "../data/experiment-results/exp-03/exp3_reduced_alph_enc_15motifs_amylogram_encoding/exp03-nSeq300-nMotifs2-frac05/"
cache_path <- "../data/experiment-drake-caches/exp03-15motifs/exp03_15motifs_nSeq300_nMotifs2_frac05/"
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
# functions
plot_times <- function(total_times) {

  melted_times <- reshape2::melt(total_times)

  g <- ggplot(melted_times, aes(L1, value))
  g + geom_boxplot() +
    theme(axis.text.x = element_text(angle=65, vjust=0.6)) +
    labs(x="Filtering method",
         y="Time (seconds in log scale)") +
    scale_y_continuous(trans='log2')
}

table_times <- function(total_times) {
  avg_time <- data.frame(`Average time` = sapply(total_times, mean), `Time sd` = sapply(total_times, sd))
  avg_time <- round(avg_time, 2)
  nice_times <- data.frame(Time = paste0(avg_time$Average.time, " +- ", avg_time$Time.sd))
  rownames(nice_times) <- rownames(avg_time)

  nice_times
}

table_nonranking <- function(nonranking_results, model, metrics_vec) {

  ans <- lapply(c("n_kmers", metrics_vec), function(metrics) {

    out <- do.call(rbind,
                   lapply(nonranking_results, function(x) {
                     data.frame(x[x$Model == model,
                                  c(paste0(metrics, "_mean"),
                                    paste0(metrics, "_std"))])
                   })
    )
    out <- round(out, 3)

    df <- data.frame(x = paste0(out[[paste0(metrics, "_mean")]],
                                " +- ",
                                out[[paste0(metrics, "_std")]]))

    colnames(df) <- ifelse(metrics == "n_kmers", "n-kmers selected", metrics)
    rownames(df) <- rownames(out)
    df
    
  })
  do.call(cbind, ans)
}

# TODO:
# table_ranking <- function(nonranking_results, model, metrics_vec) {
#
# }

plot_ranking_results <- function(ranking_results, metrics, ncol=3){

  for (filtername in names(ranking_results)){
    ranking_results[[filtername]][["Filter"]] = filtername
  }

  total_results <- do.call(rbind, ranking_results)

  new_levels <- c("Chi-squared", "QuiPT", "SU",
                  "IG", "GR", "NJMIM","MRMR",
                  "JMIM", "JMI", "DISR")

  total_results[["Filter"]] <- factor(total_results[["Filter"]],
                                      levels = new_levels)

  total_results %>%
    ggplot(aes_string(x="n_kmers", y=paste0(metrics, "_mean"), group="Filter", color="Filter")) +
    geom_line(aes(linetype=Filter), size = 1) + scale_x_continuous(name="Number of selected k-mers",  trans='log2') +
    facet_wrap(~Model, ncol=ncol) +
    scale_linetype_manual(values=rep(c("solid", "dotdash"), each=5)) +
    scale_color_manual(values=rep(c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00'), 2)) +
    ylab(metrics)
}

plot_ranking_results_w_nonranking <- function(ranking_results, nonranking_results, metrics, ncol=3){
  
  for (filtername in names(ranking_results)){
    ranking_results[[filtername]][["Filter"]] = filtername
  }
  
  for (filtername in names(nonranking_results)){
    nonranking_results[[filtername]][["Filter"]] = filtername
  }
  
  
  total_results <- do.call(rbind, ranking_results)
  total_results_nonranking <- do.call(rbind, nonranking_results)
  
  new_levels <- c("Chi-squared", "QuiPT", "SU",
                  "IG", "GR", "NJMIM","MRMR",
                  "JMIM", "JMI", "DISR")
  
  new_levels_nonranking <- c("FCBF", "QuiPT 0.01", "QuiPT 0.05", "Chi-squared 0.01", "Chi-squared 0.05")
  
  
  total_results[["Filter"]] <- factor(total_results[["Filter"]],
                                      levels = new_levels)
  
  total_results_nonranking[["Filter"]] <- factor(total_results_nonranking[["Filter"]],
                                      levels = new_levels_nonranking)
  total_results_nonranking[['n_kmers']] <- total_results_nonranking[['n_kmers_mean']]
  total_results %>%
    ggplot(aes_string(x="n_kmers", y=paste0(metrics, "_mean"))) +
    geom_line(aes_string(linetype="Filter", group="Filter", color="Filter"), size = 1) +
    geom_point(data = total_results_nonranking, mapping = aes_string(x="n_kmers", y=paste0(metrics, "_mean"), shape="Filter"), size=4) +
    facet_wrap(~Model, ncol=ncol) +
    scale_x_continuous(name="Number of selected k-mers",  trans='log2') +
    scale_linetype_manual(values=rep(c("solid", "dotdash"), each=5)) +
    scale_color_manual(values=rep(c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00'), 3)) +
    scale_shape_manual(values=c(8, 15:18)) + 
    ylab(metrics)
}



# objects
cache <- drake_cache(cache_path)
paths <- readd(paths, cache = cache)

# Exp 1
#num_reps <- grep(paste0("nMotifs_", motif), paths)

# aggregated times
ranking_times <- pblapply(ranking_methods, function(method) {
  #output_files <- dir(data_path)[grep(paste0(method, "_", num_reps, ".Rds", collapse="|"), dir(data_path))]
  output_files <- dir(data_path)[grep(paste0(method, "_", "[0-9]"), dir(data_path))]
  result_files <- paste0(data_path, output_files)
  collect_filtering_times(result_files)
})

nonranking_times <- pblapply(nonranking_methods, function(method) {
  #output_files <- dir(data_path)[grep(paste0(method, "_", num_reps, ".Rds", collapse="|"), dir(data_path))]
  output_files <- dir(data_path)[grep(paste0(method, "_", "[0-9]"), dir(data_path))]
  result_files <- paste0(data_path, output_files)
  collect_filtering_times(result_files)
})

names(ranking_times) <- ranking_methods
names(nonranking_times) <- nonranking_methods

total_times <- c(ranking_times, nonranking_times)
total_times$`Chi-squared_nonranking` <- NULL
total_times$QuiPT_nonranking <- NULL

names(total_times) <- unlist(lapply(names(total_times), function(x) nice_names[[x]]))

# aggregated results

nonranking_results <- lapply(nonranking_methods, function(method) {
  #output_files <- dir(data_path)[grep(paste0(method, "_", num_reps, ".Rds", collapse="|"), dir(data_path))]
  output_files <- dir(data_path)[grep(paste0(method, "_", "[0-9]"), dir(data_path))]
  result_files <- paste0(data_path, output_files)
  aggregate_results(parse_results(result_files), ranking = FALSE)
})

ranking_results <- lapply(ranking_methods, function(method) {
  #output_files <- dir(data_path)[grep(paste0(method, "_", num_reps, ".Rds", collapse="|"), dir(data_path))]
  output_files <- dir(data_path)[grep(paste0(method, "_", "[0-9]"), dir(data_path))]
  result_files <- paste0(data_path, output_files)
  aggregate_results(parse_results(result_files))
})

names(nonranking_results) <- nonranking_methods
names(ranking_results) <- ranking_methods

names(nonranking_results) <- unlist(lapply(names(nonranking_results), function(x) nice_names[[x]]))
names(ranking_results) <- unlist(lapply(names(ranking_results), function(x) nice_names[[x]]))

nonranking_results[["QuiPT 0.01"]] <- nonranking_results$`QuiPT (nr)`[nonranking_results$`QuiPT (nr)`$threshold == 0.01,]
nonranking_results[["QuiPT 0.05"]] <- nonranking_results$`QuiPT (nr)`[nonranking_results$`QuiPT (nr)`$threshold == 0.05,]
nonranking_results$`QuiPT (nr)` <- NULL

nonranking_results[["Chi-squared 0.01"]] <- nonranking_results$`Chi-squared (nr)`[nonranking_results$`Chi-squared (nr)`$threshold == 0.01,]
nonranking_results[["Chi-squared 0.05"]] <- nonranking_results$`Chi-squared (nr)`[nonranking_results$`Chi-squared (nr)`$threshold == 0.05,]
nonranking_results$`Chi-squared (nr)` <- NULL


