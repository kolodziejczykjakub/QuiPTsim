# Setup
experiment_name <- "Baseline (1 motif, 300 sequences)"
nSeq <- 300
motif <- 1
data_path <- "../data/experiment-results/exp01-seqLen10-nSeq300-amylogram-encoding/"
cache_path <- "../data/experiment-drake-caches/exp01_seqLen10_nSeq300_amylogram_encoding/"

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

theme_set(theme_bw())
# functions
plot_times <- function(total_times) {

  melted_times <- reshape2::melt(total_times)

  g <- ggplot(melted_times, aes(L1, value))
  g + geom_boxplot() +
    theme(axis.text.x = element_text(angle=65, vjust=0.6)) +
    labs(title="Average computation time",
         subtitle= experiment_name,
         x="Filtering method",
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

plot_ranking_results <- function(ranking_results, metrics){

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
    facet_wrap(~Model) +
    scale_linetype_manual(values=rep(c("solid", "dotdash", "dotted"), each=4)) +
    scale_color_manual(values=rep(c('#d7191c','#fdae61','#abdda4','#2b83ba'), 3)) +
    # ['#a6cee3','#1f78b4','#b2df8a','#33a02c']
    # '#d7191c','#fdae61','#abdda4','#2b83ba'
    theme_bw()
}

# objects
cache <- drake_cache(cache_path)
paths <- readd(paths, cache = cache)
num_reps <- grep(paste0("nMotifs_", motif), paths)

# aggregated times
ranking_times <- pblapply(ranking_methods, function(method) {
  result_files <- paste0(data_path, "result__", method, "_", num_reps, ".Rds")
  collect_filtering_times(result_files)
})

nonranking_times <- pblapply(nonranking_methods, function(method) {
  result_files <- paste0(data_path, "result__", method, "_", num_reps, ".Rds")
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
  aggregate_results(parse_results(paste0(data_path, "result__", method, "_", num_reps, ".Rds")), ranking = FALSE)
})

ranking_results <- lapply(ranking_methods, function(method) {
  aggregate_results(parse_results(paste0(data_path, "result__", method, "_", num_reps, ".Rds")))
})

names(nonranking_results) <- nonranking_methods
names(ranking_results) <- ranking_methods

names(nonranking_results) <- unlist(lapply(names(nonranking_results), function(x) nice_names[[x]]))
names(ranking_results) <- unlist(lapply(names(ranking_results), function(x) nice_names[[x]]))
