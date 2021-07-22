experiment_name <- "Experiment on Amylogram data"
exp_prefix <- "exp_amylogram"

data_path <- "AmylogramExp_results/"
cache_path <- "AmylogramExp/"


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

dir(data_path)

ranking_paths <- paste0(data_path, dir(data_path)[-grep("nonranking", dir(data_path))])
nonranking_paths <- paste0(data_path, dir(data_path)[grep("nonranking", dir(data_path))])
nonranking_paths <- nonranking_paths[c(1, 2, 5)]

cache <- drake_cache(cache_path)

readd(paths, cache = cache)


ranking_results <- lapply(ranking_paths, parse_results)
nonranking_results <- lapply(nonranking_paths, parse_results)

ranking_methods <- unlist(lapply(ranking_paths, function(x) substr(x, 29, nchar(x) - 4)))
nonranking_methods <- unlist(lapply(nonranking_paths, function(x) substr(x, 29, nchar(x) - 4)))

names(nonranking_results) <- nonranking_methods
names(nonranking_results) <- unlist(lapply(names(nonranking_results), function(x) nice_names[[x]]))

names(ranking_results) <- ranking_methods
names(ranking_results) <- unlist(lapply(names(ranking_results), function(x) nice_names[[x]]))

for (j in 1:length(ranking_methods)){
  ranking_results[[j]]["Filter"] <- names(ranking_results)[j] 
}

for (j in 1:length(nonranking_methods)){
  nonranking_results[[j]]["Filter"] <- names(nonranking_results)[j] 
}



nonranking_results[["QuiPT 0.01"]]  <- nonranking_results$`QuiPT (nr)`[nonranking_results$`QuiPT (nr)`$threshold == 0.01,]
nonranking_results[["QuiPT 0.05"]]  <- nonranking_results$`QuiPT (nr)`[nonranking_results$`QuiPT (nr)`$threshold == 0.05,]
nonranking_results$`QuiPT (nr)` <- NULL
nonranking_results[["QuiPT 0.01"]][["Filter"]] <- "QuiPT 0.01"
nonranking_results[["QuiPT 0.05"]][["Filter"]] <- "QuiPT 0.05"

nonranking_results[["Chi-squared 0.01"]] <- nonranking_results$`Chi-squared (nr)`[nonranking_results$`Chi-squared (nr)`$threshold == 0.01,]
nonranking_results[["Chi-squared 0.05"]] <- nonranking_results$`Chi-squared (nr)`[nonranking_results$`Chi-squared (nr)`$threshold == 0.05,]
nonranking_results$`Chi-squared (nr)` <- NULL
nonranking_results[["Chi-squared 0.01"]][["Filter"]] <- "Chi-squared 0.01"
nonranking_results[["Chi-squared 0.05"]][["Filter"]] <- "Chi-squared 0.05"

aggregated_nonranking_results <- do.call(rbind, nonranking_results)
aggregated_ranking_results <- do.call(rbind, ranking_results)

## Plots
plot_ranking_results_w_nonranking <- function(ranking_results, nonranking_results, metrics, ncol=3){
  

  total_results <- ranking_results
  total_results_nonranking <- nonranking_results
    
  new_levels <- c("Chi-squared", "QuiPT", "SU",
                  "IG", "GR", "NJMIM","MRMR",
                  "JMIM", "JMI", "DISR")

  new_levels_nonranking <- c("FCBF", "QuiPT 0.01", "QuiPT 0.05", "Chi-squared 0.01", "Chi-squared 0.05")


  total_results[["Filter"]] <- factor(total_results[["Filter"]],
                                      levels = new_levels)

  total_results_nonranking[["Filter"]] <- factor(total_results_nonranking[["Filter"]],
                                                 levels = new_levels_nonranking)
  
  
  total_results %>%
    ggplot(aes_string(x="n_kmers", y=metrics)) +
    geom_line(aes_string(linetype="Filter", group="Filter", color="Filter"), size = 1) +
    geom_point(data = total_results_nonranking, mapping = aes_string(x="n_kmers", y=metrics, shape="Filter"), size=4) +
    facet_wrap(~Model, ncol=ncol) +
    scale_x_continuous(name="Number of selected k-mers",  trans='log2') +
    scale_linetype_manual(values=rep(c("solid", "dotdash"), each=5)) +
    scale_color_manual(values=rep(c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00'), 3)) +
    scale_shape_manual(values=c(8, 15:18)) + 
    ylab(toupper(metrics))
}

## Results


