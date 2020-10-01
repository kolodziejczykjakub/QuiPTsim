library(QuiPTsim)

df <- read.csv("./reduced_alph_enc_alph4_const/alph4_const.csv")

paths <- df[df$l_seq == 10 & df$n_motifs == 2, "path"]
setup <- list(method = "FCBF")
validation_scheme <- list(type = "cv",
                          folds = 5,
                          learner = "classif.log_reg")

ranking_summary <- function(paths, setup, validation_scheme) {

  results <- pblapply(paths, function(path) {

    # read n-gram matrix
    # TODO: fraction, mixing matrices etc.
    m <- read_ngram_matrix(path)

    # scoring k-mers
    filtering_results <- filter_ngrams(m,setup[["method"]])


    # TODO: has to be done better
    if (setup[["method"]] == "FCBF") {

      X <- as.matrix(m[, filtering_results[["score"]]])
      y <- attr(m, "target")

      df <- data.frame(X, y)

      # Cross-validation

      # r = rsmp(validation_scheme[["type"]],
      #          folds=validation_scheme[["folds"]])
      filtering_results[filtering_results[["score"]], "ngram"]

      browser()
    }

  })

}

ranking_summary(paths, setup, validation_scheme)
