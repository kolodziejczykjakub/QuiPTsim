library(QuiPTsim)
library(pbapply)
library(caret)
library(pROC)

df <- read.csv("./reduced_alph_enc_alph4_const/alph4_const.csv")

paths <- df[df$l_seq == 10 & df$n_motifs == 2, "path"]
setup <- list(method = "FCBF")
validation_scheme <- list(type = "cv",
                          folds = 5)

m <- read_ngram_matrix(paths[1])
filtering_results <- filter_ngrams(m,setup[["method"]])

X <- as.matrix(m[, filtering_results[["score"]]])
y <- attr(m, "target")
df <- data.frame(X, y)

df
cvReplications <- 1
n_folds = 5
for (i in 1:cvReplications) {

}
folds <- createFolds(y = df$y, k = n_folds)

for (foldNum in 1:n_folds) {

  df_train <- df[-folds[[foldNum]], ]
  df_test <- df[folds[[foldNum]], ]
}

# linear model
logReg <- glm(y ~ ., family = binomial(link = "logit"), data = df_train)

y_true <- df_test[["y"]]
y_proba <- predict(logReg, df_test)
y_pred <- y_proba > 0.5
