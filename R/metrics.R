#' Sensitivity (true positive rate, recall)
#' @param y_true boolean vector of target variable
#' @param y_pred boolean vector of predictions

sensitivity <- function (y_true, y_pred) {

  if (!(length(y_true) == length(y_pred))) {
    stop("Target and predictions do not have the same length!")
  }

  true_positives <- sum(y_true & y_pred)
  false_negatives <- sum(y_true & !y_pred)

  true_positives / (true_positives + false_negatives)
}

#' Specificity (true negative rate)
#' @inheritParams sensitivity

specificity <- function (y_true, y_pred) {

  if (!(length(y_true) == length(y_pred))) {
    stop("Target and predictions do not have the same length!")
  }

  true_negatives <- sum(!y_true & !y_pred)
  false_positives <- sum(!y_true & y_pred)

  true_negatives / (true_negatives + false_positives)
}

#' Precision
#' @inheritParams sensitivity

precision <- function (y_true, y_pred) {

  if (!(length(y_true) == length(y_pred))) {
    stop("Target and predictions do not have the same length!")
  }

  true_positives <- sum(y_true & y_pred)
  false_positives <- sum(!y_true & y_pred)

  true_positives / (true_positives + false_positives)
}

#' Recall
#' @inheritParams sensitivity

recall <- function (y_true, y_pred) {
  sensitivity(y_true, y_pred)
}

#' F1-score (harmonic mean of precision and recall)
#' @inheritParams sensitivity
F1score <- function (y_true, y_pred) {

  if (!(length(y_true) == length(y_pred))) {
    stop("Target and predictions do not have the same length!")
  }

  prec <- precision(y_true, y_pred)
  rec <- recall(y_true, y_pred)

  (2 * prec * rec) / (prec + rec)
}

