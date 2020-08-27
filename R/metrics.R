#' Sensitivity (true positive rate, recall)
#' @param y_true boolean vector of target variable
#' @param y_pred boolean vector of predictions
#' @export

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
#' @export

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
#' @export

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
#' @export

recall <- function (y_true, y_pred) {
  sensitivity(y_true, y_pred)
}

#' F1-score (harmonic mean of precision and recall)
#' @inheritParams sensitivity
#' @export

F1score <- function (y_true, y_pred) {

  if (!(length(y_true) == length(y_pred))) {
    stop("Target and predictions do not have the same length!")
  }

  prec <- precision(y_true, y_pred)
  rec <- recall(y_true, y_pred)

  (2 * prec * rec) / (prec + rec)
}

#' Aggregated scores
#' @inheritParams sensitivity
#' @export

compute_metrics <- function(y_true, y_pred) {

  list(sensitivity = sensitivity(y_true, y_pred),
       specificity = specificity(y_true, y_pred),
       F1score = F1score(y_true, y_pred),
       precision = precision(y_true, y_pred),
       recall = recall(y_true, y_pred))

}
