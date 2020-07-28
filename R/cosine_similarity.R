#' calculates euclidean norm of a given vector
#' @param x vector
#' @return euclidean norm of `x`
#' @export
#' @examples
#' euclidean_norm(1:4)

euclidean_norm <- function(x) sqrt(sum(x^2))

#' calculates cosine similarity of two vectors
#' @param x first vector
#' @param y second vector
#' @return cosine similarity of two vectors
#' @export
#' @examples
#' cosine_similarity(1:4, 5:8)

cosine_similarity <- function(x, y) {

  stopifnot(length(x) == length(y))
  sum(x * y) / (sqrt(sum(x * x)) * sqrt(sum(y * y)))

}

#' function generates two vectors of probabilities with given cosine similarity
#' @param size expected size of probabilities' vectors
#' @param cosine_sim cosine similarity of two vectors
#' @return matrix with probabilites in rows
#' @export
#' @importFrom stats runif
#' @examples
#' generate_probs(5, 0.9)

generate_probs <- function(size, cosine_sim) {

  # generate vector of probs and norm it
  u <- runif(size, min = 0.25, max = 1)
  u <- u / sum(u)
  u_normed <- u / euclidean_norm(u)

  # pick random vector
  r <- runif(size)

  # calcualte vector perpendicular to u and norm it
  u_perp <- r - as.vector(r %*% u_normed) * u_normed
  u_perp_normed <- u_perp / euclidean_norm(u_perp)

  # calculate vector with given cosine similarity and make it vector of probabilites
  w <- cosine_sim * u_normed + sqrt(1 - cosine_sim ** 2) * u_perp_normed
  w <- w / sum(w)

  probs <- rbind(u, w)

  if (all(probs > 0)) {
    return(probs)
  } else {
    generate_probs(size, cosine_sim)
  }

}
