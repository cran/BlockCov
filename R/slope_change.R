#' This function fits to a numerical vector sorted in the non decreasing order two simple linear regressions and returns the index corresponding to the estimated change between the two regression models.
#'
#' @param  Y numerical vector sorted in the non decreasing order.
#' @return K the index corresponding to the estimated change between the two linear regression models.
#' @importFrom Matrix Matrix
#' @importFrom dplyr arrange filter mutate cummean
#' @importFrom tibble tibble rowid_to_column
#' @importFrom rlang .data
#' @examples
#' n <- 30
#' q <- 100
#' Sigma <- Simu_Sigma(q = q, diag = FALSE, equal = TRUE)
#' Matrix::image(Sigma)
#' E <- matrix(rnorm(n * q), ncol = q) %*% chol(as.matrix(Sigma))
#' corE <- cor(as.matrix(E))
#' vec_up_emp <- corE[upper.tri(corE)]
#' G <- matrix(0, ncol = (q - 1), nrow = (q - 1))
#' G[upper.tri(G, diag = TRUE)] <- vec_up_emp
#' G[lower.tri(G)] <- t(as.matrix(G))[lower.tri(t(as.matrix(G)))]
#' res_svd <- svd(G)
#' vp <- res_svd$d
#' slope_change(vp)
#' @export
slope_change <- function(Y) {
  tb <- tibble(y = Y) %>%
    arrange(desc(.data$y)) %>%
    rowid_to_column(var = "x") %>%
    mutate(
      mx = (.data$x + 1) / 2,
      my = cummean(.data$y),
      pxy = .data$x * .data$y,
      spxy = cumsum(.data$pxy),
      s2xy = .data$spxy - .data$mx * .data$x * .data$my,
      s2x = (.data$x * (.data$x + 1) * (.data$x - 1)) / 12,
      b = .data$s2xy / .data$s2x,
      a = .data$my - .data$b * .data$mx,
      y2 = cumsum(.data$y^2),
      e = .data$y2 + (.data$b^2 * (.data$x * (.data$x + 1) *
                                     (2 * .data$x + 1)) / 6) +
        2 * .data$a * .data$b * .data$x * .data$mx +
        .data$x * .data$a^2 -
        2 * .data$b * .data$spxy - 2 * .data$a * .data$my * .data$x
    ) %>%
    filter(.data$x != n())

  n <- length(Y)
  tb2 <- tibble(y = Y) %>%
    arrange(desc(.data$y)) %>%
    rowid_to_column(var = "x") %>%
    arrange(.data$y) %>%
    rowid_to_column(var = "k") %>%
    mutate(
      mx = (.data$x + n) / 2,
      my = cummean(.data$y),
      pxy = .data$x * .data$y,
      spxy = cumsum(.data$pxy),
      s2xy = .data$spxy - .data$mx * .data$k * .data$my,
      s2xi = (n * (n + 1) * (2 * n + 1) - .data$x * (.data$x - 1) *
                (2 * .data$x - 1)) / 6,
      s2x = .data$s2xi - .data$k * (n + .data$x)^2 / 4,
      b = .data$s2xy / .data$s2x,
      a = .data$my - .data$b * .data$mx,
      y2 = cumsum(.data$y^2),
      e = .data$y2 + (.data$b^2 * .data$s2xi) +
        2 * .data$a * .data$b * .data$k * .data$mx + .data$k * .data$a^2 -
        2 * .data$b * .data$spxy - 2 * .data$a * .data$my * .data$k
    ) %>%
    filter(.data$x != n())
  errors <- tb$e + rev(tb2$e)
  which.min(errors) - 1
}
