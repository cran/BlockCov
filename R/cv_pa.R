#' Title
#'
#' @param E the observation matrix such that each of its row has a block structure correlation matrix Sigma to estimate up to a permutation of its columns and rows.
#' @param v_ord the absolute value of the upper  triangular part matrix \eqn{\Gamma} (including its diagonal) order in
#' increasing order
#' @param N number of replication in the "cross-validation"
#'
#' @return  the number of non null values selected for the estimation of the covariance matrix
#' @details In order to get the treshold one must do rev(v_ord)[cv_bl(E, v_ord, N=N)]
#' @export
#'
#' @examples
#' n <- 30
#' q <- 100
#' Sigma <- Simu_Sigma(q = q, diag = FALSE, equal = TRUE)
#' Matrix::image(Sigma)
#' E <- matrix(rnorm(n * q), ncol = q) %*% chol(as.matrix(Sigma))
#' k <- 5
#' v_up <- est_up(E, k = k)
#' a_vup <- abs(v_up)
#' ord_vup <- order(a_vup)
#' v_ord <- a_vup[ord_vup]
#' N <- 10
#' nb_nn0 <- cv_bl(E, v_ord, N=N)
#' tresh <- rev(v_ord)[nb_nn0]

cv_bl<- function(E, v_ord, N){
  n <- nrow(E)
  n1 <- round(n*(1- 1/log(n)))

  r_hat <- lapply(1:N, function(i){
    s1 <- sample(seq_len(n), n1)
    s2 <- seq_len(n)[-s1]
    v1 <- est_up(E[s1, ])
    v2 <- est_up(E[s2, ])
    ord <- order(abs(v1))
    v1 <- v1[ord]
    v2 <- v2[ord]
    dif <- (v1 - v2)^2
    ajout <- v2^2 - dif
    fin <- c(sum(dif), ajout)
    r_hat <- cumsum(fin)
    reord <- findInterval(v_ord, abs(v1), left.open = TRUE)
    r_hat[c((reord + 1))]
  })
  r_hat_mean <- Reduce("+", r_hat)
  v_ord[which.min(r_hat_mean)]
  which.min(rev(r_hat_mean))
}



#' Title
#'
#' @param E the observation matrix such that each of its row has a block structure correlation matrix Sigma wich has a low rank once its diagonal is removed.
#' @param times number of random sampling
#'
#' @return the mean of the eigen values of the \code{times} sampled matrix
#' @export
#'
#' @examples
#' n <- 30
#' q <- 100
#' Sigma <- Simu_Sigma(q = q, diag = FALSE, equal = TRUE)
#' Matrix::image(Sigma)
#' E <- matrix(rnorm(n * q), ncol = q) %*% chol(as.matrix(Sigma))
#' random_eigen <- PA(E, times = 10)
PA <- function(E, times = 10){
  q<-ncol(E)
  Reduce("+",  lapply(1:times, function(lalala){
    corEs <- cor(as.matrix(apply(E,2,sample)))
    Pti_sigs <- Matrix(0, ncol = (q - 1), nrow = (q - 1))
    Pti_sigs[upper.tri(Pti_sigs, diag = TRUE)] <- corEs[upper.tri(corEs)]
    Pti_sigs[lower.tri(Pti_sigs)] <- t(as.matrix(Pti_sigs))[lower.tri(t(as.matrix(Pti_sigs)))]
    res_svd_corE <- svd(as.matrix(Pti_sigs))
    res_svd_corE$d
  })) / times
}

#' Title
#'
#' @param E the observation matrix such that each of its row has a block structure correlation matrix Sigma wich has a low rank once its diagonal is removed.
#' @param k the rank of the correlation matrix of \code{E} once its diagonal has been removed
#'
#' @return an approximation of the correlation matrix of \code{E} with its diagonal removed
#' @export
#'
#' @examples
#' n <- 30
#' q <- 100
#' Sigma <- Simu_Sigma(q = q, diag = FALSE, equal = TRUE)
#' Matrix::image(Sigma)
#' E <- matrix(rnorm(n * q), ncol = q) %*% chol(as.matrix(Sigma))
#' k <- 5
#' v_up <- est_up(E, k = k)
est_up <- function(E, k = 5){
  q <- ncol(E)
  corE <- cor(as.matrix(E))
  Pti_sig <- Matrix(0, ncol = (q - 1), nrow = (q - 1))
  Pti_sig[upper.tri(Pti_sig, diag = TRUE)] <- corE[upper.tri(corE)]
  Pti_sig[lower.tri(Pti_sig)] <- t(as.matrix(Pti_sig))[lower.tri(t(as.matrix(Pti_sig)))]
  res_svd_corE <- svd(as.matrix(Pti_sig), nu =k, nv=k)
  vp_corE <- res_svd_corE$d
  U_corE <- res_svd_corE$u
  tV_corE <- t(res_svd_corE$v)
  largest_vp_corE <- vp_corE[1:k]
  corE_aprx <- U_corE[, 1:k] %*% diag(largest_vp_corE) %*% tV_corE[1:k, ]
  return(corE_aprx [upper.tri(corE_aprx)])
}
