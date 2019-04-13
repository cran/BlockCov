#' This function computes an estimator of the covariance matrix and the square root of its inverse and permutes its rows and columns if it is necessary to make the block structure appear.
#'
#' @param  E the observation matrix such that each of its row has a block structure correlation matrix Sigma to estimate up to a permutation of its columns and rows.
#' @param  k numerical or NULL, the rank for the low rank approximation. If NULL the rank is computed using the slope_change function applied on the eigenvalues of the low rank part of Sigma. Default to NULL.
#' @param  nb_nn0 numerical or NULL, corresponds to the number of non null values to keep in the estimation of the covariance matrix.
#' If NULL the number of non null values is computed using the slope_change function to the Frobenius norm of the difference between the empirical correlation matrix and its estimation with nb_nn0 non null values. Default to NULL.
#' @param  method_k character if "Cattell" (the default) then the Cattell criterion \insertCite{cattell1966}{BlockCov} is performed on the singular values of the covariance matrix.
#' to estimate the number of rank use in the low rank approximation, while "PA" use the parrallel analysis \insertCite{Horn1965}{BlockCov}
#' wich can be more accurate if the number of rows of E is not to small but which is much slower.
#' @param times numeric the number of resampling done for the "PA" method, ignored if metod_k is different from "PA".
#' @param  method_0 character if "Elbow" (the default) then the Elbow criterion (see \insertCite{blc;textual}{BlockCov} for details) is performed
#'   to estimate the number of rank use in the low rank approximation, while "BL" use the approach proposed in
#'   \insertCite{bickel2008;textual}{BlockCov} based on cross-validation
#'    wich can be more accurate if the number of rows of E is not to small but which is much slower.
#' @param N numeric the number of fold used for the "BL" method. Ignored if method_0 is different from "BL"
#' @param  big logical, default to FALSE. If the dataset is too big the empirical correlation is calculated by crossprod(E) * 1 / n to fasten the computation
#' @param  reorder logical, default to FALSE. Whether or not the columns of E are permuted. If TRUE a hierarchical clustering is first performed and the columns are permuted according to it.
#' @param  inv_12 logical, default to FALSE Whether or not computing the square root of the inverse of the covariance matrix.
#' @return A list with the elements
#' \item{Sigma_est}{estimator of the covariance matrix}
#' \item{k}{rank of the low rank part of the covariance matrix}
#' \item{nb_nn0}{number of non null values of the upper triangular part of the covariance matrix}
#' \item{S_inv_12}{square root of the inverse of the estimated covariance matrix}
#' \item{order}{permutation to apply to the rows and the columns of the covariance to make the block structure appear}
#' @importFrom Matrix Matrix nearPD t image
#' @importFrom stats cor dist hclust runif
#' @importFrom BBmisc which.last
#'@importFrom Rdpack reprompt
#'@importFrom stats as.dist
#'@importFrom dplyr desc
#'@importFrom dplyr n
#' @examples
#' n <- 30
#' q <- 100
#' Sigma <- Simu_Sigma(q = q, diag = FALSE, equal = TRUE)
#' Matrix::image(Sigma)
#' E <- matrix(rnorm(n * q), ncol = q) %*% chol(as.matrix(Sigma))
#' res <- Sigma_estimation(E, inv_12 = TRUE)
#' Matrix::image(res$Sigma_est)
#' Matrix::image(res$S_inv_12)
#'@references
#'\insertAllCited{}
#' @export
Sigma_estimation <- function(E, k = NULL,
                             nb_nn0 = NULL,
                             big = FALSE, reorder = FALSE,
                             inv_12 = FALSE,
                             method_k = "Cattell",
                             times =10,
                             method_0="Elbow",
                             N=10) {
  ord <- NULL
  q <- ncol(E)
  n <- nrow(E)
  if (!big) {
    corE <- cor(as.matrix(E))
  } else {
    corE <- crossprod(E) * 1 / n
  }

  if (reorder) {
    clust <- hclust(as.dist(1 - corE))
    ord <- clust$order
    E <- E[, ord]
    if (!big) {
      corE <- cor(as.matrix(E))
    } else {
      corE <- crossprod(E) * 1 / n
    }
  }

  vec_up_emp <- corE[upper.tri(corE)]

  Pti_sig <- Matrix(0, ncol = (q - 1), nrow = (q - 1))
  Pti_sig[upper.tri(Pti_sig, diag = TRUE)] <- vec_up_emp
  Pti_sig[lower.tri(Pti_sig)] <- t(as.matrix(Pti_sig))[lower.tri(t(as.matrix(Pti_sig)))]


  res_svd <- svd(Pti_sig)
  vp <- res_svd$d
  u <- res_svd$u
  tv <- t(res_svd$v)

  if (is.null(k)){
    if(method_k =="Cattell")  k <- max(slope_change(vp), 2)
    if(method_k =="PA")  k <- which.last(vp > max(PA(E, times)))
  }
  large_vp <- vp[1:k]
  corE_aprx <- u[, 1:k] %*% diag(large_vp) %*% tv[1:k, ]
  vec_up <- corE_aprx[upper.tri(corE_aprx, diag = TRUE)]
  l <- length(vec_up)
  a_vup <- abs(vec_up)
  ord_vup <- order(a_vup)
  v_ord <- a_vup[ord_vup]

  if (is.null(nb_nn0)) {

    if(method_0 =="Elbow"){
    error <- c(0,cumsum(v_ord^2))
    nb_nn0 <- slope_change(rev(error))
    }
    if(method_0 == "BL"){
  nb_nn0 <- cv_bl(E, v_ord, N=N)
    }
  }
  seuil <- rev(v_ord)[nb_nn0]
  vec_up[a_vup < seuil] <- 0
  Sig_est <- Matrix(0, q, q)
  Sig_est[upper.tri(Sig_est)] <- vec_up
  Sig_est <- Sig_est + t(Sig_est)
  diag(Sig_est) <- 1


  if (reorder) {
    reord <- order(ord)
    Sig_est <- Sig_est[reord, reord]
  }
  Sig_est <- nearPD(Sig_est, corr = TRUE)$mat


  if (inv_12) {
    res_svd <- svd(Sig_est)
    vp2 <- res_svd$d
    u <- res_svd$u
    tv <- t(res_svd$v)

    sel <- which(vp2 > 0.1)
    vp[-sel] <- 0
    S_inv_12 <- u[, sel] %*% diag(sqrt(1 / vp2[sel])) %*% tv[sel, ]
  } else {
    S_inv_12 <- NULL
  }
  return(list(Sigma_est = Sig_est, k = k, nb_nn0 = nb_nn0,
              S_inv_12 = S_inv_12, order = ord))
}
