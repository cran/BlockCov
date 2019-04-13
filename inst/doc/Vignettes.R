## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message=FALSE)
library(BlockCov)
set.seed(516)

## ---- , eval =FALSE------------------------------------------------------
#  devtools::install_github("Marie-PerrotDockes/BlockCov")

## ------------------------------------------------------------------------
q <- 100
Sigma <- Simu_Sigma(q = q, diag = FALSE, equal = TRUE)

## ----fig0, fig.cap="\\label{fig:fig0}",fig.width=3.5,fig.height=3.5,echo=FALSE----
Matrix::image(Sigma)

## ------------------------------------------------------------------------
n <- 30
E <- matrix(rnorm(n * q), ncol = q) %*% chol(as.matrix(Sigma))

## ------------------------------------------------------------------------
k <- 5
nb_nn0 <- sum(Sigma[upper.tri(Sigma, diag = FALSE)] != 0)
res_known <-  Sigma_estimation(E, k = k, nb_nn0 = nb_nn0)

## ----fig1, fig.cap="\\label{fig:fig1}",fig.width=3.5,fig.height=3.5------
Matrix::image(res_known$Sigma_est)

## ----fig2, fig.cap="\\label{fig:fig2}",fig.width=3.5,fig.height=3.5------
Matrix::image(Matrix::Matrix(cor(E)))

## ----warning=FALSE-------------------------------------------------------
res <-Sigma_estimation(E, method_k = "Cattell", method_0 = "Elbow")

## ---- eval = FALSE-------------------------------------------------------
#  res <-Sigma_estimation(E)

## ------------------------------------------------------------------------
res_pabl <- Sigma_estimation(E, method_k = "PA", method_0 = "BL")

## ----fig3, fig.cap="\\label{fig:fig3}",fig.width=3.5,fig.height=3.5------
Matrix::image(res$Sigma_est)

## ----fig3pabl, fig.cap="\\label{fig:fig3pabl}",fig.width=3.5,fig.height=3.5----
Matrix::image(res_pabl$Sigma_est)

## ------------------------------------------------------------------------
res_both <- Sigma_estimation(E, method_k = "Cattell", method_0 = "Elbow", inv_12 = TRUE)

## ------------------------------------------------------------------------
samp <- sample(1:q, q, replace = FALSE)
Sigma_samp <- Sigma[samp, samp]

## ----fig4, fig.cap="\\label{fig:fig4}",fig.width=3.5,fig.height=3.5------
Matrix::image(Sigma_samp)

## ------------------------------------------------------------------------
E <- matrix(rnorm(n * q), ncol = q) %*% chol(as.matrix(Sigma_samp))
res_samp <- Sigma_estimation(E, reorder = TRUE, inv_12 = TRUE)

## ----fig5, fig.cap="\\label{fig:fig5}",fig.width=3.5,fig.height=3.5------
Matrix::image(res_samp$Sigma_est)

## ----fig6, fig.cap="\\label{fig:fig6}",fig.width=3.5,fig.height=3.5------
ord <- res_samp$order
Matrix::image(res_samp$Sigma_est[ord, ord])

## ----fig7, fig.cap="\\label{fig:fig7}",fig.width=3.5,fig.height=3.5------
Matrix::image(Sigma_samp[ord, ord])

## ----fig8, fig.cap="\\label{fig:fig8}",fig.width=3.5,fig.height=3.5------
Matrix::image(res_samp$S_inv_12 %*% Sigma_samp %*%res_samp$S_inv_12)

