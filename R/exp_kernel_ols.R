#' Linear Model Adequacy Test
#' @description Some description
#'
#' @param x parameter
#' @examples
#' Include example
#' @export
#' 
#' 

lin_test <- function(x, y, k = 20, gamma = (sqrt(5)-1)/2){
  # poly_order <- d_tuple_k_sum(x_dim = ncol(x), k = k)
  # out <- lin_test_cpp(x, y, k, gamma, t(poly_order))
  if(!is.matrix(x)){
    x <- as.matrix(x)
  }
  out_cpp <- lin_test_cpp(x, y, k, gamma)
  out <- list(stat = out_cpp$stat, p_value = out_cpp$p_value, 
              prob = out_cpp$prob, dist = list(x = out_cpp$dist[, 1], 
                                               pdf = out_cpp$dist[, 2], 
                                               cdf = out_cpp$dist[, 3]))
  return(out)
}

#' Linear Model Adequacy Test
#' @description Some description
#'
#' @param x parameter
#' @examples
#' Include example
#' @export
#' 
#' 
lin_test_mc <- function(x, y, n_boot = 100, k_grid = 500){
  
  # make sure the univariate case works
  if(!is.matrix(x)){
    x <- as.matrix(x)
  }
  out <- lin_test_mc_cpp(x, y, n_boot, k_grid, 'exp_kernel', 0.0)
  
  return(list(p_value = mean(out$stat_boot >= out$stat), prob = mean(out$stat_boot <= out$stat), 
              stat = out$stat, stat_boot = out$stat_boot, rhs = out$rhs_u))
}

#' Linear Model Adequacy Test
#' @description Some description
#'
#' @param x parameter
#' @examples
#' Include example
#' @export
#' 
#' 
lin_test_mc_bierens <- function(x, y, n_boot = 100, k_grid = 500, tau = pi/2){
  
  # make sure the univariate case works
  if(!is.matrix(x)){
    x <- as.matrix(x)
  }
  out <- lin_test_mc_cpp(x, y, n_boot, k_grid, 'Bierens', tau)
  
  return(list(p_value = mean(out$stat_boot >= out$stat), prob = mean(out$stat_boot <= out$stat), 
              stat = out$stat, stat_boot = out$stat_boot, rhs = out$rhs_u))
}


# lin_test_uni <- function(x, y, k = 20, gamma = (sqrt(5)-1)/2, export_cdf = FALSE){
#   
#   psi <- Psi_mat(x/sd(x))
# 
#   x <- scale(x, scale = FALSE)
#   y <- scale(y, scale = FALSE)
#   beta <- solve(t(x) %*% x)%*%(t(x)%*%y)
#   resid <- y - x%*%beta
#   stat <- t(resid)%*%psi%*%resid/length(x)
#   
#   sig_res <- sd(resid)
#   
#   delta_k <- delta_uni_k_ols(x, k, gamma)
#   C <- (sig_res^2)*cov(delta_k)
#   eigvals <- gamma^(0:k)
#   
#   C <- diag(sqrt(eigvals))%*%C%*%diag(sqrt(eigvals))
#   lambdas <- abs(eigen(C)$values)
#   
#   out <- pQ(stat, lambdas)
#   
#   if(export_cdf){
#     cdf <- pQ_fft_cpp(lambdas, 0.005, N = 2^13)
#   }else{
#     cdf = NULL
#   }
#   
#   return(list(stat = stat, p_value = 1 - out, prob = out, cdf = cdf))
# }
# 
# 
# delta_uni_k_ols <- function(x, k, gamma){
#   Q_m <- (t(x)%*%x)/(length(x))
#   
#   x_mat <- rep(x, k+1)
#   dim(x_mat)<- c(length(x), k+1)
#   psi_k <- psi_uni_k(x/sd(x), k, gamma)
#   mean_xpsi_k <- colMeans(x_mat*psi_k)
#   mean_xpsi_k <- rep(mean_xpsi_k, length(x))
#   dim(mean_xpsi_k)<- c(k+1, length(x))
#   mean_xpsi_k <- t(mean_xpsi_k)
#   
#   u_vec <- t(x)
#   u_vec <- rep(u_vec, k+1)
#   dim(u_vec) <- c(length(x), k+1)
#   
#   res <- psi_k - mean_xpsi_k*u_vec/as.numeric(Q_m)
#   return(res)
# }
# 
# psi_uni_k <- function(x, k, gamma){
#   u <- x
#   p <- orthopolynom::hermite.he.polynomials(k, normalized = FALSE)
#   pevalx <- orthopolynom::polynomial.values(p, u)
#   pevalx <- unlist(pevalx)
#   dim(pevalx)<- c(length(x), k+1)
#   
#   wfx <- exp(-u^2*gamma/(2*(1 + gamma)))
#   wfx <- rep(wfx, k+1)
#   dim(wfx)<- c(length(x), k+1)
#   
#   fct <- factorial(0:k)
#   fct <- rep(fct, length(x))
#   dim(fct) <- c(k+1, length(x))
#   fct <- t(fct)
#   
#   res <- ((1.0 - gamma^2)^(1/4))*pevalx*wfx/sqrt(fct)
#   return(res)
# }
# 
# lin_test_uni_boot <- function(x, y, n_boot = 100, k_grid = 500){
#   
#   x <- scale(x, scale = FALSE)
#   y <- scale(y, scale = FALSE)
#   n_obs <- length(x)
#   
#   beta <- solve(t(x) %*% x)%*%(t(x)%*%y)
#   resid <- y - x%*%beta
#   
#   # u <- rnorm(k_grid, mean = 0, sd = 1)
#   u <- qnorm(seq(0.001, 0.999, length = k_grid))
#   
#   u_mat <- rep(u, n_obs)
#   dim(u_mat) <- c(k_grid, n_obs)
#   
#   x_mat <- rep(x, k_grid)
#   dim(x_mat) <- c(n_obs, k_grid)
#   
#   x_mat <- t(x_mat)
#   x_mat_sd1 <- x_mat/sd(x)
#   
#   resid_mat <- rep(resid, k_grid)
#   dim(resid_mat) <- c(n_obs, k_grid)
#   resid_mat <- t(resid_mat)
#   scale_x <- (t(x)%*%x)/(length(x))
#   
#   # compute the integrand
#   avg_0 <- apply(x_mat*exp(1i*u_mat*x_mat_sd1), 1, 'mean')/c(scale_x)
#   avg_0 <- rep(avg_0, n_obs)
#   dim(avg_0) <- c(k_grid, n_obs)
#   psi_mat <- resid_mat*(exp(1i*u_mat*x_mat_sd1) - avg_0*x_mat)
# 
#   # stat_boot <- w_mean <- rep(NA, n_boot)
#   # for(i in 1:n_boot){
#   #   delta <- rnorm(n_obs, mean = 0, sd = 1)
#   #   w_i_n <- (psi_mat %*% delta)/sqrt(n_obs)
#   #   # center the paths under the null, otherwise the simulation does not work! 
#   #   w_i_n <- w_i_n - mean(w_i_n)
#   #   stat_boot[i] <- mean(abs(w_i_n)^2)
#   #   w_mean[i] <- mean(w_i_n)
#   # }
#   
#   w_i_n <- psi_mat %*% rep(1/sqrt(n_obs), n_obs)
#   stat <- mean(abs(w_i_n)^2)
#   
#   delta <- matrix(rnorm(n_obs*n_boot), nrow = n_obs, ncol = n_boot)
#   w_i_n <- (psi_mat %*% delta)/sqrt(n_obs)
#   w_i_n <- scale(w_i_n, scale = FALSE) # center the paths under the null, otherwise the simulation does not work! 
#   stat_boot <- colMeans(apply(w_i_n, 2, function(x){abs(x)^2}))
#   
#   return(list(p_value = mean(stat_boot >= stat), prob = mean(stat_boot < stat), stat = stat, stat_boot = stat_boot))
# }
# 
# lin_test_uni_bierens <- function(x, y, n_boot = 100, k_grid = 500, tau = pi/2){
#   
#   x <- scale(x, scale = FALSE)
#   y <- scale(y, scale = FALSE)
#   n_obs <- length(x)
#   
#   beta <- solve(t(x) %*% x)%*%(t(x)%*%y)
#   resid <- y - x%*%beta
#   
#   u <- qunif(seq(0.001, 0.999, length = k_grid), min = -tau, max = tau)
#   
#   u_mat <- rep(u, n_obs)
#   dim(u_mat) <- c(k_grid, n_obs)
#   
#   x_mat <- rep(x, k_grid)
#   dim(x_mat) <- c(n_obs, k_grid)
#   
#   x_mat <- t(x_mat)
#   # transform addtionally with atan, the support of x_std is [-pi/2, pi/2]
#   x_mat_sd1 <- atan(x_mat/sd(x))
#   
#   resid_mat <- rep(resid, k_grid)
#   dim(resid_mat) <- c(n_obs, k_grid)
#   resid_mat <- t(resid_mat)
#   scale_x <- (t(x)%*%x)/(length(x))
#   
#   # compute the integrand
#   avg_0 <- apply(x_mat*exp(1i*u_mat*x_mat_sd1), 1, 'mean')/c(scale_x)
#   avg_0 <- rep(avg_0, n_obs)
#   dim(avg_0) <- c(k_grid, n_obs)
#   psi_mat <- resid_mat*(exp(1i*u_mat*x_mat_sd1) - avg_0*x_mat)
# 
#   w_i_n <- psi_mat %*% rep(1/sqrt(n_obs), n_obs)
#   stat <- mean(abs(w_i_n)^2)
#   
#   delta <- matrix(rnorm(n_obs*n_boot), nrow = n_obs, ncol = n_boot)
#   w_i_n <- (psi_mat %*% delta)/sqrt(n_obs)
#   w_i_n <- scale(w_i_n, scale = FALSE) # center the paths under the null, otherwise the simulation does not work! 
#   stat_boot <- colMeans(apply(w_i_n, 2, function(x){abs(x)^2}))
#   
#   return(list(p_value = mean(stat_boot >= stat), prob = mean(stat_boot < stat), stat = stat, stat_boot = stat_boot))
# }
# 
# lin_test_mvt <- function(x, y, k = 20, gamma = (sqrt(5)-1)/2, export_cdf = FALSE){
#   nobs <- nrow(x)
#   
#   Sig <- cov(x)
#   eig_out <- eigen(Sig)
#   eig_out$values <- abs(eig_out$values) # some values might equal -eps, virtually zero, so take abs
#   Sig_sqrt <- eig_out$vectors%*%diag(sqrt(eig_out$values))%*%t(eig_out$vectors)
#   
#   x_std <- x%*%solve(Sig_sqrt)
#   psi <- Psi_mvt_mat(x_std)
#   
#   x <- scale(x, scale = FALSE)
#   y <- scale(y, scale = FALSE)
#   beta <- solve(t(x) %*%x)%*%(t(x)%*%y)
#   resid <- y - x%*%beta
#   
#   stat <- t(resid)%*%psi%*%resid/nobs
#   
#   sig_res <- sd(resid)
#   
#   delta_k <- delta_mvt_k_ols(x, k, x_std, gamma)
#   
#   C <- (sig_res^2)*cov(delta_k$f_val)
#   # C <- diag(sqrt(delta_k$lambdas))%*%C%*%diag(sqrt(delta_k$lambdas))
#   rescale_cov(C, delta_k$lambdas) # this function call does not return a value but modifies C directly
#   eig <- eigen(C)
#   lambdas <- abs(eig$values[1:min(nobs, nrow(C))])
# 
#   h <- 0.005
#   K <- 1 + ceiling(log(2*sum(lambdas)/(h*lambdas[1]))/log(2))
#   K <- max(K, 13)
#   
#   out <- pQ(stat, lambdas, h = h, N = 2^K)
#   
#   if(export_cdf){
#     cdf <- pQ_fft_cpp(lambdas, h = h, N = 2^K)
#   }else{
#     cdf = NULL
#   }
#   
#   return(list(stat = stat, p_value = 1 - out, prob = out, cdf = cdf))
# }
# 
# 
# lin_test_mvt_diag <- function(x, y, k = 20, gamma = (sqrt(5)-1)/2, export_cdf = FALSE){
#   nobs <- nrow(x)
#   
#   Sig <- cov(x)
#   eig_out <- eigen(Sig)
#   eig_out$values <- abs(eig_out$values) # some values might equal -eps, virtually zero, so take abs
#   Sig_sqrt <- eig_out$vectors%*%diag(sqrt(eig_out$values))%*%t(eig_out$vectors)
#   
#   x_std <- x%*%solve(Sig_sqrt)
#   psi <- Psi_mvt_mat(x_std)
#   
#   x <- scale(x, scale = FALSE)
#   y <- scale(y, scale = FALSE)
#   beta <- solve(t(x) %*%x)%*%(t(x)%*%y)
#   resid <- y - x%*%beta
#   
#   stat <- t(resid)%*%psi%*%resid/nobs
#   
#   sig_res <- sd(resid)
#   
#   delta_k <- delta_mvt_k_ols(x, k, x_std, gamma)
#   n_terms <- ncol(delta_k$f_val)
#   
#   C_rank <- min(n_terms, nobs)
#   C <- (sig_res^2)*cov(delta_k$f_val[,1:C_rank])
#   
#   rescale_cov(C, delta_k$lambdas[1:C_rank]) # this function call does not return a value but modifies C directly
#   
#   # eig <- eigen(C_lambda)
#   eig <- eigen(C)
#   coeffs_rank <- abs(eig$values)
#   if(n_terms <= nobs){
#     coeffs <- coeffs_rank
#   }else{
#     coeffs_rest <- (sig_res^2)*apply(delta_k$f_val[,(C_rank + 1):n_terms], 2, 'var')*delta_k$lambdas[(C_rank + 1):n_terms]
#     coeffs <- c(coeffs_rank, coeffs_rest)
#   }
# 
#   h <- 0.005
#   K <- 1 + ceiling(log(2*sum(coeffs)/(h*coeffs[1]))/log(2))
#   K <- max(K, 13)
#   
#   out <- pQ(stat, coeffs, h = h, N = 2^K)
#   
#   if(export_cdf){
#     cdf <- pQ_fft_cpp(coeffs, h = h, N = 2^K)
#   }else{
#     cdf = NULL
#   }
#   
#   return(list(stat = stat, p_value = 1 - out, prob = out, cdf = cdf))
# }
# 
# 
# delta_mvt_k_ols <- function(x, k, x_std, gamma){
#   d <- ncol(x)
#   nobs <- nrow(x)
#   
#   A_m <- (t(x)%*%x)/(nrow(x))
#   
#   v_all <- rep(0, d)
#   lambdas <- 1
#   for(i in 1:k){
#     v <- tuples(d, i)
#     v_all <- rbind(v_all, v)
#     lambdas <- c(lambdas, rep(gamma^i, nrow(v)))
#   }
#   
#   # psi_x <- sample_psi_he_k_cpp(x_std, v_all, gamma)
#   psi_x <- sample_psi_he_k(x_std, t(v_all), gamma)
#   
#   mean_xpsi_k <- t(x)%*%psi_x/nobs
#   u_vec <- t(x)
#   res <- psi_x - t(t(mean_xpsi_k)%*%solve(A_m)%*%u_vec)
#   
#   return(list(f_val = res, lambdas = lambdas))
# }
# 
# pQ <- function(arg, lambdas, h = 0.005, N = 2^13){
#   # out <- pQ_fft(lambdas, h, N)
#   out <- pQ_fft_cpp(lambdas, h, N)
#   indexout <- which(arg > max(out$x))
#   indexin <- which(arg <= max(out$x))
#   
#   if(length(indexout) == 0){
#     prob <- approx(x = out$x, y = out$cdf, xout = arg, method = 'linear')[[2]]
#   }else{
#     prob <- rep(1, length(arg))
#     prob[indexin] <- approx(x = out$x, y = out$cdf, xout = arg[indexin], method = 'linear')[[2]]
#   }    
#   return(prob)
# }
# 
# pQ_fft <- function(lambdas, h = 0.005, N = 2^13){
#   
#   lam_1 <- lambdas[1]
#   lambdas <- lambdas/lam_1
#   
#   s <- 1/(h*N)
#   t1 <- 1:N
#   t2 <- 2*pi*(t1 - 1 - N/2)*s
#   
#   cfvalues <- Q_chf(t2, lambdas)
#   x1 <- (-1)^(t1 - 1)*cfvalues
#   
#   # compute pdf values
#   pdf <- Re(fft(x1))*(s*(-1)^(t1 - 1 - N/2))
#   pdf <- pdf/lam_1
#   
#   x <- (t1 - 1 - N/2)*h*lam_1
#   cdf <- cumsum(abs(pdf))/sum(abs(pdf)) # needs correction because of oscillation artifacts
#   cdf <- cdf[x >= 0]
#   x <- x[x >= 0]
#   return(list(x = x, cdf = cdf))
# }
# 
# 
# 
# Q_chf <- function(arg, lambdas){
#   arg_n <- length(arg)
#   la_k <- length(lambdas)
#   arg <- rep(arg, la_k)
#   dim(arg) <- c(arg_n, la_k)
#   arg <- t(arg)
#   
#   lambdas <- rep(lambdas, arg_n)
#   dim(lambdas) <- c(la_k, arg_n)
#   
#   chf <- (1 - 2*1i*arg*lambdas)^(-1/2)
#   chf <- apply(chf, 2, 'prod')
#   return(chf)
# }
# 
# 
# norm_vec <- function(x){return(sqrt(sum(x^2)))}
# 
# 
# d_tuple_k_sum <- function(x_dim, k){
#   v_all <- rep(0, x_dim)
#   for(i in 1:k){
#     v <- tuples(x_dim, i)
#     v_all <- rbind(v_all, v)
#   }
#   return(v_all)
# }
# 
# 
# tuples <- function(size, total, tmp_tuple = NULL){
#   if(size == 1){
#     tmp_tuple[1] <- total
#     return(tmp_tuple)
#   }
#   if(is.null(tmp_tuple)){
#     tmp_tuple <- vector(mode = 'double', length = size)
#   }
#   
#   res <- NULL
#   for(i in 0:total){
#     tmp_tuple[size] <- i
#     out <- tuples(size - 1, total - i, tmp_tuple)
#     res <- rbind(res, out)
#   }
#   return(res)
# }
# 
# 
# 
# 
# 
# 
# 
# lin_test_mvt_boot <- function(x, y, n_boot = 100, k_grid = 500, u_sd = 1){
#   
#   x <- scale(x, scale = FALSE)
#   y <- scale(y, scale = FALSE)
#   n_obs <- nrow(x)
#   d_dim <- ncol(x)
#   
#   beta <- solve(t(x) %*% x)%*%(t(x)%*%y)
#   resid <- y - x%*%beta
#   
#   Sig <- cov(x)
#   eig_out <- eigen(Sig)
#   eig_out$values <- abs(eig_out$values) # some values might equal -eps, virtually zero, so take abs
#   Sig_sqrt <- eig_out$vectors%*%diag(sqrt(eig_out$values))%*%t(eig_out$vectors)
#   x_std <- x%*%solve(Sig_sqrt)
#   
#   resid_mat <- rep(resid, k_grid)
#   dim(resid_mat) <- c(n_obs, k_grid)
#   # scale_x <- (t(x)%*%x)/(length(x))
#   scale_x <- (t(x)%*%x)/n_obs
#   
#   # compute the rhs
#   scale_x_inv <- solve(scale_x)
#   rhs <- matrix(nrow = n_obs, ncol = k_grid)
#   for(i in 1:k_grid){
#     u <- rnorm(d_dim, mean = 0, sd = u_sd)
#     exp_ux <- exp(1i*(x_std %*% u))
#     exp_ux_mat <- rep(exp_ux, d_dim)
#     dim(exp_ux_mat) <- c(n_obs, d_dim)
#     rhs[,i] <- (exp_ux - x %*% scale_x_inv %*% colMeans(exp_ux_mat * x))*resid/sqrt(n_obs)
#   }
#   
#   w_i_n <- colSums(rhs)
#   stat <- mean(abs(w_i_n)^2)
#   
#   delta <- matrix(rnorm(n_obs*n_boot), nrow = n_boot, ncol = n_obs)
#   w_i_n <- t(delta %*% rhs)
#   w_i_n <- scale(w_i_n, scale = FALSE)
#   stat_boot <- colMeans(abs(w_i_n)^2)
#   
#   return(list(p_value = mean(stat_boot >= stat), prob = mean(stat_boot <= stat), stat = stat, stat_boot = stat_boot))
# }
# 
# 
# 
# lin_test_mvt_bierens <- function(x, y, n_boot = 100, k_grid = 500, tau = pi/2){
#   
#   x <- scale(x, scale = FALSE)
#   y <- scale(y, scale = FALSE)
#   n_obs <- nrow(x)
#   d_dim <- ncol(x)
#   
#   beta <- solve(t(x) %*% x)%*%(t(x)%*%y)
#   resid <- y - x%*%beta
#   
#   Sig <- cov(x)
#   eig_out <- eigen(Sig)
#   eig_out$values <- abs(eig_out$values) # some values might equal -eps, virtually zero, so take abs
#   Sig_sqrt <- eig_out$vectors%*%diag(sqrt(eig_out$values))%*%t(eig_out$vectors)
#   
#   # transform addtionally with atan, the support of x_std is [-pi/2, pi/2]
#   x_std <- atan(t(x%*%solve(Sig_sqrt)))
#   
#   u <- matrix(runif(k_grid*d_dim, min = -tau, max = tau), nrow = d_dim, ncol = k_grid)
# 
#   resid_mat <- rep(resid, k_grid)
#   dim(resid_mat) <- c(n_obs, k_grid)
#   # scale_x <- (t(x)%*%x)/(length(x))
#   scale_x <- (t(x)%*%x)/(n_obs)
#   
#   # compute the rhs
#   scale_x_inv <- solve(scale_x)
#   mu_1_mat <- matrix(nrow = d_dim, ncol = k_grid)
#   exp_ux <- matrix(nrow = n_obs, ncol = k_grid)
#   rhs <- matrix(nrow = n_obs, ncol = k_grid)
#   for(i in 1:k_grid){
#     exp_ux[,i] <- exp(1i*(t(u[,i]) %*% x_std))
#     exp_ux_mat <- rep(exp_ux[,i], d_dim)
#     dim(exp_ux_mat) <- c(n_obs, d_dim)
#     mu_1_mat[,i] <- colMeans(exp_ux_mat * x)
#     rhs[,i] <- exp_ux[,i] - x %*% scale_x_inv %*% mu_1_mat[,i]
#   }
#   
#   rhs <- resid_mat * rhs
#   
#   w_i_n <- rep(1/sqrt(n_obs), n_obs) %*% rhs
#   # stat <- mean(abs(w_i_n)^2)
#   stat <- bierens_stat_cpp(resid, t(x_std), tau)
#   
#   delta <- matrix(rnorm(n_obs*n_boot), nrow = n_boot, ncol = n_obs)
#   # rhs <- scale(t(rhs), scale = FALSE)
#   # w_i_n <- t(delta %*% t(rhs))/sqrt(n_obs)
#   w_i_n <- t(delta %*% rhs)/sqrt(n_obs)
#   w_i_n <- scale(w_i_n, scale = FALSE)
#   stat_boot <- colMeans(abs(w_i_n)^2)
#   
#   return(list(p_value = mean(stat_boot >= stat), prob = mean(stat_boot <= stat), 
#               stat = stat, stat_boot = stat_boot))
# }


# bierens_stat <- function(resid, x_transformed, tau){
#   p_j1_j2 <- function(Z_j1, Z_j2, tau){
#     delta_Z <- Z_j1 - Z_j2
#     res <-  1
#     for(i in 1:length(delta_Z)){
#       if((delta_Z[i]) != 0){
#         res <- res * sin(tau*delta_Z[i])/(tau*delta_Z[i])
#       }
#     }
#     return(res)
#   }
#   
#   int_sum <- 0
#   n_obs <- nrow(x_transformed)
#   for(j1 in 1:(n_obs-1)){
#     for(j2 in (j1 + 1):n_obs){
#       int_sum <- int_sum + resid[j1]*resid[j2]*p_j1_j2(x_transformed[j1, ], 
#                                                        x_transformed[j2, ], tau)
#     }
#   }
#   int_sum <- int_sum*2/n_obs + sum(resid^2)/n_obs
#   return(int_sum)
# }

