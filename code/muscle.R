library(snow)

one_rep <- function(new_params, current_params) {
  trash <- tryCatch({
    source("../code/datamaker_only_counts.R")
    source("../code/fit_methods.R")
    args_val <- append(current_params, new_params)
    set.seed(new_params$current_seed)

    d_out <- datamaker_counts_only(args_val)

    which_null <- d_out$meta$null

    half_null <- which_null
    half_null[half_null == 1][sample(1:sum(which_null), size = sum(which_null) / 2)] <- 0

    beta_true <- rep(0, length = args_val$Ngene)
    beta_true[!which_null] <- d_out$meta$true_log2foldchange

    X <- as.matrix(model.matrix(~d_out$input$condition))
    colnames(X) <- c("Intercept", "Treatment")
    Y <- t(log2(as.matrix(d_out$input$counts + 1)))

    num_sv <- sva::num.sv(t(Y), mod = X, method = "be")



    fit_out <- fit_all_competitors(Y = Y, X = X, num_sv = num_sv,
                                   control_genes = half_null)

    mse_vec <- colMeans((fit_out$betahat_df[half_null == 0, ] - beta_true[half_null == 0]) ^ 2)
    pi0_vec <- fit_out$pi0hat_vec
    auc_vec <- rep(NA, length = length(pi0_vec))
    if(args_val$nullpi < 1) {
        for (col_index in 1:ncol(fit_out$lfdr)) {
            if (all(is.na(fit_out$lfdr[, col_index]))) {
                auc_vec[col_index] <- NA
            } else {
                auc_vec[col_index] <- get_auc(predictor = c(fit_out$lfdr[half_null == 0, col_index]),
                                              response = which_null[half_null == 0])
            }
        }
    }

    return_vec <- c(pi0_vec, mse_vec, auc_vec)

    method_names <- colnames(fit_out$betahat_df)

    names(return_vec) <- c(paste0("pi0hat_", method_names),
                           paste0("mse_", method_names),
                           paste0("auc_", method_names))
    return_list <- list(return_vec = return_vec, lfsr_gamma = fit_out$biash_lfsr)
    TRUE
  },
  error = function(e) {NULL})
  if(is.null(trash)) {
    return_vec <- rep(NA, length = 30)
    return_list <- list(return_vec = return_vec, lfsr_gamma = NA)
  }
    return(return_list)
}


itermax <- 100

## these change
Nsamp_seq  <- c(5, 10, 20)
nullpi_seq <- c(0.5, 0.9, 1)

par_vals <- expand.grid(list(1:itermax, Nsamp_seq, nullpi_seq))
colnames(par_vals) <- c("current_seed", "Nsamp", "nullpi")
par_vals$poisthin <- TRUE
par_vals$poisthin[abs(par_vals$nullpi - 1) < 10 ^ -10] <- FALSE

par_list <- list()
for (list_index in 1:nrow(par_vals)) {
    par_list[[list_index]] <- list()
    for (inner_list_index in 1:ncol(par_vals)) {
        par_list[[list_index]][[inner_list_index]] <- par_vals[list_index, inner_list_index]
        names(par_list[[list_index]])[inner_list_index] <- colnames(par_vals)[inner_list_index]
    }
}

## these do not change
args_val              <- list()
args_val$log2foldsd   <- 1
args_val$tissue       <- "muscle"
args_val$path         <- "~/Dropbox/Gene_Set_GTEx/"
args_val$Ngene        <- 1000
args_val$log2foldmean <- 0
args_val$skip_gene    <- 0

## ## If on your own computer, use this
library(parallel)
cl <- makeCluster(detectCores()-1)
soutcomplete <- t(parSapply(cl = cl, par_list, FUN = one_rep, current_params = args_val))
stopCluster(cl)

## ## on RCC, use this
## np <- mpi.universe.size() - 1
## cluster <- makeMPIcluster(np)
## sout <- t(snow::parSapply(cl = cluster, X = par_list, FUN = one_rep, current_params = args_val))
## stopCluster(cluster)
## mpi.exit()

sout <- t(sapply(soutcomplete[, 1], c))

num_conds <- ncol(sout) / 3
pi0_mat <- cbind(par_vals, sout[, 1:num_conds])
mse_mat <- cbind(par_vals, sout[, (num_conds + 1):(2 * num_conds)])
auc_mat <- cbind(par_vals, sout[, (2 * num_conds + 1):(3 * num_conds)])

lfsr <- soutcomplete[, 2]
lfsr_mat = matrix(nrow = length(lfsr), ncol = max(sapply(lfsr, length)))
for (i in 1:nrow(lfsr_mat)){
	for (j in 1:length(lfsr[[i]])){
		lfsr_mat[i, j] = lfsr[[i]][j]
	}
}
lfsr_mat = cbind(par_vals, lfsr_mat)

write.csv(pi0_mat, file = "pi0_mat.csv", row.names = FALSE)
write.csv(mse_mat, file = "mse_mat.csv", row.names = FALSE)
write.csv(auc_mat, file = "auc_mat.csv", row.names = FALSE)
write.csv(lfsr_mat, file = "lfsr_mat.csv", row.names = FALSE)