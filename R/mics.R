#' Testing cell-type-specific mediation effects from bulk methylation data
#'
#' @param meth_data The observed methylation matrix. Rows represent CpG sites and columns are samples.
#' Its CpG site names are required to align with the CpG site names of the methylation reference when estimating cellular compositions.
#' @param meth_ref The methylation reference matrix. Rows represent CpG sites and columns are cell type annotations.
#' Its CpG site names are required to align with the CpG site names of the observed methylation matrix when estimating cellular compositions.
#' @param thd The threshold used to filter out cell types with a small proportion. After estimating the cellular composition for each sample,
#' we filter out the cell types with its proportions' 0.75 quantile across all samples less than the threshold. The default is 0.03.
#' @param S The exposure vector. Each component corresponds to a sample.
#' @param X The confounding factor matrix. Rows are samples and columns represent confounding factors.
#' Note that if there is only one confounding factor, X also needs to be in a matrix form and has only one column.
#' @param Y The outcome vector. Each component corresponds to a sample.
#'
#' @return A list containing the following items. \item{pval_mediator}{The cell-type-specific p-value matrix by regressing the mediator to the exposure, where rows represent CpG sites and columns are cell types.}
#' \item{pval_outcome}{The cell-type-specific p-value matrix by regressing the outcome to the mediator and the exposure, where rows represent CpG sites and columns are cell types.}
#' \item{pval_joint_sq}{The joint cell-type-specific p-value matrix by maximizing and squaring, where rows represent CpG sites and columns are cell types.}
#' \item{P_matr}{The cell proportion matrix, where rows are cell types and columns correspond to samples.}
#'
#' @examples
#' rm(list = ls())
#' ###########################################################
#' #Generate data used in the simulation study
#' ###########################################################
#' set.seed(03292019)
#' #Generate samples from the Dirichlet distribution
#' rDirichlet <- function(alpha_vec){
#'  num <- length(alpha_vec)
#'  temp <- rgamma(num, shape = alpha_vec, rate = 1)
#'  return(temp / sum(temp))
#' }
#' n <- 400     #number of samples
#' m <- 10000   #number of CpG sites
#' K <- 5       #underlying cell type number
#'
#' #Generate the baseline methylation profiles mu
#' #cell type 1
#' methy1 <- rbeta(m,3,3)
#'
#' #cell type 2
#' methy2 <- methy1 + rnorm(m, sd=0.01)
#' ind <- sample(1:m, m/10)
#' methy2[ind] <- rbeta(length(ind),3,3)
#'
#' #cell type 3
#' methy3 <- methy1 + rnorm(m, sd=0.01)
#' ind <- sample(1:m, m/10)
#' methy3[ind] <- rbeta(length(ind),3,3)
#'
#' #cell type 4
#' methy4 <- methy1 + rnorm(m, sd=0.01)
#' ind <- sample(1:m, m/10)
#' methy4[ind] <- rbeta(length(ind),3,3)
#'
#' #cell type 5
#' methy5 <- methy1 + rnorm(m, sd=0.01)
#' ind <- sample(1:m, m/10)
#' methy5[ind] <- rbeta(length(ind),3,3)
#'
#' mu_matr <- cbind(methy1, methy2, methy3, methy4, methy5)
#'
#' #Generate the exposure, confounders, and the outcome
#' S <- rnorm(n) #the exposure
#' D <- 2        #the confounder number
#' X <- matrix(rnorm(D*n), n, D) #confounders
#'
#' Y <- 0.5*S + 0.5*X[ ,1] + 0.5*X[ ,2] + rnorm(n, sd = 0.5)
#'
#' #Generate the coefficients in Part I
#' #the effects of the exposure on the mediators
#' beta_matr <- array(0, dim=c(m,K))
#'
#' m_common <- 10
#' max_signal <- 0.15
#' min_signal <- 0.07
#'
#' signs <- sample(c(-1,1), m_common*K, replace=T)
#' beta_matr[1:m_common,1:K] <- signs * runif(m_common*K, min=min_signal, max=max_signal)
#'
#' m_seperate <- 10
#' signs <- sample(c(-1,1), m_seperate*2, replace=T)
#' beta_matr[m_common+(1:m_seperate),1:2] <- signs * runif(m_seperate*2, min=min_signal, max=max_signal)
#'
#' for(k in 3:K){
#'   signs <- sample(c(-1,1), m_seperate, replace=T)
#'  beta_matr[m_common+(k-2)*m_seperate+(1:m_seperate),k] <- signs * runif(m_seperate, min=min_signal, max=max_signal)
#' }
#'
#' #the effects of the confounders on the mediators
#' h_vec1 <- rnorm(m, sd=0.01)
#' h_vec2 <- rnorm(m, sd=0.01)
#'
#' #generate the cellular compositions
#' P_matr <- sapply(1:n, function(i){
#'   rDirichlet(rep(2, K))
#' })
#'
#' #Generate the effect of mediators on the outcome
#' gamma_tilde_matr <- matrix(0, m, K)
#' for(k in 1:K){
#'   signs <- sample(c(-1,1), m_seperate, replace=T)
#'   gamma_tilde_matr[(k-1)*m_seperate+1:m_seperate, k] <- signs * runif(m_seperate, min=0.03, max=0.1)
#'   if(k == 3){
#'     signs <- sample(c(-1,1), m_seperate, replace=T)
#'     gamma_tilde_matr[(k-3)*m_seperate+1:m_seperate, k] <- signs * runif(m_seperate, min=0.03, max=0.1)
#'   }
#' }
#'
#' #Generate u conditional on Y
#' u <- array(NA, dim=c(n, m, K))
#' for(i in 1:n){
#'   epsilon_matr <- matrix(rnorm(m*K, sd = 0.015), m, K)
#'   #epsilon_matr <- matrix(c(rnorm(m*2, sd = 0.02), rnorm(m*2, sd = 0.01), rnorm(m, sd = 0.015)), m, K)
#'   u[i,,] <- mu_matr + beta_matr*S[i] + h_vec1*X[i,1] + h_vec2*X[i,2] + gamma_tilde_matr * Y[i] + epsilon_matr
#' }
#'
#' #Generate the observed methylation profiles
#' Ometh <- NULL
#' for(i in 1:n){
#'   tmp <- u[i,,] %*% P_matr[ ,i]
#'   Ometh <- cbind(Ometh, tmp)
#' }
#'
#' sum(Ometh > 1)
#' Ometh[Ometh > 1] <- 1
#'
#' sum(Ometh < 0)
#' Ometh[Ometh < 0] <- 0
#'
#' rownames(Ometh) <- paste0("cpg", 1:m)
#' rownames(mu_matr) <- paste0("cpg", 1:m)
#' colnames(mu_matr) <- paste0("cell_type", 1:K)
#'
#' ###########################################################
#' #Conduct MICS
#' ###########################################################
#'
#' t1 <- Sys.time()
#' out <- mics(Ometh, mu_matr, S, X, Y)
#' pval_joint_sq <- out$pval_joint_sq
#' t2 <- Sys.time()
#' print(t2 - t1)




mics <- function(meth_data, meth_ref, S, X, Y, thd = 0.03){
	############################################
	# step 1. preparation
	############################################
	Ometh <- meth_data
	ref_meth <- meth_ref

	#covert M values to beta values if Ometh collects M value.
	if(max(Ometh) > 1){
		Ometh <- 2^(Ometh) / (1+2^(Ometh))
	}

	#number of CpG sites
	m <- nrow(Ometh)

	#number of samples
	n <- ncol(Ometh)

	#number of cell types
	K <- ncol(ref_meth)

	#find the common CpG sites between Ometh and ref_meth
	#only use the common CpG sites to estimate cell proportions
	ind <- match(rownames(Ometh), rownames(ref_meth))

	#the boolean vector indicating which one is not NA
	#In this way, Ometh_comm and ref_meth_comm have the same CpG site name for each row.
	no_na <- !is.na(ind)
	Ometh_comm <- Ometh[no_na,]
	ref_meth_comm <- ref_meth[ind[no_na], ]


	##############################################
	# step 2. estimate cellular compositions
	##############################################

	#we need to use "solve.QP" function in R package "quadprog"
	#use ?solve.QP to see the following notations
	Dmat <- 2*t(ref_meth_comm) %*% ref_meth_comm

	Amat <- cbind(rep(1, K), diag(rep(1,K)))
	bvec <- c(1, rep(0, K))

	#estimate the cell proportion P matrix
	# cell types in row and samples in column
	P_matr <- sapply(1:n, function(i){
				dvec <- 2 * t(ref_meth_comm) %*% as.numeric(Ometh_comm[ ,i])
				solu <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
				#meq = 1 means that the first constraint is treated as an equal contraint
				solu$solution
			})

	#filter out rare cell types
	del_celltype <- which(apply(P_matr, 1, quantile, 0.75) < thd)

	if(length(del_celltype) > 0){
		P_matr <- P_matr[-del_celltype, ]
		P_matr <- P_matr / colSums(P_matr)
		K <- K - length(del_celltype)
		ref_meth_comm <- ref_meth_comm[ ,-del_celltype]
	}

	##############################################
	# step 3. estimate the effect of S on
	#     methylation adjusted for X
	##############################################

	#construct the design matrix x_matr (samples in rows)
	x_matr <- cbind(t(P_matr), t(P_matr)*S, X)
	x_matr <- as.matrix(x_matr)

	#x_matr divided by square root of the sum of squares
	sq_P2 <- sqrt(colSums(P_matr^2))
	x_matr <- x_matr / sq_P2

	beta_pval <- NULL

	#obtain the estimates
	for(j in 1:m){
		y_vec <- Ometh[j,] / sq_P2

		#linear regression WITHOUT the intercept
		fit.m <- lm(y_vec~-1+x_matr)
		fit.m.sum <- summary(fit.m)

		#the first K estimates are \mu_j
		#the second K estimates are \beta_j
		beta_pval <- rbind(beta_pval, fit.m.sum$coef[K+(1:K),4])
	}

	colnames(beta_pval) <- colnames(ref_meth_comm)
	rownames(beta_pval) <- rownames(Ometh)

	##############################################
	# step 4. estimate the effects of
	# methylation on Y adjusted for S and X
	##############################################

	#construct the design matrix x_matr (samples in rows)
	x_matr2 <- cbind(t(P_matr), t(P_matr)*Y, t(P_matr)*S)
	for(ell in 1:ncol(X)){
		x_matr2 <- cbind(x_matr2, t(P_matr)*X[,ell])
	}

	x_matr2 <- as.matrix(x_matr2)

	#x_matr2 divided by square root of the sum of squares
	x_matr2 <- x_matr2 / sq_P2

	gamma_tilde_pval <- NULL

	#obtain the estimates
	for(j in 1:m){
		y_vec <- Ometh[j,] / sq_P2

		#linear regression WITHOUT the intercept
		fit.m <- lm(y_vec~-1+x_matr2)
		fit.m.sum <- summary(fit.m)

		#the first K estimates are \nu_j
		#the second K estimates are \gamma_tilde_j
		gamma_tilde_pval <- rbind(gamma_tilde_pval, fit.m.sum$coef[K+(1:K),4])
	}

	colnames(gamma_tilde_pval) <- colnames(ref_meth_comm)
	rownames(gamma_tilde_pval) <- rownames(Ometh)

	########################################################
	# step 5. Joint significance method followed by squaring
	########################################################
	#pval_mediator
	pval_mediator <- beta_pval
	colnames(pval_mediator) <- colnames(ref_meth_comm)
	rownames(pval_mediator) <- rownames(Ometh)

  	#pval_outcome
	pval_outcome <- gamma_tilde_pval
	colnames(pval_outcome) <- colnames(ref_meth_comm)
	rownames(pval_outcome) <- rownames(Ometh)

	#pval_joint_sq
	pval_joint <- matrix(NA, m, K)
	for(j in 1:m){
	  for(k in 1:K){
	    pval_joint[j,k] <- max(beta_pval[j,k], gamma_tilde_pval[j,k])
	  }
	}

	pval_joint_sq <- pval_joint^2
	colnames(pval_joint_sq) <- colnames(ref_meth_comm)
	rownames(pval_joint_sq) <- rownames(Ometh)


	#return values
  	out <- list()
  	out$pval_mediator <- pval_mediator
  	out$pval_outcome <- pval_outcome
  	out$pval_joint_sq <- pval_joint_sq

 	rownames(P_matr) <- colnames(ref_meth_comm)
  	out$P_matr <- P_matr
	return(out)
}
