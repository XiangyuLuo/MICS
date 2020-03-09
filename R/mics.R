#' Testing cell-type-specific mediation effects from bulk methylation data
#'
#' @param meth_data The observed methylation matrix. Rows represent CpG sites and columns are samples.
#' Its CpG site names are required in order to automatically align with the CpG site names of the methylation reference when estimating cellular compositions.
#' @param S The exposure vector. Each component corresponds to a sample.
#' @param X The confounding factor matrix. Rows are samples and columns represent confounding factors.
#' Note that if there is only one confounding factor, X also needs to be in a matrix form and has only one column.
#' @param Y The outcome vector. Each component corresponds to a sample.
#' @param cell_prop The cellular proportion matrix, where rows correspond to cell types and columns are samples. The default is NA, in which case the cell proportions are estimated from the methylation reference.
#' @param meth_ref The methylation reference matrix. If cell_prop is NA, meth_ref must be specified, where rows represent CpG sites and columns are cell type annotations.
#' Its CpG site names are required in order to automatically align with the CpG site names of the observed methylation matrix when estimating cellular compositions.
#' @param MCP.type The type of multiple comparison procedure when implementing Heller/Sampson method. It is either FWER or FDR. The default is FDR. 
#' @param cpg_ind The index vector of the informative cpg sites. The default is NA, which means that all cpg sites are used to estimate
#' the cellular compositions. Otherwise, only the informative cpg sites are used.
#' @param thd The threshold used to filter out cell types with a small proportion. After estimating the cellular composition for each sample,
#' we filter out the cell types with its proportions' 0.75 quantile across all samples less than the threshold. The default is 0, which means that no cell type would be removed.
#' @param maxp_sq A Boolean type argument indicating whether outputs combined p-values using maximum and squaring procedure. The default is FALSE.
#'
#' @return A list containing the following items. \item{pval_mediator}{The cell-type-specific p-value matrix by regressing the mediator to the exposure, where rows represent CpG sites and columns are cell types.}
#' \item{pval_outcome}{The cell-type-specific p-value matrix by regressing the outcome to the mediator and the exposure, where rows represent CpG sites and columns are cell types.}
#' \item{pval_joint_MultiMed}{The joint cell-type-specific p-value matrix using the Heller/Sampson procedure, where rows represent CpG sites and columns are cell types.}
#' \item{pval_joint_sq}{Only when the argument maxp_sq is TRUE, the returned value is the joint cell-type-specific p-value matrix by maximizing and squaring, where rows represent CpG sites and columns are cell types. Otherwise, it is NA.}
#' \item{P_matr}{The cell proportion matrix, where rows are cell types and columns correspond to samples.}
#'
#' @examples
#' rm(list = ls())
#' set.seed(03092020)
#' 
#' #set the working directory
#' ############################################################################
#' # Generate Data 
#' ############################################################################
#' #Generate samples from the Dirichlet distribution
#' rDirichlet <- function(alpha_vec){
#'   num <- length(alpha_vec)
#'   temp <- rgamma(num, shape = alpha_vec, rate = 1)
#'   return(temp / sum(temp))
#' }
#' 
#' n <- 400     #number of samples 
#' m <- 10000   #number of CpG sites
#' K <- 5       #underlying cell type number
#' 
#' ####################################################
#' #Generate the exposure, confounders, and the outcome
#' S <- rnorm(n) #the exposure
#' D <- 2        #the confounder number
#' X <- matrix(rnorm(D*n), n, D) #confounders
#' 
#' ##############################################	
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
#' colnames(mu_matr) <- c("celltype 1", "celltype 2", "celltype 3",
#'                        "celltype 4", "celltype 5")
#' 
#' ################################################
#' #Generate the coefficients
#' #the effects of the exposure on the mediators
#' beta_matr <- array(0, dim=c(m,K))
#' 
#' #the cpg sites affected by exposure for each cell type
#' cpg_type1 <- 1:75
#' cpg_type2 <- 26:100
#' cpg_type3 <- 51:125
#' cpg_type4 <- 76:150
#' cpg_type5 <- 101:175
#' 
#' #coefficients
#' max_signal <- 0.25
#' min_signal <- 0.1
#' 
#' signs <- sample(c(-1,1), 75, replace=T)
#' beta_matr[cpg_type1, 1] <- signs * runif(75, min=min_signal, max=max_signal)
#' 
#' signs <- sample(c(-1,1), 75, replace=T)
#' beta_matr[cpg_type2, 2] <- signs * runif(75, min=min_signal, max=max_signal)
#' 
#' signs <- sample(c(-1,1), 75, replace=T)
#' beta_matr[cpg_type3, 3] <- signs * runif(75, min=min_signal, max=max_signal)
#' 
#' signs <- sample(c(-1,1), 75, replace=T)
#' beta_matr[cpg_type4, 4] <- signs * runif(75, min=min_signal, max=max_signal)
#' 
#' signs <- sample(c(-1,1), 75, replace=T)
#' beta_matr[cpg_type5, 5] <- signs * runif(75, min=min_signal, max=max_signal)
#' 
#' #the effects of the confounders on the mediators
#' h_vec1 <- rnorm(m, sd=0.01)
#' h_vec2 <- rnorm(m, sd=0.01)
#' 
#' 
#' ##################################################
#' #Generate the effect of mediators on the outcome
#' gamma_matr <- matrix(0, m, K)
#' max_signal_gamma <- 0.035
#' min_signal_gamma <- 0.025
#' 
#' #the cpg sites affected by mediators for each cell type
#' cpg_type1_m <- 1:50
#' cpg_type2_m <- 26:75
#' cpg_type3_m <- 51:100
#' cpg_type4_m <- 76:125
#' cpg_type5_m <- 101:150
#' 
#' signs <- sample(c(-1,1), 50, replace=T)
#' gamma_matr[cpg_type1_m, 1] <- signs * runif(50, min=min_signal_gamma, max=max_signal_gamma)
#' 
#' signs <- sample(c(-1,1), 50, replace=T)
#' gamma_matr[cpg_type2_m, 2] <- signs * runif(50, min=min_signal_gamma, max=max_signal_gamma)
#' 
#' signs <- sample(c(-1,1), 50, replace=T)
#' gamma_matr[cpg_type3_m, 3] <- signs * runif(50, min=min_signal_gamma, max=max_signal_gamma)
#' 
#' signs <- sample(c(-1,1), 50, replace=T)
#' gamma_matr[cpg_type4_m, 4] <- signs * runif(50, min=min_signal_gamma, max=max_signal_gamma)
#' 
#' signs <- sample(c(-1,1), 50, replace=T)
#' gamma_matr[cpg_type5_m, 5] <- signs * runif(50, min=min_signal_gamma, max=max_signal_gamma)
#' 
#' 
#' #######################################################
#' #Generate the mediators and the outcome 
#' 
#' #the effects of the confounders on the outcome
#' eta1 <- 0.5
#' eta2 <- 0.5
#' 
#' #the effects of the exposure on the outcome
#' eta0 <- 0.5
#' 
#' #jointly generating u and Y is equivalent to first generating
#' # Y and then generating u conditional on Y
#' 
#' Y <- 0.5*S + 0.5*X[ ,1] + 0.5*X[ ,2] + rnorm(n, sd = 0.5)
#' 
#' 
#' #Generate u conditional on Y
#' u <- array(NA, dim=c(n, m, K))
#' for(i in 1:n){
#'   epsilon_matr <- matrix(rnorm(m*K, sd = 0.01), m, K)
#'   u[i,,] <- mu_matr + beta_matr*S[i] + h_vec1*X[i,1] + h_vec2*X[i,2] + 
#'   gamma_matr*Y[i] + epsilon_matr  
#' }
#' 
#' #############################################
#' #Generate the cellular compositions 
#' P_matr <- sapply(1:n, function(i){
#'   rDirichlet(c(4, 1, 3, 1, 1))
#' })
#' 
#' #############################################
#' #Generate the observed methylation profiles 
#' Ometh <- NULL
#' for(i in 1:n){
#'   tmp <- u[i,,] %*% P_matr[ ,i]
#'   Ometh <- cbind(Ometh, tmp)
#' }
#' 
#' #sum(Ometh > 1)
#' Ometh[Ometh > 1] <- 1
#' 
#' #sum(Ometh < 0)
#' Ometh[Ometh < 0] <- 0
#' 
#' rownames(Ometh) <- paste0("cpg", 1:m)
#' rownames(mu_matr) <- paste0("cpg", 1:m)
#' 
#' 
#' ############################################################################
#' # Conduct MICS 
#' ############################################################################
#' 
#' 
#' #install MICS package first
#' library(MICS)
#' 
#' out <- mics(meth_data = Ometh, S = S, X = X, Y = Y, meth_ref = mu_matr, 
#'             MCP.type = "FDR", maxp_sq = TRUE)
#' 
#' out$pval_joint_sq[1:100, 1]
#' 
#' 
#' pval_MultiMed <- out$pval_joint_MultiMed
#' K <- 5
#' m <- 10000
#' 
#' fdr_thred <- 0.2
#' 
#' #the true cpg sites through the mediation pathway 
#' cpg_type1 <- 1:50
#' cpg_type2 <- 26:75
#' cpg_type3 <- 51:100
#' cpg_type4 <- 76:125
#' cpg_type5 <- 101:150
#' 
#' ###cell type 1
#' #power 
#' sum(pval_MultiMed[cpg_type1, 1] < fdr_thred) / 50
#' #false positive rate
#' sum(pval_MultiMed[setdiff(1:m, cpg_type1), 1] < fdr_thred) / (m-50)
#' 
#' ###cell type 2
#' #power 
#' sum(pval_MultiMed[cpg_type2, 2] < fdr_thred) / 50
#' #false positive rate
#' sum(pval_MultiMed[setdiff(1:m, cpg_type2), 2] < fdr_thred) / (m-50)
#' 
#' ###cell type 3
#' #power 
#' sum(pval_MultiMed[cpg_type3, 3] < fdr_thred) / 50
#' #false positive rate
#' sum(pval_MultiMed[setdiff(1:m, cpg_type3), 3] < fdr_thred) / (m-50)
#' 
#' ###cell type 4
#' #power 
#' sum(pval_MultiMed[cpg_type4, 4] < fdr_thred) / 50
#' #false positive rate
#' sum(pval_MultiMed[setdiff(1:m, cpg_type4), 4] < fdr_thred) / (m-50)
#' 
#' ###cell type 5
#' #power 
#' sum(pval_MultiMed[cpg_type5, 5] < fdr_thred) / 50
#' #false positive rate
#' sum(pval_MultiMed[setdiff(1:m, cpg_type5), 5] < fdr_thred) / (m-50)
#' 
#' 


mics <- function(meth_data, S, X, Y, cell_prop = NA, meth_ref = NA, MCP.type = "FDR", cpg_ind = NA, thd = 0, maxp_sq = FALSE){
  if(is.na(cell_prop)[1]&is.na(meth_ref)[1]){
    stop("Either cell_prop or meth_ref should be specified!\n")
  }


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
	if(is.na(cell_prop)[1]){
	  K <- ncol(ref_meth)
	}else{
	  K <- nrow(cell_prop)
	}

	if(is.na(cell_prop)[1]){

	  #If the informative cpg sites are available
	  if(!is.na(cpg_ind)[1]){
      Ometh_tmp <- Ometh[cpg_ind, ]
	  }else{
	    Ometh_tmp <- Ometh
	  }
	  #find the common CpG sites between Ometh and ref_meth
	  #only use the common CpG sites to estimate cell proportions
	  ind <- match(rownames(Ometh_tmp), rownames(ref_meth))

	  #the boolean vector indicating which one is not NA
	  #In this way, Ometh_comm and ref_meth_comm have the same CpG site name for each row.
	  no_na <- !is.na(ind)
	  Ometh_comm <- Ometh_tmp[no_na,]
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
	  rownames(P_matr) <- colnames(ref_meth_comm)
	  colnames(P_matr) <- colnames(Ometh_comm)
	}else{
	  P_matr <- cell_prop
	}

	#filter out rare cell types
	del_celltype <- which(apply(P_matr, 1, quantile, 0.75) < thd)

	if(length(del_celltype) > 0){
	  P_matr <- P_matr[-del_celltype, ]
	  P_matr <- P_matr / colSums(P_matr)
	  K <- K - length(del_celltype)
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

	colnames(beta_pval) <- rownames(P_matr)
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

	colnames(gamma_tilde_pval) <- rownames(P_matr)
	rownames(gamma_tilde_pval) <- rownames(Ometh)

	########################################################
	# step 5. Joint significance method followed by squaring
	########################################################
	#pval_mediator
	pval_mediator <- beta_pval
	colnames(pval_mediator) <- rownames(P_matr)
	rownames(pval_mediator) <- rownames(Ometh)

	#pval_outcome
	pval_outcome <- gamma_tilde_pval
	colnames(pval_outcome) <- rownames(P_matr)
	rownames(pval_outcome) <- rownames(Ometh)


	if(maxp_sq == TRUE){
		#pval_joint_sq
		pval_joint <- matrix(NA, m, K)
		for(j in 1:m){
			for(k in 1:K){
				pval_joint[j,k] <- max(beta_pval[j,k], gamma_tilde_pval[j,k])
			}
		}
		pval_joint_sq <- pval_joint^2
		colnames(pval_joint_sq) <- rownames(P_matr)
		rownames(pval_joint_sq) <- rownames(Ometh)
	}else{
		pval_joint_sq <- NA
	}

	if(MCP.type == "FDR"){
		pval_joint_MultiMed <- NULL
		for(k in 1:K){
			tmp <- medTest.SBMH(as.vector(pval_mediator[ ,k]), as.vector(pval_outcome[ ,k]), MCP.type = "FDR")
			tmp <- as.vector(tmp)
			pval_joint_MultiMed <- cbind(pval_joint_MultiMed, tmp)
		}
		colnames(pval_joint_MultiMed) <- rownames(P_matr)
		rownames(pval_joint_MultiMed) <- rownames(Ometh)
	}else if(MCP.type == "FWER"){
		pval_joint_MultiMed <- NULL
		for(k in 1:K){
			tmp <- medTest.SBMH(as.vector(pval_mediator[ ,k]), as.vector(pval_outcome[ ,k]), MCP.type = "FWER")
			tmp <- as.vector(tmp)
			pval_joint_MultiMed <- cbind(pval_joint_MultiMed, tmp)
		}
		colnames(pval_joint_MultiMed) <- rownames(P_matr)
		rownames(pval_joint_MultiMed) <- rownames(Ometh)
	}else{
	    stop("MCP.type must be specified as either \"FDR\" or \"FWER\"!\n")
	}
	

	#return values
	out <- list()
	out$pval_mediator <- pval_mediator
	out$pval_outcome <- pval_outcome
	out$pval_joint_sq <- pval_joint_sq
	out$pval_joint_MultiMed <- pval_joint_MultiMed
	
	out$P_matr <- P_matr
	return(out)
}
