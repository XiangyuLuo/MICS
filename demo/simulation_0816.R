rm(list = ls())
#install.packages("MICS_0.1.2.tar.gz", repos = NULL)

#############################################################################
#Simulate data
#############################################################################
set.seed(03292019)

#Generate samples from the Dirichlet distribution
rDirichlet <- function(alpha_vec){
	num <- length(alpha_vec)
	temp <- rgamma(num, shape = alpha_vec, rate = 1)
	return(temp / sum(temp))
}

n <- 400     #number of samples 
m <- 10000   #number of CpG sites
K <- 3       #underlying cell type number
	
#Generate the baseline methylation profiles mu
#cell type 1
methy1 <- rbeta(m,3,3)

#cell type 2 
methy2 <- methy1 + rnorm(m, sd=0.01)
ind <- sample(1:m, m/10) 
methy2[ind] <- rbeta(length(ind),3,3)

#cell type 3
methy3 <- methy1 + rnorm(m, sd=0.01)
ind <- sample(1:m, m/10) 
methy3[ind] <- rbeta(length(ind),3,3)

mu_matr <- cbind(methy1, methy2, methy3)

#Generate the exposure, confounders
S <- rnorm(n) #the exposure
D <- 2        #the confounder number
X <- matrix(rnorm(D*n), n, D) #confounders


###Generate the cell-type-specific methylation values
Marray <- array(NA, dim = c(m, n, K))
ind_mediator <- matrix(c( 1, 5,
				  4, 8,
				  6, 10), K, 2, byrow = TRUE)
for(k in 1:K){
	for(i in 1:n){
		Marray[ ,i,k] <- mu_matr[ ,k] + 0.01*X[i,1] + 0.01*X[i,2] + rnorm(m, sd = 0.02)
		Marray[ind_mediator[k,1]:ind_mediator[k,2],i,k] <- 
			Marray[ind_mediator[k,1]:ind_mediator[k,2],i,k] + 0.03 * S[i]
	}
}

###Generate the outcome variable
ind_outcome <- matrix(c(  4, 5,
				  7, 8,
				  9, 10), K, 2, byrow = TRUE)

Y <- 0.5*S + 0.3*X[ ,1] + 0.3*X[ ,2] + 0.5 * colSums(Marray[ind_outcome[1,1]:ind_outcome[1,2], ,1]) +
						0.5 * colSums(Marray[ind_outcome[2,1]:ind_outcome[2,2], ,2]) +
						0.6 * colSums(Marray[ind_outcome[3,1]:ind_outcome[3,2], ,3]) 


#Generate the cellular proportions
P_matr <- sapply(1:n, function(i){
				rDirichlet(rep(2, K))
			})


#Generate the observed methylation profiles 
Ometh <- NULL
for(i in 1:n){
	tmp <- Marray[ ,i, ] %*% P_matr[ ,i]
	Ometh <- cbind(Ometh, tmp)
}

#sum(Ometh > 1)
#sum(Ometh < 0)
Ometh[Ometh > 1] <- 1
Ometh[Ometh < 0] <- 0

rownames(Ometh) <- paste0("cpg", 1:m)
rownames(mu_matr) <- paste0("cpg", 1:m)
colnames(mu_matr) <- paste0("cell_type", 1:K)
rownames(P_matr) <- paste0("cell_type", 1:K)


#############################################################################
#Conduct MICS
#############################################################################
library(MICS)

#out <- mics(meth_data = Ometh, S=S, X=X, Y=Y, cell_prop = P_matr)
out <- mics(meth_data = Ometh, S=S, X=X, Y=Y, meth_ref = mu_matr)

library(qqman)
png(file = "pval_mediator.png", width=5000,height=2000,res=600)
par(mfrow = c(1,3))
for(k in 1:K){
  qq(out$pval_mediator[,k], main = paste0("Cell type ", k))
}
dev.off()

png(file = "pval_outcome.png", width=5000,height=2000,res=600)
par(mfrow = c(1,3))
for(k in 1:K){
  qq(out$pval_outcome[,k], main = paste0("Cell type ", k))
}
dev.off()

png(file = "pval_joint_sq.png", width=5000,height=2000,res=600)
par(mfrow = c(1,3))
for(k in 1:K){
  qq(out$pval_joint_sq[,k], main = paste0("Cell type ", k))
}
dev.off()
