rm(list = ls())
set.seed(03292019)

#Generate samples from the Dirichlet distribution
rDirichlet <- function(alpha_vec){
	num <- length(alpha_vec)
	temp <- rgamma(num, shape = alpha_vec, rate = 1)
	return(temp / sum(temp))
}

n <- 400     #number of samples 
m <- 10000   #number of CpG sites
K <- 5       #underlying cell type number
	
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

#cell type 4
methy4 <- methy1 + rnorm(m, sd=0.01)
ind <- sample(1:m, m/10) 
methy4[ind] <- rbeta(length(ind),3,3)

#cell type 5
methy5 <- methy1 + rnorm(m, sd=0.01)
ind <- sample(1:m, m/10) 
methy5[ind] <- rbeta(length(ind),3,3)

mu_matr <- cbind(methy1, methy2, methy3, methy4, methy5)

#Generate the exposure, confounders, and the outcome
S <- rnorm(n) #the exposure
D <- 2        #the confounder number
X <- matrix(rnorm(D*n), n, D) #confounders

Y <- 0.5*S + 0.5*X[ ,1] + 0.5*X[ ,2] + rnorm(n, sd = 0.5)

#Generate the coefficients in Part I
#the effects of the exposure on the mediators
beta_matr <- array(0, dim=c(m,K))

m_common <- 10
max_signal <- 0.15
min_signal <- 0.07

signs <- sample(c(-1,1), m_common*K, replace=T)
beta_matr[1:m_common,1:K] <- signs * runif(m_common*K, min=min_signal, max=max_signal)

m_seperate <- 10
signs <- sample(c(-1,1), m_seperate*2, replace=T)
beta_matr[m_common+(1:m_seperate),1:2] <- signs * runif(m_seperate*2, min=min_signal, max=max_signal)

for(k in 3:K){
	signs <- sample(c(-1,1), m_seperate, replace=T)
	beta_matr[m_common+(k-2)*m_seperate+(1:m_seperate),k] <- signs * runif(m_seperate, min=min_signal, max=max_signal)
}

#the effects of the confounders on the mediators
h_vec1 <- rnorm(m, sd=0.01)
h_vec2 <- rnorm(m, sd=0.01)

#generate the cellular compositions 
P_matr <- sapply(1:n, function(i){
				rDirichlet(rep(2, K))
			})

#Generate the effect of mediators on the outcome
gamma_tilde_matr <- matrix(0, m, K)
for(k in 1:K){
	signs <- sample(c(-1,1), m_seperate, replace=T)
	gamma_tilde_matr[(k-1)*m_seperate+1:m_seperate, k] <- signs * runif(m_seperate, min=0.03, max=0.1)
	if(k == 3){
		signs <- sample(c(-1,1), m_seperate, replace=T)
		gamma_tilde_matr[(k-3)*m_seperate+1:m_seperate, k] <- signs * runif(m_seperate, min=0.03, max=0.1)		
	}
} 

#Generate u conditional on Y
u <- array(NA, dim=c(n, m, K))
for(i in 1:n){
	epsilon_matr <- matrix(rnorm(m*K, sd = 0.015), m, K)
	#epsilon_matr <- matrix(c(rnorm(m*2, sd = 0.02), rnorm(m*2, sd = 0.01), rnorm(m, sd = 0.015)), m, K)
	u[i,,] <- mu_matr + beta_matr*S[i] + h_vec1*X[i,1] + h_vec2*X[i,2] + gamma_tilde_matr * Y[i] + epsilon_matr  
 
}


#Generate the observed methylation profiles 
Ometh <- NULL
for(i in 1:n){
	tmp <- u[i,,] %*% P_matr[ ,i]
	Ometh <- cbind(Ometh, tmp)
}

Ometh[Ometh > 1] <- 1

Ometh[Ometh < 0] <- 0

rownames(Ometh) <- paste0("cpg", 1:m)
rownames(mu_matr) <- paste0("cpg", 1:m)
colnames(mu_matr) <- paste0("cell_type", 1:K)

#############################################################################
#Conduct MICS
#############################################################################
library(MICS)
pval_joint_sq <- mics(Ometh, mu_matr, S, X, Y)

library(qqman)
pdf("pval_joint_significance_sq.pdf")
for(k in 1:K){
  qq(pval_joint_sq[,k], main = paste0("Cell type ", k))
}
dev.off()

