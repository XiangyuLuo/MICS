rm(list = ls())
set.seed(03092020)

#set the working directory
############################################################################
# Generate Data 
############################################################################
#Generate samples from the Dirichlet distribution
rDirichlet <- function(alpha_vec){
  num <- length(alpha_vec)
  temp <- rgamma(num, shape = alpha_vec, rate = 1)
  return(temp / sum(temp))
}

n <- 400     #number of samples 
m <- 10000   #number of CpG sites
K <- 5       #underlying cell type number

####################################################
#Generate the exposure, confounders, and the outcome
S <- rnorm(n) #the exposure
D <- 2        #the confounder number
X <- matrix(rnorm(D*n), n, D) #confounders

##############################################	
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
colnames(mu_matr) <- c("celltype 1", "celltype 2", "celltype 3",
                       "celltype 4", "celltype 5")

################################################
#Generate the coefficients
#the effects of the exposure on the mediators
beta_matr <- array(0, dim=c(m,K))

#the cpg sites affected by exposure for each cell type
cpg_type1 <- 1:75
cpg_type2 <- 26:100
cpg_type3 <- 51:125
cpg_type4 <- 76:150
cpg_type5 <- 101:175

#coefficients
max_signal <- 0.25
min_signal <- 0.1

signs <- sample(c(-1,1), 75, replace=T)
beta_matr[cpg_type1, 1] <- signs * runif(75, min=min_signal, max=max_signal)

signs <- sample(c(-1,1), 75, replace=T)
beta_matr[cpg_type2, 2] <- signs * runif(75, min=min_signal, max=max_signal)

signs <- sample(c(-1,1), 75, replace=T)
beta_matr[cpg_type3, 3] <- signs * runif(75, min=min_signal, max=max_signal)

signs <- sample(c(-1,1), 75, replace=T)
beta_matr[cpg_type4, 4] <- signs * runif(75, min=min_signal, max=max_signal)

signs <- sample(c(-1,1), 75, replace=T)
beta_matr[cpg_type5, 5] <- signs * runif(75, min=min_signal, max=max_signal)

#the effects of the confounders on the mediators
h_vec1 <- rnorm(m, sd=0.01)
h_vec2 <- rnorm(m, sd=0.01)




##################################################
#Generate the effect of mediators on the outcome
gamma_matr <- matrix(0, m, K)
max_signal_gamma <- 0.035
min_signal_gamma <- 0.025

#the cpg sites affected by mediators for each cell type
cpg_type1_m <- 1:50
cpg_type2_m <- 26:75
cpg_type3_m <- 51:100
cpg_type4_m <- 76:125
cpg_type5_m <- 101:150

signs <- sample(c(-1,1), 50, replace=T)
gamma_matr[cpg_type1_m, 1] <- signs * runif(50, min=min_signal_gamma, max=max_signal_gamma)

signs <- sample(c(-1,1), 50, replace=T)
gamma_matr[cpg_type2_m, 2] <- signs * runif(50, min=min_signal_gamma, max=max_signal_gamma)

signs <- sample(c(-1,1), 50, replace=T)
gamma_matr[cpg_type3_m, 3] <- signs * runif(50, min=min_signal_gamma, max=max_signal_gamma)

signs <- sample(c(-1,1), 50, replace=T)
gamma_matr[cpg_type4_m, 4] <- signs * runif(50, min=min_signal_gamma, max=max_signal_gamma)

signs <- sample(c(-1,1), 50, replace=T)
gamma_matr[cpg_type5_m, 5] <- signs * runif(50, min=min_signal_gamma, max=max_signal_gamma)


#######################################################
#Generate the mediators and the outcome 

#the effects of the confounders on the outcome
eta1 <- 0.5
eta2 <- 0.5

#the effects of the exposure on the outcome
eta0 <- 0.5

#jointly generating u and Y is equivalent to first generating
# Y and then generating u conditional on Y

Y <- 0.5*S + 0.5*X[ ,1] + 0.5*X[ ,2] + rnorm(n, sd = 0.5)


#Generate u conditional on Y
u <- array(NA, dim=c(n, m, K))
for(i in 1:n){
  epsilon_matr <- matrix(rnorm(m*K, sd = 0.01), m, K)
  u[i,,] <- mu_matr + beta_matr*S[i] + h_vec1*X[i,1] + h_vec2*X[i,2] + 
            gamma_matr*Y[i] + epsilon_matr  
}

#############################################
#Generate the cellular compositions 
P_matr <- sapply(1:n, function(i){
  rDirichlet(c(4, 1, 3, 1, 1))
})

#############################################
#Generate the observed methylation profiles 
Ometh <- NULL
for(i in 1:n){
  tmp <- u[i,,] %*% P_matr[ ,i]
  Ometh <- cbind(Ometh, tmp)
}

#sum(Ometh > 1)
Ometh[Ometh > 1] <- 1

#sum(Ometh < 0)
Ometh[Ometh < 0] <- 0

rownames(Ometh) <- paste0("cpg", 1:m)
rownames(mu_matr) <- paste0("cpg", 1:m)

save(list=c("S", "Ometh", "Y", "X", "mu_matr"), file = "simulated_data.RData")
save(list=c("u", "P_matr", "beta_matr", "gamma_matr"), file = "simulation_truth.RData")


############################################################################
# Conduct MICS 
############################################################################

rm(list = ls())
load("simulated_data.RData")

#install MICS package first
library(MICS)

out <- mics(meth_data = Ometh, S = S, X = X, Y = Y, meth_ref = mu_matr, 
            MCP.type = "FDR", maxp_sq = TRUE)

out$pval_joint_sq[1:100, 1]


pval_MultiMed <- out$pval_joint_MultiMed
K <- 5
m <- 10000

fdr_thred <- 0.2

#the true cpg sites through the mediation pathway 
cpg_type1 <- 1:50
cpg_type2 <- 26:75
cpg_type3 <- 51:100
cpg_type4 <- 76:125
cpg_type5 <- 101:150

###cell type 1
#power 
sum(pval_MultiMed[cpg_type1, 1] < fdr_thred) / 50
#false positive rate
sum(pval_MultiMed[setdiff(1:m, cpg_type1), 1] < fdr_thred) / (m-50)

###cell type 2
#power 
sum(pval_MultiMed[cpg_type2, 2] < fdr_thred) / 50
#false positive rate
sum(pval_MultiMed[setdiff(1:m, cpg_type2), 2] < fdr_thred) / (m-50)

###cell type 3
#power 
sum(pval_MultiMed[cpg_type3, 3] < fdr_thred) / 50
#false positive rate
sum(pval_MultiMed[setdiff(1:m, cpg_type3), 3] < fdr_thred) / (m-50)

###cell type 4
#power 
sum(pval_MultiMed[cpg_type4, 4] < fdr_thred) / 50
#false positive rate
sum(pval_MultiMed[setdiff(1:m, cpg_type4), 4] < fdr_thred) / (m-50)

###cell type 5
#power 
sum(pval_MultiMed[cpg_type5, 5] < fdr_thred) / 50
#false positive rate
sum(pval_MultiMed[setdiff(1:m, cpg_type5), 5] < fdr_thred) / (m-50)



