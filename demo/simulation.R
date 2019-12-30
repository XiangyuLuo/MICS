rm(list = ls())
set.seed(12022019)
setwd("E:/EWAS_mediation/simulation/simulation_v3/")

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

####################################################
#Generate the exposure, confounders, and the outcome
S <- rnorm(n) #the exposure
D <- 2        #the confounder number
X <- matrix(rnorm(D*n), n, D) #confounders

Y <- 0.5*S + 0.5*X[ ,1] + 0.5*X[ ,2] + rnorm(n, sd = 0.5)


################################################
#Generate the coefficients in Part I
#the effects of the exposure on the mediators
beta_matr <- array(0, dim=c(m,K))

m_common <- 3
max_signal <- 0.025
min_signal <- 0.02

signs <- sample(c(-1,1), m_common*K, replace=T)
beta_matr[1:m_common,1:K] <- signs * runif(m_common*K, min=min_signal, max=max_signal)

m_seperate <- 3
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


##################################################
#Generate the effect of mediators on the outcome
gamma_tilde_matr <- matrix(0, m, K)
max_signal_gamma <- 0.032
min_signal_gamma <- 0.027
for(k in 1:K){
	signs <- sample(c(-1,1), m_seperate, replace=T)
	gamma_tilde_matr[(k-1)*m_seperate+1:m_seperate, k] <- signs * runif(m_seperate, min=min_signal_gamma, 
													max=max_signal_gamma)
	#if(k == 3){
	#	signs <- sample(c(-1,1), m_seperate, replace=T)
	#	gamma_tilde_matr[(k-3)*m_seperate+1:m_seperate, k] <- signs * runif(m_seperate, min=min_signal_gamma,
	#												max=max_signal_gamma)		
	#}
} 



####################################################
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

#sum(Ometh > 1)
Ometh[Ometh > 1] <- 1

#sum(Ometh < 0)
Ometh[Ometh < 0] <- 0

rownames(Ometh) <- paste0("cpg", 1:m)
rownames(mu_matr) <- paste0("cpg", 1:m)

save(list=c("S", "Ometh", "Y", "X", "mu_matr"), file = "simulated_data.RData")
save(list=c("u", "P_matr", "beta_matr", "gamma_tilde_matr"), file = "simulation_truth.RData")


############################################################################
# Conduct MICS 
############################################################################

rm(list = ls())
load("simulated_data.RData")

library(qqman)
library(gplots)
library(MICS)

out <- mics(meth_data = Ometh, S = S, X = X, Y = Y, meth_ref = mu_matr)
pval_joint_sq <- out$pval_joint_sq
K <- 5
m <- 10000

pdf("pval_joint_significance_sq.pdf")
for(k in 1:K){
  qq(pval_joint_sq[,k], main = paste0("Cell type ", k))
}
dev.off()


#########################
#control FDR via q-values
#########################
library(qvalue)
qvalue_matr <- matrix(NA, nrow(pval_joint_sq), ncol(pval_joint_sq))

for(k in 1:K){
	qvalue_matr[ ,k] <- qvalue(pval_joint_sq[ ,k])$qvalues
}


for(k in 1:K){
	print(sum(qvalue_matr[,k] < 0.1))
}

for(k in 1:K){
	print(which(qvalue_matr[,k] < 0.1))
}


#TPR and FPR
ind_tmp <- which(qvalue_matr[,1] < 0.1)
sum(ind_tmp <= 3) / 3
sum(ind_tmp > 3) / (m-3)

ind_tmp <- which(qvalue_matr[,2] < 0.1)
sum((ind_tmp > 3)&(ind_tmp <= 6)) / 3
sum((ind_tmp <= 3)|(ind_tmp > 6)) / (m-3)

ind_tmp <- which(qvalue_matr[,3] < 0.1)
sum((ind_tmp > 6)&(ind_tmp <= 9)) / 3
sum((ind_tmp <= 6)&(ind_tmp > 9)) / (m-3)

ind_tmp <- which(qvalue_matr[,4] < 0.1)
sum((ind_tmp > 9)&(ind_tmp <= 12)) / 3
sum((ind_tmp <= 9)|(ind_tmp > 12)) / (m-3)

ind_tmp <- which(qvalue_matr[,5] < 0.1)
sum((ind_tmp > 12)&(ind_tmp <= 15)) / 3
sum((ind_tmp <= 12)|(ind_tmp > 15)) / (m-3)

#aggregated level
ind <- which(rowSums(qvalue_matr < 0.1) >0 )
sum(ind <= 15) / 15
sum(ind > 15) / (m-15)

#########################
#heatmaps
#########################
load("simulation_truth.RData")
colors_func <- colorRampPalette(c("grey", "black"))
colors <- colors_func(5)

pdf("beta_matr.pdf")
pattern_truth <- abs(beta_matr[1:50, ])
pattern_truth[pattern_truth>0] <- 1
main_title <- "   Effects of the exposure \n on the methylation"
heatmap.2(pattern_truth, col = c("grey", "black"), scale = "none", 
            key = TRUE, key.xlab = "\n  0=no effect\n1= with effect", Colv = FALSE, 
            Rowv = FALSE, density.info = "none", trace = "none", 
            dendrogram = "none", ylab = "", 
            xlab = "", margins = c(7, 1), main = "", 
            labCol =  1:K, labRow = FALSE, 
            cexRow = 5, srtCol = 0, cexCol = 2.5, adjCol = c(NA,1), keysize=1.5, 
		key.par=list(cex.main=1.8, cex.lab=1.8, cex.axis=1.8),
		lmat = rbind(c(4,3),c(2,1)), lhei=c(1.35,5), lwid=c(1.35,5), sepcolor="white",
            colsep=1:K)
mtext("Cell types", side=1, line=3, cex=2.5)
mtext(paste0("First ", 50, " CpG sites"), side=2, line=0.1, cex=2.5)
title(main_title, line= -1.2, cex.main = 2.5)
dev.off()



beta_pval <- out$pval_mediator

pdf("beta_matr_pval.pdf")
pval_matr <- beta_pval[1:50, ]
pval_matr[pval_matr == 0] <- 10^(-100)
main_title <- "        The minus log10 p-values \n from MICS part I"
heatmap.2(-log10(pval_matr), col = colors, scale = "none", 
            key = TRUE, key.xlab = "-log10(p-value)", Colv = FALSE, 
            Rowv = FALSE, density.info = "none", trace = "none", 
            dendrogram = "none", ylab = "", 
            xlab = "", margins = c(7, 1), main = "", 
            labCol = 1:K, labRow = FALSE, 
            cexRow = 5, srtCol = 0, cexCol = 2.5, adjCol = c(NA,1), keysize=1.5, 
		key.par=list(cex.main=1.8, cex.lab=1.8, cex.axis=1.8),
		lmat = rbind(c(4,3),c(2,1)), lhei=c(1.35,5), lwid=c(1.35,5), sepcolor="white",
            colsep=1:K)
mtext("Cell types", side=1, line=3, cex=2.5)
mtext(paste0("First ", 50, " CpG sites"), side=2, line=0.1, cex=2.5)
title(main_title, line= -1.2, cex.main = 2.5)
dev.off()



pdf("gamma_tilde_matr.pdf")
pattern_truth <- abs(gamma_tilde_matr[1:50, ])
pattern_truth[pattern_truth>0] <- 1
main_title <- "       Effects of the methylation \n on the outcome"
heatmap.2(pattern_truth, col = c("grey", "black"), scale = "none", 
            key = TRUE, key.xlab = "\n  0=no effect\n1= with effect", Colv = FALSE, 
            Rowv = FALSE, density.info = "none", trace = "none", 
            dendrogram = "none", ylab = "", 
            xlab = "", margins = c(7, 1), main = "", 
            labCol =  1:K, labRow = FALSE, 
            cexRow = 5, srtCol = 0, cexCol = 2.5, adjCol = c(NA,1), keysize=1.5, 
		key.par=list(cex.main=1.8, cex.lab=1.8, cex.axis=1.8),
		lmat = rbind(c(4,3),c(2,1)), lhei=c(1.35,5), lwid=c(1.35,5), sepcolor="white",
            colsep=1:K)
mtext("Cell types", side=1, line=3, cex=2.5)
mtext(paste0("First ", 50, " CpG sites"), side=2, line=0.1, cex=2.5)
title(main_title, line= -1.2, cex.main = 2.5)
dev.off()


colors <- colors_func(5)

gamma_tilde_pval <- out$pval_outcome

pdf("gamma_tilde_matr_pval.pdf")
pval_matr <- gamma_tilde_pval[1:50, ]
pval_matr[pval_matr == 0] <- 10^(-100)
main_title <- "        The minus log10 p-values \n from MICS part II"
heatmap.2(-log10(pval_matr), col = colors, scale = "none", 
            key = TRUE, key.xlab = "-log10(p-value)", Colv = FALSE, 
            Rowv = FALSE, density.info = "none", trace = "none", 
            dendrogram = "none", ylab = "", 
            xlab = "", margins = c(7, 1), main = "", 
            labCol = 1:K, labRow = FALSE, 
            cexRow = 5, srtCol = 0, cexCol = 2.5, adjCol = c(NA,1), keysize=1.5, 
		key.par=list(cex.main=1.8, cex.lab=1.8, cex.axis=1.8),
		lmat = rbind(c(4,3),c(2,1)), lhei=c(1.35,5), lwid=c(1.35,5), sepcolor="white",
            colsep=1:K)
mtext("Cell types", side=1, line=3, cex=2.5)
mtext(paste0("First ", 50, " CpG sites"), side=2, line=0.1, cex=2.5)
title(main_title, line= -1.2, cex.main = 2.5)
dev.off()



