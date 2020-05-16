MICS is short for "Mediation In a Cell-type-Specific fashion" and is implemented as an R package here to test cell-type-specific mediation effects in the epigenome wide mediation analysis.

The R code to run MICS can be seen in the example of the function mics.

Installation instruction:
```
devtools::install_github("XiangyuLuo/MICS")
```

Example code:
```
library(MICS)

#load example data
data(example_data)

#perform mics
out <- mics(meth_data = Ometh, S = S, X = X, Y = Y, cell_prop = P_matr, 
                       MCP.type = "FDR", maxp_sq = TRUE)

pval_MultiMed <- out$pval_joint_MultiMed

#FDR threshold
fdr_thred <- 0.2

ind1 <- which(pval_MultiMed[,1] < fdr_thred)

ind2 <- which(pval_MultiMed[,2] < fdr_thred)

ind3 <- which(pval_MultiMed[,3] < fdr_thred)

#detected CpG sites in cell type 1
ind1

#detected CpG sites in cell type 2
ind2

#detected CpG sites in cell type 3
ind3
```

Check how to use the mics function using:

```
?mics
```


Once you have any questions, please email Xiangyu Luo at xiangyuluo@ruc.edu.cn or Zhonghua Liu at zhhliu@hku.hk.
