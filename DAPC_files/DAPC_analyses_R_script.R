# Script for Discriminant Analysis of Principal Components # 

# Author: Courtney H Babin
# Date: 20230720

#Libraries
library(adegenet)

#Read in setup table and create object (example using HIF1A_setup)

HIF1A_setup <- read.table("HIF1A_setup.txt", header = F, row.names = 1)

#Replace NA values in data with the mean of the corresponding column

data <- HIF1A_setup                                              # Duplicate data frame
for(i in 1:ncol(HIF1A_setup)) {                                   # Replace NA in all columns
  data[ , i][is.na(data[ , i])] <- mean(data[ , i], na.rm = TRUE)
}

#clustering HIF1A_setup using min, retain all PCs in graph
grp1 <- find.clusters(data, criterion = "min", choose.n.clust=FALSE, n.iter=1000000, n.start=10000)
50
grp1

#perform DAPC on PCs and LDs, retain number of PCs that explain 90% of variance from graph, then retain total number of linear discriminants from DA eigenvalues graph
dapc1 <- dapc(data, grp1$grp)
19
3

dapc1

#create scatterplot with positioning for DA and PCA eigenvalues (change colors with myCol = c("blue", "orange", "green")
# Add clab = 0 if you want cluster labels removed
# Add xax=1, yax=3 if you want to plot axes (LDs) other than 1 and 2 against each other
#if you only want to represent ONE LD, use: scatter(dapc1,scree.da=FALSE)
scatter(dapc1, posi.da="topright", bg="white", scree.pca=TRUE, posi.pca="bottomright")

#make summary data from dapc1 a data frame object
dapcinfo=summary.data.frame(dapc1)

#write summary data object as csv
write.csv(dapcinfo, file = "dapcinfo.csv")

#make ind.coord LDs data frame as object
lds=data.frame(dapc1$ind.coord)

#write ind.coord as csv
write.csv(lds, file="lds.csv")

#make var.contr data frame as object
var_cont=data.frame(dapc1$var.contr)

#write var.contr as csv
write.csv(var_cont, file="var_cont.csv")

#loading plots for each axis (each LD) showing variable contribution with threshold of 0.90)
#change axis = number for whichever LD you want variable contributions 
contrib <- loadingplot(dapc1$var.contr, threshold=quantile(dapc1$var.contr,0.90), axis = 1)

contrib

#write var.values as csv
var_above_90=data.frame(contrib$var.values)
write.csv(var_above_90, file="var_above_90.csv")

# saving session information with all packages versions for reproducibility purposes
sink("Data_analyses_Supplemental_Figures_S1_through_S7_R_session.txt")
sessionInfo()
sink()

# ------------------------------------------------------------------------------
################################################################################
############################ END OF SCRIPT #####################################
################################################################################