# PGLS script using HIF2A as example # 

#Libraries
library (ape)
library(geiger)
library(caper)
library(tidyverse)
library(phytools)
# ------------------------------------------------------------------------------
#load new species with Pcrit data (duplicated "gene_species_ID" and add new column labeled "tiplabel" before importing)
d <- read.csv("HIF2A_data.csv", header = TRUE, row.names = 1)

# have a look at the data
head(d)
d$DAPC_group<-as.factor(d$DAPC_group) #convert DAPC_group to factor
d$clade<-as.factor(d$clade) #convert clade to factor
d$TSD<-as.factor(d$TSD) #convert TSD to factor (coded gar as "c" for placeholder purposes only)
table(d$species)

#Reading in the tree
HIF2Atree <- read.newick("HIF2A_tree.newick")

check <- name.check(phy = HIF2Atree, data = d, 
                    data.names = d$tiplabel)
check

# correct differences in spelling
rownames(d)[-which(rownames(d)%in%HIF2Atree$tip.label)]

HIF2Atree$tip.label
HIF2Atree$tip.label[15]<-"hif2Aa_s1_Oncorhynchus_tshawytscha_4583"
HIF2Atree$tip.label[18]<-"hif2Aa_s2_Oncorhynchus_tshawytscha_0055"
HIF2Atree$tip.label[34]<-"hif2Aa_Cyprinidon_variegatus_0627"

tre <- drop.tip(HIF2Atree, HIF2Atree$tip.label[-which(HIF2Atree$tip.label%in%rownames(d))])
tre#42 species

#convert to a dataframe

d <- as.data.frame(d)
class(d)

###########################################################################
# Analyzing the data
############################################################################

#may need to remove node labels with tre$node.label<-NULL
HIF2A <- comparative.data(phy = tre, data = d, 
                         names.col = tiplabel, vcv = TRUE, 
                         na.omit = FALSE, warn.dropped = TRUE)

#Run full model then reduce (retain variables using threshold p <= 0.1)
#Repeat for Pc24 and Pc28, and lambda = 0.001)
model.pgls <- pgls(Pc15 ~ LD1+LD2+LD3, 
                   data = HIF2A, lambda = "ML")

anova(model.pgls)

#save results of (reduced) ANOVA
sink("HIF2A_PGLS_Pc15_ANOVA_results.txt")
anova(model.pgls)
sink()

#save results of (reduced) PGLS model summary
sink("HIF2A_PGLS_Pc15_summary_results.txt")
summary(model.pgls)
sink()

#get model diagnostic plots for (reduced) PGLS
par(mfrow = c(2, 2))
plot(model.pgls)

# Create a likelihood profile of the lambda estimate (from reduced PGLS model)
lambda.profile <- pgls.profile(model.pgls, "lambda")
# Plot the likelihood profile
plot(lambda.profile)

# Extract the confidence intervals on lambda
pgls.confint(model.pgls, "lambda")$ci.val

#save results of pgls.confint
sink("HIF2A_PGLS_Pc15_lambda_conf_int.txt")
pgls.confint(model.pgls, "lambda")$ci.val
sink()

#plot results - change x = for each LD, y = for each Pcrit temperature, labs x = expression for italic label,
#slope coefficients to match bracket number, and adjust number of shapes to represent TSD for each gene
ggplot(d, aes(x = LD1, y = Pc15, color = clade, shape = TSD)) + 
  geom_point(size = 5) + ylim(0,15) + labs(x = expression(italic("HIF2A")*" linear discriminant 1"), y = "Critical oxygen tension (kPa)", color = "Clade", shape = "Paralog") + 
  geom_abline(slope = coefficients(model.pgls)[2], intercept = coefficients(model.pgls)[1]) + 
  theme_bw()+ theme(text = element_text(size = rel(4)), axis.text = element_text(size = 15), legend.position = "none") + 
  scale_shape_manual(values = c(19, 17, 3)) + 
  scale_color_manual(values = c("#000000", "#8B7355", "#0000FF", "#458B00", "#FFA500"))

# saving session information with all packages versions for reproducibility purposes
sink("Data_analyses_Figure_4_R_session.txt")
sessionInfo()
sink()

# ------------------------------------------------------------------------------
################################################################################
############################ END OF SCRIPT #####################################
################################################################################