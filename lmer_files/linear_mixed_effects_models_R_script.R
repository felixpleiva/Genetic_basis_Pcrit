# Script to relate Pcrit to HIF sequence data # 

# Author: Wilco Verberk
# Date: 20230611
# Modifications: Felix P Leiva (20230628) and Courtney H Babin (20230710)

#clear the work environment
rm(list=ls()) 

getwd()#to check

#Libraries
library (lme4)
library(lmerTest)
library(visreg)
library(car)
library(bbmle)
library(r2glmm)
# ------------------------------------------------------------------------------
#load new species with Pcrit data
d <- read.csv("HIF2A_long_forms_only_LD_Pcrit_data_matrix.csv")

# have a look at the data
head(d)
d$DAPC_group<-as.factor(d$DAPC_group) #convert DAPC_group to factor
d$clade<-as.factor(d$clade) #convert clade to factor
d$TSD<-as.factor(d$TSD) #convert TSD to factor (coded gar as "c" for placeholder purposes only)
table(d$species)

# convert to a single Pcrit data column, and add temperature as a separate column, yielding one row for each of the temperature
d15<-d[,!names(d) %in% c("Pc24", "Pc28")]
d15$temp<-15
names(d15)[which(names(d15)=="Pc15")]<-"Pcrit"

d24<-d[,!names(d) %in% c("Pc15", "Pc28")]
d24$temp<-24
names(d24)[which(names(d24)=="Pc24")]<-"Pcrit"

d28<-d[,!names(d) %in% c("Pc24", "Pc15")]
d28$temp<-28
names(d28)[which(names(d28)=="Pc28")]<-"Pcrit"

d<-rbind(d15,d24,d28) # this should result in 3 x the number of obs.

###########################################################################
# Analyzing the data, compare models and apply the correction for degrees of
# freedom for the model with the highest support
############################################################################

m0 <- lmer(Pcrit ~ (1|gene_species_ID), data = d) #null model
summary(m0)

m1 <- lmer(Pcrit ~ temp + (1|gene_species_ID), data = d)
summary(m1)

m2 <- lmer(Pcrit~LD1+(1|gene_species_ID),data=d)
summary(m2)

m3 <- lmer(Pcrit~LD1+temp+(1|gene_species_ID),data=d)
summary(m3)

m4 <- lmer(Pcrit~LD2+(1|gene_species_ID),data=d)
summary(m4)

m5 <- lmer(Pcrit~LD2+temp+(1|gene_species_ID),data=d)
summary(m5)

m6 <- lmer(Pcrit~LD3+(1|gene_species_ID),data=d)
summary(m6)

m7 <- lmer(Pcrit~LD3+temp+(1|gene_species_ID),data=d)
summary(m7)

m8 <- lmer(Pcrit~LD1+LD2+(1|gene_species_ID),data=d)
summary(m8)

m9 <- lmer(Pcrit~LD1+LD2+temp+(1|gene_species_ID),data=d)
summary(m9)

m10 <- lmer(Pcrit~LD1+LD3+(1|gene_species_ID),data=d)
summary(m10)

m11 <- lmer(Pcrit~LD1+LD3+temp+(1|gene_species_ID),data=d)
summary(m11)

m12 <- lmer(Pcrit~LD2+LD3+(1|gene_species_ID),data=d)
summary(m12)

m13 <- lmer(Pcrit~LD2+LD3+temp+(1|gene_species_ID),data=d)
summary(m13)

m14 <- lmer(Pcrit~LD1+LD2+LD3+(1|gene_species_ID),data=d)
summary(m14)

m15 <- lmer(Pcrit~LD1+LD2+LD3+temp+(1|gene_species_ID),data=d)
summary(m15)

#-------------------------------------------------------------------------------
#Table 1
#-------------------------------------------------------------------------------
fit.list.table.1 <- list(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15)

fit.names.table.1 <-c("Null",
                      "temp + (1|gene_species_ID)",
                      "LD1+(1|gene_species_ID)",
                      "LD1+temp+(1|gene_species_ID)",
                      "LD2+(1|gene_species_ID)",
                      "LD2+temp+(1|gene_species_ID)",
                      "LD3+(1|gene_species_ID)",
                      "LD3+temp+(1|gene_species_ID)",
                      "LD1+LD2+(1|gene_species_ID)",
                      "LD1+LD2+temp+(1|gene_species_ID)",
                      "LD1+LD3+(1|gene_species_ID)",
                      "LD1+LD3+temp+(1|gene_species_ID)",
                      "LD2+LD3+(1|gene_species_ID",
                      "LD2+LD3+temp+(1|gene_species_ID)",
                      "LD1+LD2+LD3+(1|gene_species_ID)",
                      "LD1+LD2+LD3+temp+(1|gene_species_ID)"
                      
)

#compare by using AICc
HIF2A.fit.table <- AICctab(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,mnames = fit.names.table.1, base = TRUE, weights = TRUE, logLik = TRUE)
HIF2A.fit.table

# Export tables to .csv
write.csv(HIF2A.fit.table,"HIF2A_Table_1_Model_comparison.csv",row.names = TRUE)

#save results of lmer summary for best model
sink("HIF2A_m9_lmer_results.txt")
summary(m9)
sink()

#obtain R-squared for best model
r2beta(m9, method = "nsj")

#save results for r2beta
sink("HIF2A_m9_r2_results.txt")
r2beta(m9, method = "nsj")
sink()

#correct p-values for each effect in best model (repeat for each effect)
t<- 7.604 #value from m9 best lmer model results for (Intercept)
degf<-91.95555 #value from df column of best model for (Intercept)
2*pt(q=t,df=degf,lower.tail=FALSE) # p-value formula to compare to model p-value for (Intercept)
2*pt(q=t,df=24,lower.tail=FALSE) # p-value formula to correct for 28 species-1-3 fixed effects for (Intercept)

#obtain 95% confidence intervals for parameters
confint(m9)

#save results of confint(m9)
sink("HIF2A_m9_conf_int.txt")
confint(m9)
sink()

################################################################################
# Figure
################################################################################
#border around plot
par(mai=c(1,1,0.1,0.1))

#plotn(REPEAT FOR ALL SIGNIFICANT LDs)
xlab<- expression(italic("HIF2A")*" linear discriminant 1")

#use pch = c(1,3) for models that only have two LDs
visreg(m9, "LD1", points = list(cex = 1.8, pch = c(1, 2, 3)[d$TSD], col = c("black", "burlywood4", "blue", "chartreuse4", "orange")[d$clade], lwd=2),
       xlab=xlab, ylab="Critical oxygen tension (kPa)", 
       cex.axis=1.3,cex.lab=1.3, top="points")

# saving session information with all packages versions for reproducibility purposes
sink("Data_analyses_Table_1_and_Figure_1_R_session.txt")
sessionInfo()
sink()

# ------------------------------------------------------------------------------
################################################################################
############################ END OF SCRIPT #####################################
################################################################################