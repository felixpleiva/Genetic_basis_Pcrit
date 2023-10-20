################################################################################
# Script to fit phylogenetic multilevel model using brms and extract corrected
# Pcrit for a given list of species. It is an modified version of the code published in 
#------------------------------------------------------------------------------
# Cleaning working space
rm(list = ls())
today <- format(Sys.Date(), "%Y%m%d")
# ------------------------------------------------------------------------------
# set the working directory to the folder containing this script:
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# ------------------------------------------------------------------------------
# check directory
getwd()
# ------------------------------------------------------------------------------
# libraries
library(brms)
library(ape)
library(dplyr)
library(phytools)
library(tidybayes)
library(bayestestR)
library(ggplot2)
library(corrplot)
# ------------------------------------------------------------------------------
set.seed(6955)# we need this to replicate the results
# ------------------------------------------------------------------------------
# General STAN specifications
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# ------------------------------------------------------------------------------
#load new species with new Pcrit data (3 add species)
data <- read.csv("Observed Pcrit data.csv")

#load the new tree. This tree object contains information on the
#relationship between species.

tree <- read.tree("Phylogenetic tree for 198 species of fish with Pcrit data.tre")
tree$tip.label#198 species

#-------------------------------------------------------------------------------
# apply complete.clases (exclude NAs) for the columns of interest
d<-data[complete.cases(data[,c("Estimated_C.Value",
                               "Temp",
                               "logMass",
                               "Pcrit..kPa.",
                               "Salt",
                               "animal",
                               "residMR")]),]

# explore whether our tree covers all the species we wanted it to include, and
# making sure that the species names in our database match those in the tree. We
# use the following code.
# 
setdiff(d$animal, as.character(tree$tip.label))#listed in our database but not in the tree
setdiff(as.character(tree$tip.label),d$animal)# listed in the tree but not in our database
# -------------------------------------------------------------------------------
# 
# There are 24 species not included in the data frame. These species does not
# have information on the variables selected in Line 60

# Lets now drop these species from the phylogenetic tree.
tree<-drop.tip(tree, setdiff(tree$tip.label, d$animal))
# 
# -------------------------------------------------------------------------------
# check again
setdiff(d$animal, as.character(tree$tip.label))#listed in our database but not in the tree
setdiff(as.character(tree$tip.label),d$animal)# listed in the tree but not in our database

# uffff, now we have the same species sin the tree and the dataframe

# ------------------------------------------------------------------------------
#count number of species to be included in the model
length(unique(d$animal))

# ------------------------------------------------------------------------------
# check if we have the new three species in the pruned tree

is_species_present <- function(tree, species_name) {
  tip_labels <- tree$tip.label
  return(species_name %in% tip_labels)
}

is_species_present(tree, "Esox_lucius")
is_species_present(tree, "Clupea_harengus")
is_species_present(tree, "Lepisosteus_oculatus")
# ------------------------------------------------------------------------------
# Using this information, we can construct a covariance matrix of species (Hadfield & Nakagawa, 2010).
# The power 0.4 was the most informative value used in Verberk et al 2022 (GCB)
B <- ape::vcv.phylo(compute.brlen(tree, power=0.4))  
# ------------------------------------------------------------------------------
# Plot the tree
pdf("Phylogenetic tree for 174 species of fish with Pcrit data.pdf",width =10,height = 10,useDingbats = FALSE)
plotTree(tree, fsize = 0.6,ftype = "i", type = "fan")
dev.off()

# save the tree for further analyses

write.tree(tree,"Phylogenetic tree for 174 species of fish with Pcrit data.tre")
write.tree(tree,"Phylogenetic tree for 174 species of fish with Pcrit data.nwk")
#-------------------------------------------------------------------------------
#use simplest names
names(d)
d$c_value_log <- log10(d$Estimated_C.Value)
d<-rename(d,temp_test=Temp)
d<-rename(d,pcrit_kpa=Pcrit..kPa.)
d<-rename(d,sal_test=Salt)
d<-rename(d,temp_acclim=Taccl)
d<-rename(d,res_metab_rate=residMR)
d<-rename(d,mass_log=logMass)

# and create a new column "species"
d$species <- d$animal
# ------------------------------------------------------------------------------
# 
# #Define priors
priors1=c(
  prior(normal(0,10), "b"),
  prior(normal(0,50), "Intercept"),
  prior(student_t(3,0,20), "sd"),
  prior(student_t(3,0,20), "sigma")
)
# 
# #-----------------------------------------------------------------------------
# fit and run the model using the same predictors as in Verberk et al 2022.
# Running time is about 4 hours!!!!!


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#  If you do NOT want to wait so long and load the .rds file below
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

best_model_final <- brm(
 pcrit_kpa~ c_value_log * temp_test  +
   mass_log * temp_test  +
   sal_test +
   res_metab_rate +
   (1|gr(animal, cov = B)) +
   (1|species),
 data = d,
 family = gaussian(),
 data2 = list(B = B),
 prior = priors1,
 sample_prior = TRUE,
 chains = 3, cores = 2,
 iter = 150, warmup = 75,
 control = list(adapt_delta = 0.999, max_treedepth = 20),
 save_pars = save_pars(all = TRUE))
 
 
# # save the model run
# saveRDS(best_model,"best_model.rds")

# Models are saved so that they don't need to be rerun for every session. 
# It may takes several hours

best_model_final <- readRDS("best_model.rds")

#-------------------------------------------------------------------------------
#extract the summary of the model
summary(best_model_final)

# extract fixed effect estimates of the model
fixedEstimates                       <-  brms::fixef(best_model_final, estimate = 'mean')
intercept                            <-  fixedEstimates['Intercept', 'Estimate']
c_value_log_Slope                    <-  fixedEstimates['c_value_log', 'Estimate']
temp_test_Slope                      <-  fixedEstimates['temp_test', 'Estimate']
mass_log_Slope                       <-  fixedEstimates['mass_log', 'Estimate']
sal_test_Slope                       <-  fixedEstimates['sal_test', 'Estimate']
res_metab_rate_Slope                 <-  fixedEstimates['res_metab_rate', 'Estimate']
temp_test_and_c_value_log_Slope      <-  fixedEstimates['c_value_log:temp_test', 'Estimate']
temp_test_and_mass_log_Slope         <-  fixedEstimates['temp_test:mass_log', 'Estimate']

# extract random effect estimates of the model
animal_Random  <-  brms::ranef(best_model_final, estimate = 'mean')$animal[, , 'Intercept'][as.character(d$animal), 'Estimate']
species_Random     <-  brms::ranef(best_model_final, estimate = 'mean')$species[, , 'Intercept'][d$species, 'Estimate']
#-------------------------------------------------------------------------------
# Calculate predicted response

d$pcrit_predicted   <-  intercept+(temp_test_Slope * (d$temp_test)) +
  (mass_log_Slope * (d$mass_log)) + 
  (c_value_log_Slope * (d$c_value_log)) + 
  (res_metab_rate_Slope * (d$res_metab_rate)) +
  (sal_test_Slope * mean(d$sal_test)) +
  (temp_test_and_mass_log_Slope * (d$temp_test)* (d$mass_log)) +
  (temp_test_and_c_value_log_Slope * (d$temp_test) * (d$c_value_log)) + 
  species_Random + 
  animal_Random
plot(pcrit_kpa~pcrit_predicted,data=d) # good correlation
abline(a=0,b=1)
plot(I(pcrit_kpa-pcrit_predicted)~pcrit_predicted,data=d) # plotting residuals vs predicted values

summary(lm(animal_Random~species_Random))# random effects are correlated


# check te test temp we have
summary(d$temp_test) #15,24,28


# 15 centigrade
d$pcrit_predicted15   <-  intercept+(temp_test_Slope * (15)) +
  (mass_log_Slope * (d$mass_log)) + 
  (c_value_log_Slope * (d$c_value_log)) + 
  (res_metab_rate_Slope * (d$res_metab_rate)) +
  (sal_test_Slope * mean(d$sal_test)) +
  (temp_test_and_mass_log_Slope * (15)* (d$mass_log)) +
  (temp_test_and_c_value_log_Slope * (15) * (d$c_value_log)) + 
  species_Random + 
  animal_Random

# 24 centigrade
d$pcrit_predicted24   <-  intercept+(temp_test_Slope * (24)) +
  (mass_log_Slope * (d$mass_log)) + 
  (c_value_log_Slope * (d$c_value_log)) + 
  (res_metab_rate_Slope * (d$res_metab_rate)) +
  (sal_test_Slope * mean(d$sal_test)) +
  (temp_test_and_mass_log_Slope * (24)* (d$mass_log)) +
  (temp_test_and_c_value_log_Slope * (24) * (d$c_value_log)) + 
  species_Random + 
  animal_Random

# 28 centigrade
d$pcrit_predicted28   <-  intercept+(temp_test_Slope * (28)) +
  (mass_log_Slope * (d$mass_log)) + 
  (c_value_log_Slope * (d$c_value_log)) + 
  (res_metab_rate_Slope * (d$res_metab_rate)) +
  (sal_test_Slope * mean(d$sal_test)) +
  (temp_test_and_mass_log_Slope * (28)* (d$mass_log)) +
  (temp_test_and_c_value_log_Slope * (28) * (d$c_value_log)) + 
  species_Random + 
  animal_Random

d$pcrit_predicted15   <-  d$pcrit_predicted15 + d$pcrit_kpa-d$pcrit_predicted
d$pcrit_predicted24   <-  d$pcrit_predicted24 + d$pcrit_kpa-d$pcrit_predicted
d$pcrit_predicted28   <-  d$pcrit_predicted28 + d$pcrit_kpa-d$pcrit_predicted

plot(d$pcrit_predicted15 ~ d$pcrit_predicted28)
abline(a=0,b=1)


d$Phylogeny_Random <- species_Random + animal_Random
pcritsummary<-aggregate(cbind(pcrit_predicted15,pcrit_predicted24,pcrit_predicted28,Phylogeny_Random)~animal,
                  data=d,median)

# Export the file with corrected Pcrit values per species
write.csv(pcritsummary,"standardized Pcrit for 174 species.csv", row.names = F)

#-------------------------------------------------------------------------------
#saving session information with all packages versions for reproducibility
#purposes
sink("./phylogenetic_hierarchical_models_R_session.txt")
sessionInfo()
sink()
################################################################################
############################ END OF SCRIPT #####################################
################################################################################