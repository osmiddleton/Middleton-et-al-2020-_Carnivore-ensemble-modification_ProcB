#---------------------------------------------------------------------------------------------------------------------------------------#
#       
#       ------------------------------------------------------------------------------#
#       
#       Middleton, O., Scharlemann, J.P.W. & Sandom, C., 2020. Global restrictions
#       in the trophic complexity of carnivorous mammal ensembles by modern humans.
#       Proceedings of the Royal Society B.
#       
#       ------------------------------------------------------------------------------#
#       
#       Analyses and figure production for the above manuscript investigating the
#       difference in carnivorous mammal community structure today, following the 
#       long-term impacts of modern-humans interacting with changing climates, compared
#       to a counter-factual scenario today (referred to throughout as the
#       'present-natural') had humans not migrated out of Africa (~100,000 years ago)
#       and interacted with changing climates to cause widespread extinctions, and range
#       contractions and restrictions of carnivorous mammals. 
#       
#       ------------------------------------------------------------------------------#
#              
#       Data accessibility: 
#       The main data used for this analyses is available from the Phylogenetic Atlas of Mammals
#       (Phylacine) Version 1.2, available at: https://github.com/MegaPast2Future/PHYLACINE_1.2
#       
#       Notes on data manipulation:
#       
#       Functional trait data: None. Access from Phylacine.
#       
#       Spatial data: Raster data in Phylacine was converted into a species-site (cell ID) matrix,
#       with attributes for each cell describing:
#       1. Continent
#       2. Biogeographic Realm
#       3. Present on the mainland during the Last Glacial Maximum (~110m below current sea level).
#       These species-site matrices were made available for analyses.
#       
#       The species-site matrices are available in DRYAD to re-run the
#       following analyses and produce the manuscript figures
#       
#       ------------------------------------------------------------------------------#
#       
#       Code content:
#       
#       # Main section
#       
#       1. Loading data, including functional traits and pre-processed species-site
#       matrices for (1) present-natural, and(2) current. 
#           - Data provided as a species-site matrix with raster cell information
#             including coordinates, and the biogeographical realm and continent
#             it falls into.
#              
#       2. Species selection from Phylacine. 
#           - Terrestrial species which consume >= 5% vertebrate meat.
#         
#       3. Perform Bayesian phylogenetic mixed models (BBPMM)
#           - Run using the mulTree package with MCMCglmm class.
#           - Select best models and predict on species list.
#           - Plot raw data with model predictions.
#           - This is performed on the whole species list, as well as smaller data
#             set using only those species which fit the systematic review of change.
#       
#       4. Calculate functional diversity of carnivore ensembles for (1) continents
#          (ESRI, 2002) and (2) biogeographical realms (as a sensitivity analysis).
#               a) Current, using two methods:
#                    (i)  Present-absence only 
#                    (ii) Including species range sizes in the analysis 
#               b) Present-Natural, using two methods:
#                    (i)  Present-absence only 
#                    (ii) Including species range sizes in the analysis 
#               
#       5. Global functional trait space.
#       
#       6. Continental ensemble structural difference between the two scenarios, using the
#          two methods (bubble-plots).
#       
#       7. Functional diversity differences between the present-natural and current.
#       
#       8. PC1 and PC2 shift in assemblage centroid between the two scenarios: Only
#          using continental assemblages
#       
#       
#       Senstivity analyses:
#       
#       8. Look at variation in change in FRic and FDis across the world by looking per
#          cell - No longer doing this
#       
#       
#       
#               
#---------------------------------------------------------------------------------------------------------------------------------------#

#————————————————————————————————————————————————————————————————————————————————————————————####
# ---- 1. Data preparation ----
#————————————————————————————————————————————————————————————————————————————————————————————####

# ---- Load libraries ----
# Install 'pacman' package if needed: install.packages("pacman")
library(pacman)

# Note: Some of these may be redundant throughout the script.
pacman::p_load(tidyverse, vegan, sf, raster, spatial.tools,FD,gridExtra,grid,viridis,ggfortify,
               scatterplot3d,scico,SDMTools,gtable, robustbase,caret,wosr,rgdal,bibliometrix,
               units,ggrepel,plotly,dunn.test,tiff,ggpubr,effects,compositions,lme4,nloptr,nlme,
               MuMIn,AICcmodavg,ape,devtools,ggplot2,pez,picante,DMwR,metafor,
               caTools,pscl,cowplot,merTools,MCMCglmm, mulTree) 
# Note: You might need to download some packages using 'devtools' and downloading directly from
# Github repositories.

#---- Prepare colour schemes ----
good.colours <- c("#FFCC99","#99CCFF","#CCCC99","#CC99CC","#FFCC00","#66CC99")
good.colours.2 <- c("#66CC99","#FFCC00","#CC99CC","#CCCC99","#99CCFF","#FFCC99")


#---- Load data ----

# Species-site matrices: current and present-natural

# Download these datasets from DRYAD.

# Note: These analyses were performed using Phylacine Version 1.2. If using this code for future analyses,
# downloading the most recent version of Phylacine would be best as some geographic ranges may have been
# updated.

# The files are too large to upload to GitHub, so save in a folder behind the current working directory
# in the directories: DRYAD Data > Phylacine_Range_Matrices to mensure this code runs.

# Load species-site matrices 
current.range.full <- read.csv("../DRYAD Data/Phylacine_Range_Matrices/Current_Species_Matrix.csv",
                               header = TRUE)
presnat.range.full <- read.csv("../DRYAD Data/Phylacine_Range_Matrices/PN_Species_Matrix.csv",
                               header = TRUE)

# Phylacine functional trait data (pre-downloaded and stored in GitHub)
# Note: Phylacine has since been updated and may have different functional traits now.

# Load trait data
traits <- read.csv("./Data/Phylacine_Traits/Trait_data.csv", header = TRUE)
traits <- traits[c(1:9,12,18,20:22,17)] # Subset traits 

#————————————————————————————————————————————————————————————————————————————————————————————####
# ---- 2. Species selection -----
#————————————————————————————————————————————————————————————————————————————————————————————####

# ---- Spatial filtering ----

# Species on islands attached the Last Glacial Maximum mainland (~110m below sea level)
current.range.tidy <- current.range.full %>%
                            filter(!is.na(Mainland))
presnat.range.tidy <- presnat.range.full %>%
                            filter(!is.na(Mainland))

# Remove cells which did not occur in continent shapefiles.
current.range.tidy <- current.range.tidy %>%
  filter(!is.na(Continent))

presnat.range.tidy <- presnat.range.tidy %>%
  filter(!is.na(Continent))

# Potential species occurring in continents and islands previously attached to the mainland.
potential.species <- cbind(presnat.range.tidy[c(1:8)],
                           subset(presnat.range.tidy[9:5839], 
                                  select = names(which(colSums(presnat.range.tidy[9:5839]) != 0))))
potential.species <- as.data.frame(colnames(potential.species[9:4943]))
colnames(potential.species) <- c("Binomial.1.2")

# ---- Remove humans ---- 
traits <- traits[!traits$Genus.1.2 == "Homo",] # Remove humans

# ---- Identify terrestrial carnivorous mammals ----
species <- traits[traits$Marine == 0,]
species <- species[species$Aerial == 0,]
species <- species[species$Freshwater == 0,]
species <- species[species$Terrestrial == 1,]
species <- species[species$Diet.Vertebrate >= 5,]

# Which species remaining fit our criteria
keep <- which(species$Binomial.1.2 %in% potential.species$Binomial.1.2)
species <- species[keep,]

# Taxonomic information of selected species.
species$Order.1.2 <- factor(species$Order.1.2)
levels(species$Order.1.2) # 15 Orders
length(species$Binomial.1.2) # 1081 total species

# How many species in each order?
(order <- species %>% group_by(Order.1.2) %>% summarize(sum = length(Order.1.2)))

# How many extinct extinct species?
count(species[species$IUCN.Status.1.2 == "EX" | # Extinct recent
                species$IUCN.Status.1.2 == "EP",]) # Extinct Pleistocene
# 31 extinct species
# 1050 extant species

# Bar plot of species in each category
# Make colour scale with 16 colours
colours <- c("#CCCCCC","#99CC00","#999999","#339900","#CC9900","#FFCC00","#FF6600","#993300",
             "#FFCCCC","#CC6699","#FF99FF","#663399","#9999CC","#3399FF","#CCFFFF","#003333")

# Update Order levels to show the number of species in each order.
species2 <- species
species2$Order.1.2 <- plyr::revalue(species2$Order.1.2, c("Afrosoricida" = "Afrosoricida (2)",
                                            "Carnivora" = "Carnivora (228)",
                                            "Cetartiodactyla" = "Cetartiodactyla (30)",
                                            "Cingulata" = "Cingulata (5)",
                                            "Dasyuromorphia" = "Dasyuromorphia (66)",
                                            "Didelphimorphia" = "Didelphimorphia (59)",
                                            "Diprotodontia" = "Diprotodontia (20)",
                                            "Eulipotyphla" = "Eulipotyphla (345)",
                                            "Lagomorpha" = "Lagomorpha (1)",
                                            "Macroscelidea" = "Macroscelidea (1)",
                                            "Paucituberculata" = "Paucituberculata (6)",
                                            "Peramelemorphia" = "Peramelemorphia (9)",
                                            "Primates" = "Primates (98)",
                                            "Rodentia" = "Rodentia (210)",
                                            "Scandentia" = "Scandentia (1)"))

# Interquartile range of percentage of vertebrate meat consumed in species' diets
quantile(species2$Diet.Vertebrate)
# Median = 30%
# Lower Quartile = 10%
# Upper Quartile = 50%

# ---- Figure S2: Species vertebrate consumption histogram ----
tiff("./Results/Supplementary materials/Figures/Species carnivory distribution.tif", height = 4500, width = 8150, res = 800)
ggplot(species2, aes(x = Diet.Vertebrate, colour = Order.1.2, fill = Order.1.2)) +
  geom_histogram(breaks=seq(0, 100, by=10)) +
  # scale_fill_brewer(palette = "Paired") +
  # scale_colour_brewer(palette = "Paired") +
  scale_x_continuous(name = "Vertebrates in diet (%)", breaks = seq(0,100, 20)) +
  scale_fill_manual(values = colours,name = "Order") +
  geom_vline(xintercept = 30, size = 1,alpha = .5) +
  geom_vline(xintercept = 10, size = 1, linetype = 2,alpha = .5) +
  geom_vline(xintercept = 50, size = 1, linetype = 2,alpha = .5) +
  scale_colour_manual(values = colours, name = "Order") +
  theme_bw() + xlab("Vertebrates in diet (%)") + ylab("Number of species") + 
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12)) 
dev.off()

#————————————————————————————————————————————————————————————————————————————————————————————####
# ---- 3. Bayesian phylogenetic mixed modelling (BPMM) ----
#————————————————————————————————————————————————————————————————————————————————————————————####

# Note: The results of MCMCglmms can vary each time, and so re-running these models would produce 
#       slightly different values of effect size and CIs to what is reported in the manuscript.

# An example of the model predictions from these models is available on DRYAD. 

# ---- Data preparation ----

# Get current and present-natural range sizes for all species
# Range size = Number of raster cells.
specs <- as.character(species$Binomial.1.2)
ranges <- matrix(ncol = 2, nrow = nrow(species))

for(i in 1:length(specs)) {
  x <- specs[[i]]
  a <- which(colnames(current.range.full) %in% x)
  y <- sum(current.range.full[c(a)])
  z <- sum(presnat.range.full[c(a)])
  ranges[i,1] <- z
  ranges[i,2] <- y
  print(i)
}

# Add range size to the Phylacine Trait Matrix.
ranges <- as.data.frame(ranges)
colnames(ranges) <- c("PresnatRange","CurrentRange")
species2 <- cbind(species, ranges)

# Add proportional loss
species2$Prop_RangeContraction <- (1 - (species2$CurrentRange/species2$PresnatRange))

# Tidy predictor variables for models (log, if necessary, and scale variables)
species2$log.bodymass <- log10(species2$Mass.g) # Log body mass
scaled.variables <- scale(species2[c(19,13)]) # Scale diet vertebrate and body mass
colnames(scaled.variables) <- c("scaled.mass","scaled.diet")
species2 <- cbind(species2, scaled.variables)
colnames(species2)
species2 <- species2[c(1:3,11,13,16:21)]

# Four species have had range expansions:
# Coyote, North African white-toothed shrew, least weasel, lesser white-toothed shrew.
# Convert proportional range difference to 0 (i.e. NO range loss)
species2[species2$Prop_RangeContraction < 0,]$Prop_RangeContraction <- 0
species2 <- species2[!is.na(species2$Binomial.1.2),]


# ---- Identify non-systematically reviewed species ----

# In Faurby & Svenning (2015), species current IUCN ranges were systematically
# inspected for creation of a present-natural range which were:

# 1. Large-bodied species (>=1kg) 
species2$LargeSpecies <- NA
species2$mass.g <- 10^(species2$log.bodymass)
species2[species2$mass.g >= 1000,]$LargeSpecies <- 1

# 2. Species which were classified by the IUCN as vulnerable, endangered, critically endangered,
#    extinct, extinct in the wild, or data deficient.
species2$systematic.species <- NA
species2[species2$mass.g >= 1000,]$systematic.species <- 1
species2[!species2$IUCN.Status.1.2 == "LC" &
         !species2$IUCN.Status.1.2 == "NT",]$systematic.species <- 1

# Small-bodied species (<1kg) and those that are currently not classified as being
# threatened were assumed to currently be at their full potential geographic range.

# All analyses were performed twice to ensure our results were still consistent for only
# those that were able to have differences in geographic ranges between the current and
# present-natural.


# ---- Bayesian phylogenetic mixed models (BPMM) ----

# ---- Prepare data ----
colnames(species2)
spec <- species2
spec
nrow(spec %>% filter(IUCN.Status.1.2 == "EP")) # 26 species
nrow(spec %>% filter(IUCN.Status.1.2 == "EX")) # 5 species
# 31 extinct species in total.


# Download phylogenies from Phylacine and save in the same folder as the folder containing 
# DRYAD Data folder.

# Load species' phylogeny.
trees <- read.nexus("../Complete_phylogeny.nex")
tree_samp <- sample(trees, size = 100) # This will be random every time.

# Prune tree to only the species list being used in this study
pruned.forest <- lapply(tree_samp, 
                        function(x) {
                          drop.tip(x, x$tip.label[!x$tip.label %in% levels(factor(spec$Binomial.1.2))])
                        }
)
class(pruned.forest) <- "multiPhylo"

# Check a single tree
plot(pruned.forest[[1]]) 

## Number of tips in both
unlist(lapply(pruned.forest, Ntip)) # Should be 1,081 tips for all species.

# Change species which have undergone range exapnsions as having
# PN ranges as an equal size to the current range.
expand <- spec[which(spec$PresnatRange < spec$CurrentRange),]
expand$CurrentRange <- expand$PresnatRange
non.expand <- spec[which(spec$PresnatRange >= spec$CurrentRange),]
# rbind the two dataframes back together:
spec <- rbind(non.expand,expand)

# Identify range loss

spec$Rangeloss <- spec$PresnatRange - spec$CurrentRange

# The multinomial2 model requires predictor variables as integers
is.integer(spec$CurrentRange) # TRUE
is.integer(spec$Rangeloss) # FALSE - not sure why but change to be integer...
hist(spec$Rangeloss)
spec$Rangeloss <- as.integer(spec$Rangeloss)
is.integer(spec$Rangeloss) # TRUE

# ---- Main results: Full species list ----

# Running MCMCglmm using multiple sampled trees - takes into consideration phylogenetic uncertainty
# Save as mulTree object
mulTree_data <- as.mulTree(data = spec,
                           tree = pruned.forest[1:100],
                           taxa = "Binomial.1.2")

str(mulTree_data)
show(mulTree_data)

# MCMC parameters
my_parameters <- c(200000, # Number of MCMCglmm iterations,
                   100, # Thinning interval (Default = 10)
                   10000) # Burnin period (Discard first 10,000 iterations)
                                    
## MCMCglmm priors
## Using default inverse-Wishart uninformative priors for:
## (R) the model residuals
## (G) the random effects [phylogenetic covariance matrix]
my_priors <- list(R = list(V = 1/2, nu = 0.002),            # Covariance matrix of model residuals
                  G = list(G1 = list(V = 1/2, nu = 0.002))) # Random effects priors: V = scale matrix, nu = df
                                                            # nu = 0.002 
# These priors were selected based on papers by Healy.

## This object is classified as "mulTree" and contains different elements to be
## passed to the mulTree function
class(mulTree_data) ; names(mulTree_data)


# ---- Hierarchical model running ----

# Note: These models do take hours to run.

# Create list to store the DIC outputs
DIC.list <- list()

# 1. Global model
mulTree(cbind(Rangeloss, CurrentRange) ~ scaled.mass * scaled.diet, 
        mulTree.data = mulTree_data,
        priors = my_priors,
        parameters = my_parameters, 
        output = "./Results/MCMCglmm chain outputs/Model1/range_loss", # May have added on an extra number to make sure it can be separated from the others
        family = "multinomial2",
        ESS = 1000, chains = 2) # Using the default assessment of effective sample size

# Get average DIC value and interquartile range from all of the models
DICs <- matrix(ncol = 1, nrow = length(pruned.forest[1:100]))
setwd("./Results/MCMCglmm chain outputs/Model1")
for(i in 1:length(pruned.forest[1:100])) {
  x <- read.mulTree(paste0("range_loss-tree",i,"_chain1"),
                    model = TRUE)
  DICs[i,1] <- x$DIC
  print(i)
}
DIC.list[[1]] <- DICs

# Re-set working directory.
setwd("../../../")


# 2. No interaction model
mulTree(cbind(Rangeloss, CurrentRange) ~ scaled.mass + scaled.diet, 
        mulTree.data = mulTree_data,
        priors = my_priors,
        parameters = my_parameters, 
        output = "./Results/MCMCglmm chain outputs/Model2/range_loss",
        family = "multinomial2",
        ESS = 1000, chains = 2)

# Get average DIC value and interquartile range from all of the models
DICs <- matrix(ncol = 1, nrow = length(pruned.forest[1:100]))
setwd("./Results/MCMCglmm chain outputs/Model2")
for(i in 1:length(pruned.forest[1:100])) {
  x <- read.mulTree(paste0("range_loss-tree",i,"_chain1"),
                    model = TRUE)
  DICs[i,1] <- x$DIC
  print(i)
}

DIC.list[[2]] <- DICs

# Re-set working directory.
setwd("../../../")


# 3. Only  bodymass
mulTree(cbind(Rangeloss, CurrentRange) ~ scaled.mass, 
        mulTree.data = mulTree_data,
        priors = my_priors,
        parameters = my_parameters, 
        output = "./Results/MCMCglmm chain outputs/Model3/range_loss",
        family = "multinomial2",
        ESS = 1000, chains = 2)

# Get average DIC value and interquartile range from all of the models
DICs <- matrix(ncol = 1, nrow = length(pruned.forest[1:100]))
setwd("./Results/MCMCglmm chain outputs/Model3")
for(i in 1:length(pruned.forest[1:100])) {
  x <- read.mulTree(paste0("range_loss-tree",i,"_chain1"),
                    model = TRUE)
  DICs[i,1] <- x$DIC
  print(i)
}
DIC.list[[3]] <- DICs

# re-set working directory
setwd("../../../")


# 4. Only diet
mulTree(cbind(Rangeloss, CurrentRange) ~ scaled.diet, 
        mulTree.data = mulTree_data,
        priors = my_priors,
        parameters = my_parameters, 
        output = "./Results/MCMCglmm chain outputs/Model4/range_loss2",
        family = "multinomial2",
        ESS = 1000, chains = 2)

# Get average DIC value and interquartile range from all of the models
DICs <- matrix(ncol = 1, nrow = length(pruned.forest[1:100]))
setwd("./Results/MCMCglmm chain outputs/Model4")
for(i in 1:length(pruned.forest[1:100])) {
  x <- read.mulTree(paste0("range_loss-tree",i,"_chain1"),
                    model = TRUE)
  DICs[i,1] <- x$DIC
  print(i)
}

DIC.list[[4]] <- DICs

# The model with the interaction between diet and body mass was the best-supported model.
# The DIC value can be interpreted in a similar way to an AIC value
lapply(DIC.list, quantile)
lapply(DIC.list,mean)


# ---- Diagnostic plots -----

# re-set working directory
setwd("../../..")

setwd("./Results/MCMCglmm chain outputs/Model1")
chain.1 <- read.mulTree(paste0("range_loss-tree","20","_chain1"),
                  model = TRUE)
chain.2 <- read.mulTree(paste0("range_loss-tree","20","_chain2"),
                        model = TRUE)

# Diagnostic check for how well the chains converged on the same result for posterior
# distributions of the fixed effects
sol <- bind_rows(as.data.frame(chain.1$Sol[, 1:4]), 
                 as.data.frame(chain.2$Sol[, 1:4]))

sol["chain"] <- gl(2, 1900) # Number of chains, number of iterations(?)
sol["sample"] <- rep(1:1900, 2)
sol <- gather(sol, key = "variable", value = "value", -chain, -sample)

left <- ggplot(sol, aes(x = sample, y = value, col = chain)) + # Should look like a fuzzy caterpillar
  geom_line() + 
  facet_wrap(~ variable, scales = "free", nrow = 4) +
  theme(legend.position="none") + 
  ylab("")
right <- ggplot(sol, aes(x = value, col = chain)) + # Distribution of fixed effects parameter estimates
  geom_density() +
  geom_rug() +
  facet_wrap(~ variable, scales = "free", nrow = 4) +
  theme(legend.position="none") + 
  labs(x = "", y = "")
grid.arrange(left, right, nrow = 1)

# Need to also get the distribution of the random effects parameters
sol <- bind_rows(as.data.frame(chain.1$VCV[, 1:2]), 
                 as.data.frame(chain.2$VCV[, 1:2]))

sol["chain"] <- gl(2, 1900) # Number of chains, number of iterations(?)
sol["sample"] <- rep(1:1900, 2)
sol <- gather(sol, key = "variable", value = "value", -chain, -sample)

left <- ggplot(sol, aes(x = sample, y = value, col = chain)) + # Should look like a fuzzy caterpillar
  geom_line() + 
  facet_wrap(~ variable, scales = "free", nrow = 4) +
  theme(legend.position="none") + 
  ylab("")
right <- ggplot(sol, aes(x = value, col = chain)) + # Distribution of fixed effects parameter estimates
  geom_density() +
  geom_rug() +
  facet_wrap(~ variable, scales = "free", nrow = 4) +
  theme(legend.position="none") + 
  labs(x = "", y = "")
grid.arrange(left, right, nrow = 1)

# Checking convergence for our fixed factors (close to one is good)
class(chain.2$Sol)
gelman.diag(mcmc.list(as.mcmc(chain.1$Sol[, 1:4]),
                      as.mcmc(chain.2$Sol[, 1:4])),
            autoburnin = FALSE)

# Checking convergence for our random terms
gelman.diag(mcmc.list(chain.1$VCV,
                      chain.2$VCV),
            autoburnin = FALSE)

### Checking only the first chain
effectiveSize(chain.1$Sol[, 1:4])
effectiveSize(chain.1$VCV)

# acf plot for the fixed estimates
acf(chain.1$Sol[, 1], lag.max = 20)
acf(chain.1$Sol[, 2], lag.max = 20)
acf(chain.1$Sol[, 3], lag.max = 20)
acf(chain.1$Sol[, 4], lag.max = 20)

# acf plot for the first random term in our model (the phyl term)
acf(chain.1$VCV[, 1], lag.max = 20)
acf(chain.1$VCV[, 2], lag.max = 20)


# ---- Model predictions on raw data ----

# To get predictions for on raw data, follow this code.
# Skip to line 645 if you want to recreate the figure and look at model output example.

# Reading all the models to perform the MCMCglmm analysis on multiple trees
all_models <- read.mulTree("range_loss") 
str(all_models)
class(all_models)

# Summarising the results by estimating the highest density regions 
# and their associated 95 and 50 confidence intervals (default) 
summarised_results <- summary(all_models, use.hdr = FALSE)
summarised_results

# Predict range loss using the model outputs:

# To recreate the figure for the manuscript, skip to line 641.

# Create new dataframe for predicting
colnames(spec)
new.df <- spec[c(1,7,8,10:11,15)]
colnames(new.df)[1] <- "animal" # Change Binomial.1.2 to be 'animal' for predicting.

# Create matrix to store model predictions:
model.predicts <- list()

# Get model estimates from each of the models (takes a while to run)
for(i in 1:length(tree_samp)) {
one_model <- read.mulTree(paste0("range_loss-tree",i,"_chain1"),
                          model = TRUE) ## This model is a normal MCMCglmm object that has been ran on one single tree class(one_model) ; names(one_model)

# Get model estimates from each of the models
test <- predict(one_model, new.df,
                type = "response", interval = "confidence")
model.predicts[[i]] <- test
print(i)
}

predict.matrix <- matrix(nrow = 1081, ncol = 100)
lower.matrix <- matrix(nrow = 1081, ncol = 100)
upper.matrix <- matrix(nrow = 1081, ncol = 100)

# Get average for the 100 models
for(i in 1:100) {
  
  x <- model.predicts[[i]]
  predict.matrix[,i] <- x[,1]
  lower.matrix[,i] <- x[,2]
  upper.matrix[,i] <- x[,3]
  print(i)
}

spec$PredictedLoss <- apply(predict.matrix,1,mean)
spec$PredictedLossLower <- apply(lower.matrix,1,mean)
spec$PredictedLossUpper <- apply(upper.matrix,1,mean)
spec$ProportionalPredictedLoss <- spec$PredictedLoss/spec$PresnatRange
spec$ProportionalPredictedLossLower <- spec$PredictedLossLower/spec$PresnatRange
spec$ProportionalPredictedLossUpper <- spec$PredictedLossUpper/spec$PresnatRange

# If it needs to be saved, save the output
#write.csv(spec,"../Model1 predictions on raw data.csv")

# Read the model predictions in
setwd("../../../")
spec <- read.csv("../DRYAD Data/BBPMM_Predictions/Model1 predictions on raw data.csv", header = TRUE)

# Create new factor of carnivore categories
spec$vert.cat <- NA
spec[spec$Diet.Vertebrate >= 70,]$vert.cat <- "1" #  Hypercarnivore
spec[spec$Diet.Vertebrate < 70,]$vert.cat <- "2" # Non hypercarnivore

# ---- Figure 1: Predictions of propotional range difference ----

# Creating two figures separately:
 
# 1. Predictions and raw data for hypercarnivores (>=70%)
# 2. Predictions and raw data for non-hypercarnivores (5-69%)

# 1. Hypercarnivores (>=70%)
(g1 <- ggplot(data = spec[spec$vert.cat == "1",]) +
  geom_point(aes(x = log.bodymass, y = Prop_RangeContraction*-1),
             alpha = 0.3) +
  geom_smooth(aes(x = log.bodymass, y = ProportionalPredictedLossLower*-1),
              method = "loess", se = FALSE, linetype = 0, size = 0.5, colour = "black") +
  geom_smooth(aes(x = log.bodymass, y = ProportionalPredictedLossUpper*-1),
              method = "loess", se = FALSE, linetype = 0, alpha = 0.1, size = 0.5, colour = "black") +
  #scale_colour_brewer(palette = "Paired") + 
  scale_x_continuous(breaks = c(0:6),
                     labels = c("0","0.01", "0.1", "1", "10", "100", "1000"),
                     limits = c(0,6)) +
  xlab("Body mass (kg)") +
  ylab("Proportional range difference") +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()))

# build plot object for rendering 
gg1 <- ggplot_build(g1)

# extract data for the loess lines from the 'data' slot
df2 <- data.frame(x = gg1$data[[2]]$x,
                  ymin = gg1$data[[2]]$y,
                  ymax = gg1$data[[3]]$y) 

# use the loess data to add the 'ribbon' to plot 
(gg1 <- g1 +
  geom_smooth(aes(x = log.bodymass, y = ProportionalPredictedLoss*-1),
              method = "loess", se = FALSE, colour = "black") +
  geom_ribbon(data = df2, aes(x = x, ymin = ymin, ymax = ymax),
              fill = "grey", alpha = 0.4, linetype = 0))

# 2. Non-hypercarnivores (5-69%)
(g2 <- ggplot(data = spec[spec$vert.cat == "2",]) +
  geom_point(aes(x = log.bodymass, y = Prop_RangeContraction*-1),
             alpha = 0.3) +
  geom_smooth(aes(x = log.bodymass, y = ProportionalPredictedLossLower*-1),
              method = "loess", se = FALSE, linetype = 0, size = 0.5, colour = "black") +
  geom_smooth(aes(x = log.bodymass, y = ProportionalPredictedLossUpper*-1),
              method = "loess", se = FALSE, linetype = 0, alpha = 0.1, size = 0.5, colour = "black") +
  #scale_colour_brewer(palette = "Paired") + 
  scale_x_continuous(breaks = c(0:6),
                     labels = c("0","0.01", "0.1", "1", "10", "100", "1000"),
                     limits = c(0,6)) +
  xlab("Body mass (kg)") +
  ylab("Proportional range difference") +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()))

# build plot object for rendering 
gg2 <- ggplot_build(g2)

# extract data for the loess lines from the 'data' slot
df2 <- data.frame(x = gg2$data[[2]]$x,
                  ymin = gg2$data[[2]]$y,
                  ymax = gg2$data[[3]]$y) 

# use the loess data to add the 'ribbon' to plot 
(gg2 <- g2 +
  geom_smooth(aes(x = log.bodymass, y = ProportionalPredictedLoss*-1),
              method = "loess", se = FALSE, colour = "black") +
  geom_ribbon(data = df2, aes(x = x, ymin = ymin, ymax = ymax),
              fill = "grey", alpha = 0.4, linetype = 0))



pdf("./Results/Figures/BPMM_Results.pdf",
     height = 4, width = 9)
plot_grid(gg2, gg1,
          align = "v", labels = "auto")
dev.off()


#————————————————————————————————————————————————————————————————————————————————————————————####
# ----  4. Functional diversity of continental carnivorous mammal enemsbles ----
#————————————————————————————————————————————————————————————————————————————————————————————####
#      Under two scenarios:
#               a) Current, using two methods:
#                    (i)  Present-absence only 
#                    (ii) Including species range sizes in the analysis
#                     
#               b) Present-Natural, using two methods:
#                    (i)  Present-absence only 
#                    (ii) Including species range sizes in the analysis 

# Load in continent shape files (Behrmann Projection) using sf
continents <- st_read("./Data/continentbehrmann/continentprojected.shp")
plot(continents$geometry)

# Tidy up community matrices
c <- which(colnames(current.range.tidy[c(9:5839)]) %in% species2$Binomial.1.2)
current.output <- cbind(current.range.tidy[2:8],
                        current.range.tidy[9:5839][c(c)])

c <- which(colnames(presnat.range.tidy[c(9:5839)]) %in% species2$Binomial.1.2)
present.natural.output <- cbind(presnat.range.tidy[2:8],presnat.range.tidy[9:5839][c(c)])


#---- Carnivore continental range size ----
# Combine current and present-natural outputs
global.scenarios.merged <- rbind(current.output,
                                 present.natural.output)
global.scenarios.merged$Scenario <- as.factor(global.scenarios.merged$Scenario)
global.scenarios.merged$Continent <- as.factor(global.scenarios.merged$Continent)
global.scenarios.merged$realm <- as.factor(global.scenarios.merged$realm)
str(global.scenarios.merged)

colnames(global.scenarios.merged)

scenario <- levels(global.scenarios.merged$Scenario)
continents <- levels(global.scenarios.merged$Continent)

# Sum all species per continent 
carnivore.continent.range <- data.frame()

for (i in 1:length(scenario)) {
  x <- global.scenarios.merged[global.scenarios.merged$Scenario == scenario[[i]],]
  continent.range <- matrix(ncol = length(colnames(current.output[-c(1:7)])), nrow = 6)
  
  for (j in 1:length(continents)) {
    x <- x[!is.na(x$Continent),]
    y <- x[x$Continent == continents[[j]],]
    y <- y[c(8:length(current.output))]
    z <- colSums(y)
    continent.range[j,] <- z
    print(j)
  }
  
  continent.range <- as.data.frame(continent.range)
  carnivore.continent.range <- rbind(continent.range, carnivore.continent.range)
  print(i)
}


# Tidy up carnivore continent range size dataframe
colnames(carnivore.continent.range) <- colnames(current.output[c(8:length(current.output))])
carnivore.continent.range$Scenario <- c("Present-Natural","Present-Natural","Present-Natural","Present-Natural","Present-Natural","Present-Natural",
                                        "Current","Current","Current","Current","Current","Current")
carnivore.continent.range$Continent <- c(continents,continents)
carnivore.continent.range <- carnivore.continent.range[c((length(carnivore.continent.range)-1):(length(carnivore.continent.range)),
                                                         1:(length(carnivore.continent.range)-2))]


# ----  a) Functional diversity WITH geographic range size----

# 1. Community Matrix
row.names(carnivore.continent.range) <- paste(carnivore.continent.range$Scenario, carnivore.continent.range$Continent, sep = "_")
carnivore.continent.range <- carnivore.continent.range[-c(1:2)]
carnivore.not.on.continents <- subset(carnivore.continent.range, select = names(which(colSums(carnivore.continent.range) == 0))) # which species are we losing here
carnivore.not.on.continents <- as.data.frame(colnames(carnivore.not.on.continents))
colnames(carnivore.not.on.continents) <- "Binomial.1.2"
carnivore.not.on.continents <- merge(carnivore.not.on.continents,species,by = "Binomial.1.2") # Species being lost due to ranges not falling into continental shapefiles...

carnivore.continent.range <- subset(carnivore.continent.range, select = names(which(colSums(carnivore.continent.range) != 0))) 
str(carnivore.continent.range)

# Quickly check a species:
carnivore.continent.range$Urocyon_cinereoargenteus

# 2. Traits
new.species <- as.data.frame(colnames(carnivore.continent.range))
colnames(new.species)[1] <- "Binomial.1.2"
traits.sorted <- merge(species, new.species)
row.names(traits.sorted) <- traits.sorted$Binomial.1.2

str(traits.sorted)
row.names(traits.sorted) <- traits.sorted$Binomial.1.2
traits.sorted <- traits.sorted[c(10,12:14)]

# Transform body mass 
traits.sorted$Mass.g <- log10(traits.sorted$Mass.g)

# Scale continuous traits
traits.sorted <- scale(traits.sorted)

# Put matrix and traits into a list to call dbFD
joined.world.function <- list(traits.sorted, carnivore.continent.range)
names(joined.world.function)[1] <- "traits"
names(joined.world.function)[2] <- "presence"
str(joined.world.function$traits)
str(joined.world.function$presence)

# Set weights
w <- c(1,1/3,1/3,1/3)

# Call dbFD
global.FD.continents.GR <- dbFD(joined.world.function$traits, 
                                joined.world.function$presence, 
                                w,
                                m = 4,
                                stand.FRic = TRUE, # Standardise functional richness of sites by the global trait space
                                calc.CWM = TRUE,
                                print.pco = TRUE) # 0.73 R2 value

# Extract the eigenvalues
x <- global.FD.continents.GR$x.values

# Calculate the importance of the first four principal coordinate axes:
sum(x[1:4])/sum(x)

# Calculate the importance of the first four axes individually:
x[1]/sum(x) # 39.5%
x[2]/sum(x) # 19.2%
x[3]/sum(x) # 8.5%
x[4]/sum(x) # 5.8%

#   Extract  results of differences between scenarios at the level of continents
FD.Continents.GR <- as.data.frame(cbind(global.FD.continents.GR$nbsp, global.FD.continents.GR$FRic,global.FD.continents.GR$FDis))

FD.Continents.GR$Scenario <- c("Present-Natural","Present-Natural","Present-Natural","Present-Natural","Present-Natural","Present-Natural",
                               "Current","Current","Current","Current","Current","Current")

FD.Continents.GR$Region.names <- c(continents,continents)
str(FD.Continents.GR)
FD.Continents.GR$Scenario <- as.factor(FD.Continents.GR$Scenario)
FD.Continents.GR$Scenario <- relevel(FD.Continents.GR$Scenario, "Present-Natural")
FD.Continents.GR$Region.names <- as.factor(FD.Continents.GR$Region.names)
FD.Continents.GR$Method <- "Geographic Range"
FD.Continents.GR$Region <- "Continent"
colnames(FD.Continents.GR)[1:3] <- c("SR","FRic","FDis")


# ----- b) Functional diversity WITHOUT geographic range size ----
# 1. change community matrix
str(carnivore.continent.range)
carnivore.continent.presence <- carnivore.continent.range
carnivore.continent.presence[carnivore.continent.presence > 1] <- 1
str(carnivore.continent.presence)

# 2. Traits
new.species <- as.data.frame(colnames(carnivore.continent.presence))
colnames(new.species)[1] <- "Binomial.1.2"
traits.sorted <- merge(species, new.species)
row.names(traits.sorted) <- traits.sorted$Binomial.1.2

str(traits.sorted)
row.names(traits.sorted) <- traits.sorted$Binomial.1.2
traits.sorted <- traits.sorted[c(10,12:14)]
str(traits.sorted)

# Transform body mass
traits.sorted$Mass.g <- log10(traits.sorted$Mass.g)

# Scale traits
traits.sorted <- scale(traits.sorted)

# PUT MATRIX AND TRAITS INTO LIST FOR WHICH TO CALL FD
joined.world.function <- list(traits.sorted, carnivore.continent.presence)
names(joined.world.function)[1] <- "traits"
names(joined.world.function)[2] <- "presence"
str(joined.world.function$traits)
str(joined.world.function$presence)

# SET WEIGHTS
w <- c(1,1/3,1/3,1/3)

# CALL FD - OTHER PARAMETERS ARE OPTIONS: REMEMBER! STAND.FRIC TO COMPARE ALL COMMUNITIES TO GLOBAL POOL FOR COMPARABILITY.
global.FD.continents.no.GR <- dbFD(joined.world.function$traits, 
                                   joined.world.function$presence, w, m = 4,
                                   stand.FRic = TRUE,
                                   calc.CWM = TRUE,
                                   print.pco = TRUE) # 73%

#     Extract and look at results of differences between scenarios at the level of continents
FD.Continents.no.GR <- as.data.frame(cbind(global.FD.continents.no.GR$nbsp,
                                           global.FD.continents.no.GR$FRic,
                                           global.FD.continents.no.GR$FDis))

FD.Continents.no.GR$Scenario <- c("Present-Natural","Present-Natural","Present-Natural","Present-Natural","Present-Natural","Present-Natural",
                                  "Current","Current","Current","Current","Current","Current")

FD.Continents.no.GR$Region.names <- continents
str(FD.Continents.no.GR)
FD.Continents.no.GR$Scenario <- as.factor(FD.Continents.no.GR$Scenario)
FD.Continents.no.GR$Scenario <- relevel(FD.Continents.no.GR$Scenario, "Present-Natural")
FD.Continents.no.GR$Region.names <- as.factor(FD.Continents.no.GR$Region.names)
FD.Continents.no.GR$Method <- "P/A"
FD.Continents.no.GR$Region <- "Continent"
colnames(FD.Continents.no.GR)[1:3] <- c("SR","FRic","FDis")


#————————————————————————————————————————————————————————————————————————————————————————————####
# ---- 5. Global functional trait space ----
#————————————————————————————————————————————————————————————————————————————————————————————####
# Get taxonomic information
traits <- read.csv("./Data/Phylacine_Traits/Trait_data.csv", header = TRUE)

# Get the PCoA axes
functional.axes <- global.FD.continents.GR$x.axes[c(1:4)] # All of these will be the same, regardless of the scenario
functional.axes$Binomial.1.2 <- row.names(functional.axes)
functional.axes <- functional.axes[c(5,1:4)]
functional.axes.and.traits <- merge(functional.axes, traits, by.y = "Binomial.1.2")

# ---- Figure S6: Global trait space ----
(global.plot <- ggplot(functional.axes.and.traits, aes(x=A1, y = A2)) +
  geom_point(aes(fill = Order.1.2, label = Binomial.1.2), pch = 21, colour = "black") +
  geom_density_2d(size = .1, alpha = .6, colour = "black") +
  theme_classic() + 
  scale_fill_manual(values = colours) +
  guides(colour = guide_legend(title="Order", title.position = "top")) +
  theme(legend.text=element_text(size=12),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14)) +
  expand_limits(y = c(-0.4,0.5), x = c(-0.55,0.6)) +
  xlab("A1 (39.5%)") + ylab("A2 (19.2%)"))

tiff("./Results/Figures/Global Trait Space_Updated.tif", height = 4500, width = 6000, res = 800)
global.plot
dev.off()
ggplotly(global.plot)

# ---- Figure S5: Axes correlations with traits ----
# Correlations of A1, A2, A3, and A4 with diet traits
summary(lm(data = functional.axes.and.traits, A1~log10(Mass.g)))
cor.test(functional.axes.and.traits$A1,log10(functional.axes.and.traits$Mass.g))
A1.mass <- ggplot(functional.axes.and.traits,aes(x=A1,y=log10(Mass.g))) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("") + ylab(~Log[10]~'(body mass)')

summary(lm(data = functional.axes.and.traits, A1~Diet.Plant))
cor.test(functional.axes.and.traits$A1,functional.axes.and.traits$Diet.Plant)
A1.plant <- ggplot(functional.axes.and.traits,aes(x=A1,y=Diet.Plant)) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("") + ylab("Plant (%)")

summary(lm(data = functional.axes.and.traits, A1~Diet.Vertebrate))
cor.test(functional.axes.and.traits$A1,functional.axes.and.traits$Diet.Vertebrate)
A1.vert <- ggplot(functional.axes.and.traits,aes(x=A1,y=Diet.Vertebrate)) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("") + ylab("Vertebrate (%)")

summary(lm(data = functional.axes.and.traits, A1~Diet.Invertebrate))
cor.test(functional.axes.and.traits$A1,functional.axes.and.traits$Diet.Invertebrate)
A1.invert <- ggplot(functional.axes.and.traits,aes(x=A1,y=Diet.Invertebrate)) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("A1 (39.5%)") + ylab("Invertebrate(%)")

summary(lm(data = functional.axes.and.traits, A2~log10(Mass.g)))
cor.test(functional.axes.and.traits$A2,log10(functional.axes.and.traits$Mass.g))
A2.mass <- ggplot(functional.axes.and.traits,aes(x=A2,y=log10(Mass.g))) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("") + ylab("")

summary(lm(data = functional.axes.and.traits, A2~Diet.Plant))
cor.test(functional.axes.and.traits$A2,functional.axes.and.traits$Diet.Plant)
A2.plant <- ggplot(functional.axes.and.traits,aes(x=A2,y=Diet.Plant)) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("") + ylab("")

summary(lm(data = functional.axes.and.traits, A2~Diet.Vertebrate))
cor.test(functional.axes.and.traits$A2,functional.axes.and.traits$Diet.Vertebrate)
A2.vert <-ggplot(functional.axes.and.traits,aes(x=A2,y=Diet.Vertebrate)) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("") + ylab("") 

summary(lm(data = functional.axes.and.traits, A2~Diet.Invertebrate))
cor.test(functional.axes.and.traits$A2,functional.axes.and.traits$Diet.Invertebrate)
A2.invert <- ggplot(functional.axes.and.traits,aes(x=A2,y=Diet.Invertebrate)) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("A2 (19.2%)") + ylab("")

summary(lm(data = functional.axes.and.traits, A3~log10(Mass.g)))
cor.test(functional.axes.and.traits$A3, log10(functional.axes.and.traits$Mass.g))
A3.mass <- ggplot(functional.axes.and.traits,aes(x=A3,y=log10(Mass.g))) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("") + ylab("")

summary(lm(data = functional.axes.and.traits, A3~Diet.Plant))
cor.test(functional.axes.and.traits$A3, functional.axes.and.traits$Diet.Plant)
A3.plant <- ggplot(functional.axes.and.traits,aes(x=A3,y=Diet.Plant)) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("") + ylab("")

summary(lm(data = functional.axes.and.traits, A3~Diet.Vertebrate))
cor.test(functional.axes.and.traits$A3, functional.axes.and.traits$Diet.Vertebrate)
A3.vert <- ggplot(functional.axes.and.traits,aes(x=A3,y=Diet.Vertebrate)) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("") + ylab("")

summary(lm(data = functional.axes.and.traits, A3~Diet.Invertebrate))
cor.test(functional.axes.and.traits$A3, functional.axes.and.traits$Diet.Invertebrate)
A3.invert <-ggplot(functional.axes.and.traits,aes(x=A3,y=Diet.Invertebrate)) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("A3 (8.5%)") + ylab("")

summary(lm(data = functional.axes.and.traits, A4~log10(Mass.g)))
cor.test(functional.axes.and.traits$A4, log10(functional.axes.and.traits$Mass.g))
A4.mass<- ggplot(functional.axes.and.traits,aes(x=A4,y=log10(Mass.g))) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("") + ylab("")

summary(lm(data = functional.axes.and.traits, A4~Diet.Plant))
cor.test(functional.axes.and.traits$A4, functional.axes.and.traits$Diet.Plant)
A4.plant <- ggplot(functional.axes.and.traits,aes(x=A4,y=Diet.Plant)) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("") + ylab("")

summary(lm(data = functional.axes.and.traits, A4~Diet.Vertebrate))
cor.test(functional.axes.and.traits$A4, functional.axes.and.traits$Diet.Vertebrate)
A4.vert <- ggplot(functional.axes.and.traits,aes(x=A4,y=Diet.Vertebrate)) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("") + ylab("")

summary(lm(data = functional.axes.and.traits, A4~Diet.Invertebrate))
cor.test(functional.axes.and.traits$A4, functional.axes.and.traits$Diet.Invertebrate)
A4.invert <- ggplot(functional.axes.and.traits,aes(x=A4,y=Diet.Invertebrate)) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("A4 (5.8%)") + ylab("")

tiff("./Results/Figures/As_Traits_Update.tif", height = 6000, width = 6000, res = 800)
grid.arrange(A1.mass, A2.mass, A3.mass, A4.mass,
             A1.plant,A2.plant, A3.plant, A4.plant,
             A1.vert, A2.vert, A3.vert, A4.vert,
             A1.invert, A2.invert, A3.invert, A4.invert,
             nrow = 4, ncol = 4)
dev.off()

# A1 = Larger-bodied species (+ve) or smaller-bodied species (-ve)
# A2 = Carnivores (+ve) or Herbivores (-ve)
# A3 = larger insectivores (+ve) or medium/small carnivorous/herbivorous (-ve)?
# A4 = Medium-bodied with mixed diets

#————————————————————————————————————————————————————————————————————————————————————————————####
# ----  6. Continental trait space and centroid identity ----
#————————————————————————————————————————————————————————————————————————————————————————————####
#   Produce two figures:
#    a) continental with geographic range
#    b) continental without geographic range

# ----  1. Continental with geographic range ----

# Now to merge it and plot it in trait space
str(functional.axes.and.traits)
functional.axes.and.traits <- functional.axes.and.traits[c(1:7)]

# Get species range per continent
str(carnivore.continent.range)
carnivore.continent.range$Scenario <- c("Present-Natural","Present-Natural","Present-Natural","Present-Natural","Present-Natural","Present-Natural",
                                        "Current","Current","Current","Current","Current","Current")
carnivore.continent.range$Continent <- c(continents,continents)
carnivore.continent.range <- carnivore.continent.range[c((length(carnivore.continent.range)-1):(length(carnivore.continent.range)),
                                                         1:(length(carnivore.continent.range)-2))]
str(carnivore.continent.range)

# gather to go from wide to long
carnivore.continent.range_long <- gather(carnivore.continent.range, Binomial.1.2, Range.Size,3:length(carnivore.continent.range),factor_key = TRUE)
carnivore.range.pca.taxonomy <- merge(carnivore.continent.range_long, functional.axes.and.traits, by.y = "Binomial.1.2", all = TRUE)
str(carnivore.range.pca.taxonomy)

# Remove absences
continental.carnivore.trait.space <- carnivore.range.pca.taxonomy[carnivore.range.pca.taxonomy$Range.Size > 0,]
continental.carnivore.trait.space$Continent <- factor(continental.carnivore.trait.space$Continent, 
                                                      levels = c("N. America","Europe","Asia",
                                                                 "S. America","Africa","Australia"))
# Get PCoA axes too
plotly.1 <- ggplot(continental.carnivore.trait.space,
                   aes(x=A1, y = A2, colour = Family.1.2, labels = Binomial.1.2)) +
  geom_point()
ggplotly(plotly.1)



# ---- FDis: Centroid Locations ----
continental.carnivore.trait.space.test <- continental.carnivore.trait.space

continental.carnivore.trait.space.test$centroid.A1 <- continental.carnivore.trait.space.test$Range.Size*continental.carnivore.trait.space.test$A1
continental.carnivore.trait.space.test$centroid.A2 <- continental.carnivore.trait.space.test$Range.Size*continental.carnivore.trait.space.test$A2
continental.carnivore.trait.space.test$centroid.A3 <- continental.carnivore.trait.space.test$Range.Size*continental.carnivore.trait.space.test$A3
continental.carnivore.trait.space.test$centroid.A4 <- continental.carnivore.trait.space.test$Range.Size*continental.carnivore.trait.space.test$A4

centroid.location.A1 <- continental.carnivore.trait.space.test %>% group_by(Scenario,Continent) %>% summarise(centroid.A1 = sum(centroid.A1)/sum(Range.Size))
centroid.location.A2 <- continental.carnivore.trait.space.test %>% group_by(Scenario,Continent) %>% summarise(centroid.A2 = sum(centroid.A2)/sum(Range.Size))
centroid.location.A3 <- continental.carnivore.trait.space.test %>% group_by(Scenario,Continent) %>% summarise(centroid.A3 = sum(centroid.A3)/sum(Range.Size))
centroid.location.A4 <- continental.carnivore.trait.space.test %>% group_by(Scenario,Continent) %>% summarise(centroid.A4 = sum(centroid.A4)/sum(Range.Size))

centroids <- cbind(centroid.location.A1, 
                   centroid.location.A2,
                   centroid.location.A3,
                   centroid.location.A4)

centroids <- centroids[c(1:3,6,9,12)]
colnames(centroids)[3:6] <- c("A1","A2","A3","A4")



# ---- Relative geographic range loss per species ----
continental.carnivore.trait.space.test.2 <- data.frame()

for (i in 1:length(scenario)) {
  x <- continental.carnivore.trait.space.test[continental.carnivore.trait.space.test$Scenario == scenario[[i]],]
  
  for (j in 1:length(continents)) {
    y <- x[x$Continent == continents[[j]],]
    y$Relative.Range.Size <- (y$Range.Size/sum(y$Range.Size)) * 100
    continental.carnivore.trait.space.test.2 <- rbind(y,continental.carnivore.trait.space.test.2)
    print(j)
  }
}

continental.carnivore.trait.space.test <- continental.carnivore.trait.space.test.2



# ---- FRic: Minimum convex hull ----
find_hull <- function(df) df[chull(df$A1, df$A2), ]

master.hulls <- list()

for (i in 1:length(scenario)) {
  continental.hulls <- data.frame()
  hull.data <- data.frame()
  scenarionew <- scenario[[i]]
  list.number <- scenario[[i]]
  x <- continental.carnivore.trait.space.test[continental.carnivore.trait.space.test$Scenario == scenario[[i]],]
  
  
  for (j in 1:length(continents)) {
    y <- x[x$Continent == continents[[j]],]
    y.hulls <- plyr::ddply(y,"Continent", find_hull)
    hull.data <- rbind(y.hulls, hull.data)
    print(j)
  }
  
  hull.data$scenarios <- scenarionew
  master.hulls[[list.number]] <- hull.data
}

# --- Figure 2: Continental functional space ----
master.hullage <- rbind(master.hulls$Current,
                        master.hulls$`Present-Natural`)
master.hullage$Scenario <- as.factor(master.hullage$Scenario)
str(master.hullage$Scenario)
master.hullage$Scenario <- relevel(master.hullage$Scenario, "Present-Natural")

master.hullage$ContinentSpecs <- plyr::revalue(master.hullage$Continent,
                                         c("N. America" = "N. America (385 sp)",
                                           "Europe" = "Europe (125 sp)",
                                           "Asia" = "Asia (636 sp)",
                                           "S. America" = "S. America (369 sp)",
                                           "Africa" = "Africa (737 sp)",
                                           "Australia" = "Australia (134 sp)"))

# Make a new variable which has the continent names and number of species
continental.carnivore.trait.space.test$ContinentSpecs <- continental.carnivore.trait.space.test$Continent

table(continental.carnivore.trait.space.test$Continent)

continental.carnivore.trait.space.test$ContinentSpecs <- plyr::revalue(continental.carnivore.trait.space.test$ContinentSpecs,
                                               c("N. America" = "N. America (385 sp)",
                                                 "Europe" = "Europe (125 sp)",
                                                 "Asia" = "Asia (636 sp)",
                                                 "S. America" = "S. America (369 sp)",
                                                 "Africa" = "Africa (737 sp)",
                                                 "Australia" = "Australia (134 sp)"))

# Rename continents in Centroids too
centroids$ContinentSpecs <- centroids$Continent
centroids$ContinentSpecs <- plyr::revalue(centroids$ContinentSpecs,
                                                                       c("N. America" = "N. America (385 sp)",
                                                                         "Europe" = "Europe (125 sp)",
                                                                         "Asia" = "Asia (636 sp)",
                                                                         "S. America" = "S. America (369 sp)",
                                                                         "Africa" = "Africa (737 sp)",
                                                                         "Australia" = "Australia (134 sp)"))

# Save figure
pdf("./Results/Figures/Continental_Trait Space_GeographicRange.pdf", height = 4, width = 7)
  p <- ggplot() +
    geom_point(data = continental.carnivore.trait.space.test,
               aes(x=A1, y = A2,
                   fill = Scenario,
                   size = Relative.Range.Size),
               alpha = 0.15, pch = 21, colour = "black") +
    geom_density_2d(data = continental.carnivore.trait.space.test, aes(x=A1, y = A2,
                                                                colour = Scenario,
                                                                size = Relative.Range.Size),
                    size = .1, alpha = .5) +
    scale_fill_manual(values = c("#339900","#660000")) +
    scale_colour_manual(values = c("#339900","#660000")) +
    geom_point(data = centroids, aes(x = A1, y = A2, colour = Scenario), 
               size = 2, alpha = 1, shape = 3, stroke = 2) +
    geom_polygon(data = master.hullage, aes(x = A1, y = A2, colour = scenarios), 
                 alpha = 0, size = 0.6, linetype = 2) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid=element_blank(),
        strip.text = element_text(face="bold", size=9,lineheight=5.0),
        strip.background = element_rect(fill= c("black","white"), colour="black",
                                        size=1)) +
  facet_wrap(.~ContinentSpecs, ncol = 3) +
  xlab("A1 (39.5%)") + ylab("A2 (19.2%)")

g <- ggplot_gtable(ggplot_build(p))
strip_t <- which(grepl('strip-t', g$layout$name))
fills <- c("#FFCC00","#CCCC99","#FFCC99","#66CC99","#CC99CC","#99CCFF")
k <- 1
for (i in strip_t) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)
dev.off()

pdf("./Results/Figures/Updated continent trait space_UpdatedFINAL.pdf", height = 4, width = 7)
ggarrange(g,
          labels = c("a"))
dev.off()


# ---- 2. Continental without Geographic Range ----
str(carnivore.continent.presence)
carnivore.continent.presence$Scenario <- c("Present-Natural","Present-Natural","Present-Natural","Present-Natural","Present-Natural","Present-Natural",
                                           "Current","Current","Current","Current","Current","Current")
continents <- c("Africa","Asia","Australia","Europe","N. America","S. America")
carnivore.continent.presence$Continent <- c(continents,continents)
carnivore.continent.presence <- carnivore.continent.presence[c((length(carnivore.continent.presence)-1):(length(carnivore.continent.presence)),
                                                               1:(length(carnivore.continent.presence)-2))]

# gather to go from wide to long
carnivore.continent.presence.long <- gather(carnivore.continent.presence, Binomial.1.2, Range.Size,3:length(carnivore.continent.presence), factor_key = TRUE)
carnivore.range.presence.pca.taxonomy <- merge(carnivore.continent.presence.long, functional.axes.and.traits,
                                      by.y = "Binomial.1.2", all = TRUE)
str(carnivore.range.pca.taxonomy)

carnivore.range.presence.pca.taxonomy <- carnivore.range.presence.pca.taxonomy[carnivore.range.presence.pca.taxonomy$Range.Size > 0,]
carnivore.range.presence.pca.taxonomy$Continent <- factor(carnivore.range.presence.pca.taxonomy$Continent, levels = c("N. America","Europe","Asia",
                                                                                                                                "S. America","Africa","Australia"))


# ---- FDis: Centroid locations ----
centroid.location.A1 <- carnivore.range.presence.pca.taxonomy %>% group_by(Scenario,Continent) %>% summarise(centroid.A1 = sum(A1)/sum(Range.Size))
centroid.location.A2 <- carnivore.range.presence.pca.taxonomy %>% group_by(Scenario,Continent) %>% summarise(centroid.A2 = sum(A2)/sum(Range.Size))
centroid.location.A3 <- carnivore.range.presence.pca.taxonomy %>% group_by(Scenario,Continent) %>% summarise(centroid.A3 = sum(A3)/sum(Range.Size))
centroid.location.A4 <- carnivore.range.presence.pca.taxonomy %>% group_by(Scenario,Continent) %>% summarise(centroid.A4 = sum(A4)/sum(Range.Size))

centroids.No.GR <- cbind(centroid.location.A1,
                         centroid.location.A2,
                         centroid.location.A3,
                         centroid.location.A4)

centroids.No.GR <- centroids.No.GR[c(1:3, 6, 9,12)]
colnames(centroids.No.GR)[3:6] <- c("A1","A2","A3","A4")

# ---- Figure S9: Continental functional space ----
tiff("./Results/Supplementary materials/Figures/Continental_Trait Space.No.GR.tif", height = 3700, width = 6830, res = 800)
p <- ggplot(carnivore.range.presence.pca.taxonomy, aes(x=A1, y = A2, colour = Scenario, size = Range.Size)) +
  geom_point(alpha = 0.15, size = 3) +
  geom_density_2d(size = .1, alpha = .5) +
  scale_colour_manual(values = c("#339900","#660000")) +
  geom_point(data = centroids.No.GR, size = 2, alpha = 1, shape = 3, stroke = 2) +
  geom_polygon(data = master.hullage, aes(x = A1, y = A2, fill = scenarios), 
               alpha = 0, size = 0.6, linetype = 2) +
  scale_fill_discrete(guide = FALSE) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        strip.text = element_text(face="bold", size=9,lineheight=5.0),
        strip.background = element_rect(fill="white", colour="black",
                                        size=1)) +
  facet_wrap(.~Continent, ncol = 3) +
  xlab("A1 (39.5%)") + ylab("A2 (19.2%)")

g.pa <- ggplot_gtable(ggplot_build(p))
strip_t <- which(grepl('strip-t', g$layout$name))
fills <- c("#FFCC00","#CCCC99","#FFCC99","#66CC99","#CC99CC","#99CCFF")
k <- 1
for (i in strip_t) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g.pa)
dev.off()

ggarrange(g.pa)

#————————————————————————————————————————————————————————————————————————————————————————————####
# ----  7. Functional diversity changes ----
#————————————————————————————————————————————————————————————————————————————————————————————####
# First step, tidy data which I apparently have not done
Master.Global.FD.Regions <- rbind(FD.Continents.GR,
                                  FD.Continents.no.GR)

# Paired Wilcoxon signed-rank test for difference between PN and current for:
#  1. Functional richness
#  2a. Functional dispersion with geographic range
#  2b. Functional dispersion without geographic range

# --- 1. Functional richness (FRic; not influenced by weighting) ----
fric.change <- Master.Global.FD.Regions[Master.Global.FD.Regions$Method == "P/A",]
fric.change <- fric.change[fric.change$Scenario == "Current",]$FRic - fric.change[fric.change$Scenario == "Present-Natural",]$FRic
median(fric.change)

# Wilcox signed-rank test:
fric.change <- Master.Global.FD.Regions[Master.Global.FD.Regions$Method == "P/A",]

test <- data.frame(scenario = c("PN","PN","PN","PN","PN","PN", "C","C","C","C","C","C"),
           value = c(fric.change[fric.change$Scenario == "Present-Natural",]$FRic,
                     fric.change[fric.change$Scenario == "Current",]$FRic))

wilcoxsign_test(value ~ scenario,
                data = test,
                alternative = "less")

# Get the relative change in functional richness (% lost from PN)
fric.change <- Master.Global.FD.Regions[Master.Global.FD.Regions$Method == "P/A",]
fric.change.rel <- (fric.change[fric.change$Scenario == "Current",]$FRic - fric.change[fric.change$Scenario == "Present-Natural",]$FRic)/fric.change[fric.change$Scenario == "Present-Natural",]$FRic
median(fric.change.rel)
quantile(fric.change.rel)


# ---- 2a. Functional dispersion with geographic range ----
fdis.change.gr <- Master.Global.FD.Regions[Master.Global.FD.Regions$Method == "Geographic Range",]
fdis.change.gr <- fdis.change.gr[fdis.change.gr$Scenario == "Current",]$FDis - fdis.change.gr[fdis.change.gr$Scenario == "Present-Natural",]$FDis
median(fdis.change.gr)

# Wilcox signed-rank test:
fdis.change.gr <- Master.Global.FD.Regions[Master.Global.FD.Regions$Method == "Geographic Range",]

test2 <- data.frame(scenario = c("PN","PN","PN","PN","PN","PN", "C","C","C","C","C","C"),
                   value = c(fdis.change.gr[fdis.change.gr$Scenario == "Present-Natural",]$FDis,
                             fdis.change.gr[fdis.change.gr$Scenario == "Current",]$FDis))
wilcoxsign_test(value ~ scenario, data = test2, alternative = "less")

# Relative change
fdis.change.gr <- Master.Global.FD.Regions[Master.Global.FD.Regions$Method == "Geographic Range",]
fdis.change.gr <- (fdis.change.gr[fdis.change.gr$Scenario == "Current",]$FDis - fdis.change.gr[fdis.change.gr$Scenario == "Present-Natural",]$FDis)/fdis.change.gr[fdis.change.gr$Scenario == "Present-Natural",]$FDis
median(fdis.change.gr)
quantile(fdis.change.gr)


# ---- 2b. Functional dispersion without geographic range ----
fdis.change.pa <- Master.Global.FD.Regions[Master.Global.FD.Regions$Method == "P/A",]
fdis.change.pa <- fdis.change.pa[fdis.change.pa$Scenario == "Current",]$FDis - fdis.change.pa[fdis.change.pa$Scenario == "Present-Natural",]$FDis

# Wilcox signed-rank test:
fdis.change.pa <- Master.Global.FD.Regions[Master.Global.FD.Regions$Method == "P/A",]

test3 <- data.frame(scenario = c("PN","PN","PN","PN","PN","PN", "C","C","C","C","C","C"),
                   value = c(fdis.change.pa[fdis.change.pa$Scenario == "Present-Natural",]$FDis,
                             fdis.change.pa[fdis.change.pa$Scenario == "Current",]$FDis))
wilcoxsign_test(value ~ scenario, data = test3, alternative = "less")


fdis.change.pa <- Master.Global.FD.Regions[Master.Global.FD.Regions$Method == "P/A",]
fdis.change.pa <- (fdis.change.pa[fdis.change.pa$Scenario == "Current",]$FDis - fdis.change.pa[fdis.change.pa$Scenario == "Present-Natural",]$FDis)/fdis.change.pa[fdis.change.pa$Scenario == "Present-Natural",]$FDis
median(fdis.change.pa)
quantile(fdis.change.pa)


# ---- Figure 2 (b,c):  Functional diversity differences ----
FRic.Change.Continents <- Master.Global.FD.Regions[Master.Global.FD.Regions$Region == "Continent",]
FRic.Change.Continents.GR <- Master.Global.FD.Regions[Master.Global.FD.Regions$Method == "Geographic Range",]
FRic.Change.Continents.PA <- Master.Global.FD.Regions[Master.Global.FD.Regions$Method == "P/A",]

FRic.Change.Continents.GR$Region.names <- factor(FRic.Change.Continents.GR$Region.names, 
                                                 levels = c("N. America","S. America",
                                                            "Europe","Africa",
                                                            "Asia","Australia"))

good.colours <- c("#FFCC99","#99CCFF","#CCCC99","#CC99CC","#FFCC00","#66CC99")
good.colours <- c("#66CC99","#FFCC00","#CC99CC","#CCCC99","#99CCFF","#FFCC99")


levels(FRic.Change.Continents.GR$Scenario)[levels(FRic.Change.Continents.GR$Scenario) == "Present-Natural"] <- "PN"
levels(FRic.Change.Continents.GR$Scenario)[levels(FRic.Change.Continents.GR$Scenario) == "Current"] <- "C"

# Functional richness change
(p1 <- ggplot(FRic.Change.Continents.GR,
              aes(y = FRic, x = Scenario, fill = Region.names, group = Scenario)) +
    geom_line(aes(group = Region.names), colour = "black", linetype = "dashed") +
    geom_point(size = 4, aes(fill = Region.names), pch = 21, colour = "black") +
  scale_fill_manual(values = good.colours, name = "Continent") +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()))

FDis.Change.Continents <- Master.Global.FD.Regions[Master.Global.FD.Regions$Region == "Continent",]
FDis.Change.Continents$Scenario <- factor(FDis.Change.Continents$Scenario, levels = c("Present-Natural","Current"))

levels(FDis.Change.Continents$Scenario)[levels(FDis.Change.Continents$Scenario) == "Present-Natural"] <- "PN"
levels(FDis.Change.Continents$Scenario)[levels(FDis.Change.Continents$Scenario) == "Current"] <- "C"

FDis.Change.Continents$Region.names <- factor(FDis.Change.Continents$Region.names, 
                                              levels = c("N. America","S. America",
                                                         "Europe","Africa",
                                                         "Asia","Australia"))

# Functional dispersion
(p2 <- ggplot(FDis.Change.Continents,
              aes(y = FDis, x = Scenario, fill = Region.names, group = Scenario)) +
    geom_line(aes(group = Region.names), colour = "black", linetype = "dashed") +
    geom_point(size = 4, aes(fill = Region.names, shape = Method), colour = "black") +
  scale_fill_manual(values = good.colours, name = "Continent") +
    scale_shape_manual(values = c(21,24)) +
  facet_wrap(. ~ Method) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank()))

# Combine functional richness and functional dispersion
pdf("./Results/Figures/Continental.FRic.FDis.Change_UpdatedFINAL.pdf", height = 2.6, width = 7)
ggarrange(p1 + theme(legend.position="none",
                     axis.title.x = element_blank()),
          # axis.title = element_text(size = 12),
          #axis.text = element_text(size = 10)),
          p2 + theme(axis.title.x = element_blank(),
                     #axis.title = element_text(size = 12),
                     #axis.text = element_text(size = 10),
                     legend.position = "none"),
          labels = c("b","c"),
          ncol = 2, nrow = 1)
dev.off()


#————————————————————————————————————————————————————————————————————————————————————————————####
#  ---- 8. Centroid magnitude shift from PN to current ----
#————————————————————————————————————————————————————————————————————————————————————————————####

# ---- Centroid shift with geographic range ----
str(centroids)
continents <- c("N. America","Europe","Asia","S. America","Africa","Australia")

centroids$ContinentSpecs <- NULL

centroids.2 <- cbind(centroids[1:6,], centroids[7:12,])
centroids.2 <- centroids.2[c(1:6,9:12)]
centroids.2$Change.A1 <- centroids.2$A1-centroids.2$A11
centroids.2$Change.A2 <- centroids.2$A2-centroids.2$A21
centroids.2$Change.A3 <- centroids.2$A3-centroids.2$A31
centroids.2$Change.A4 <- centroids.2$A4-centroids.2$A41

centroid.position.change <- centroids.2[c(2,11:14)]
centroid.position.change$Method <- "Geographical Range"

centroid.change.continent.GR <- as.data.frame(cbind(centroids.2$Continent, centroids.2$Change.A1))
colnames(centroid.change.continent.GR) <- c("Continent","Change")
centroid.change.continent.GR$Continent <- continents
centroid.change.continent.GR$A <- "1"

centroid.change.continent.GR.2 <- as.data.frame(cbind(centroids.2$Continent, centroids.2$Change.A2))
colnames(centroid.change.continent.GR.2) <- c("Continent","Change")
centroid.change.continent.GR.2$Continent <- continents
centroid.change.continent.GR.2$A <- "2"

centroid.change.continent.GR.3 <- as.data.frame(cbind(centroids.2$Continent, centroids.2$Change.A3))
colnames(centroid.change.continent.GR.3) <- c("Continent","Change")
centroid.change.continent.GR.3$Continent <- continents
centroid.change.continent.GR.3$A <- "3"

centroid.change.continent.GR.4 <- as.data.frame(cbind(centroids.2$Continent, centroids.2$Change.A4))
colnames(centroid.change.continent.GR.4) <- c("Continent","Change")
centroid.change.continent.GR.4$Continent <- continents
centroid.change.continent.GR.4$A <- "4"

centroid.change.continent.GR <- rbind(centroid.change.continent.GR, 
                                      centroid.change.continent.GR.2,
                                      centroid.change.continent.GR.3,
                                      centroid.change.continent.GR.4)

centroid.change.continent.GR$Method <- "Geographic range"

# ---- b. Centroid shift without geographic range ----
str(centroids)

centroids.2 <- cbind(centroids.No.GR[1:6,], centroids.No.GR[7:12,])
centroids.2 <- centroids.2[c(1:6,9:12)]
centroids.2$Change.A1 <- centroids.2$A1-centroids.2$A11
centroids.2$Change.A2 <- centroids.2$A2-centroids.2$A21
centroids.2$Change.A3 <- centroids.2$A3-centroids.2$A31
centroids.2$Change.A4 <- centroids.2$A4-centroids.2$A41

centroid.position.change.2 <- centroids.2[c(2,11:14)]
centroid.position.change.2$Method <- "Presence/Absence"
Master.centroid.position.change <- rbind(centroid.position.change,
                                         centroid.position.change.2)

centroid.change.continent.no.GR <- as.data.frame(cbind(centroids.2$Continent, centroids.2$Change.A1))
colnames(centroid.change.continent.no.GR) <- c("Continent","Change")
centroid.change.continent.no.GR$Continent <- continents
centroid.change.continent.no.GR$A <- "1"

centroid.change.continent.no.GR.2 <- as.data.frame(cbind(centroids.2$Continent, centroids.2$Change.A2))
colnames(centroid.change.continent.no.GR.2) <- c("Continent","Change")
centroid.change.continent.no.GR.2$Continent <- continents
centroid.change.continent.no.GR.2$A <- "2"

centroid.change.continent.no.GR.3 <- as.data.frame(cbind(centroids.2$Continent,centroids.2$Change.A3))
colnames(centroid.change.continent.no.GR.3) <- c("Continent","Change")
centroid.change.continent.no.GR.3$Continent <- continents
centroid.change.continent.no.GR.3$A <- "3"

centroid.change.continent.no.GR.4 <- as.data.frame(cbind(centroids.2$Continent,centroids.2$Change.A4))
colnames(centroid.change.continent.no.GR.4) <- c("Continent","Change")
centroid.change.continent.no.GR.4$Continent <- continents
centroid.change.continent.no.GR.4$A <- "4"

centroid.change.continent.no.GR <- rbind(centroid.change.continent.no.GR,
                                         centroid.change.continent.no.GR.2,
                                         centroid.change.continent.no.GR.3,
                                         centroid.change.continent.no.GR.4)

centroid.change.continent.no.GR$Method <- "Presence/Absence"

Master.Centroid.Shift <- rbind(centroid.change.continent.GR, centroid.change.continent.no.GR)

Master.Centroid.Shift$Continent <- factor(Master.Centroid.Shift$Continent, levels = c("N. America","S. America","Europe","Africa","Asia","Australia"))
Master.Centroid.Shift$Method <- factor(Master.Centroid.Shift$Method, levels = c("Presence/Absence","Geographic range"))

Master.centroid.position.change$Method <- factor(Master.centroid.position.change$Method, levels = c("Presence/Absence","Geographical Range"))


# ---- Statistics ---- 

# 1. KW Test: change between principle coordinate axes (pool methods)

# Unweighted analyses (without geographic range)
Master.Centroid.Shift.PA <- Master.Centroid.Shift[Master.Centroid.Shift$Method == "Presence/Absence",]
Master.Centroid.Shift.PA$A <- as.factor(Master.Centroid.Shift.PA$A)
kruskal.test(Master.Centroid.Shift.PA$Change, Master.Centroid.Shift.PA$A)
dunn.test(Master.Centroid.Shift.PA$Change,Master.Centroid.Shift.PA$A)

# Weighted analyses (with geographic range)
Master.Centroid.Shift.GR <- Master.Centroid.Shift[Master.Centroid.Shift$Method == "Geographic range",]
Master.Centroid.Shift.GR$A <- as.factor(Master.Centroid.Shift.GR$A)
kruskal.test(Master.Centroid.Shift.GR$Change, Master.Centroid.Shift.GR$A)
dunn.test(Master.Centroid.Shift.GR$Change,Master.Centroid.Shift.GR$A)


# 2. MW Tests: changes within principle components - comparing methods 

# PCoA 1
wilcox.test(Master.centroid.position.change$Change.A1 ~ Master.centroid.position.change$Method)
wilcox_test(Master.centroid.position.change$Change.A1 ~ Master.centroid.position.change$Method)
tapply(Master.centroid.position.change$Change.A1, Master.centroid.position.change$Method, median)

# Find the U-value
# U = n1n2  +  n2(n2+1)/2  – R2 
(u1 <- (6*6) + (6*(6+1))/2 - sum(rank(Master.centroid.position.change$Change.A1)[1:6]))
(u2 <- (6*6) + (6*(6+1))/2 - sum(rank(Master.centroid.position.change$Change.A1)[7:12]))


# PCoA 2
wilcox.test(Master.centroid.position.change$Change.A2 ~ Master.centroid.position.change$Method)
wilcox_test(Master.centroid.position.change$Change.A2 ~ Master.centroid.position.change$Method)
tapply(Master.centroid.position.change$Change.A2, Master.centroid.position.change$Method, median)

(u1 <- (6*6) + (6*(6+1))/2 - sum(rank(Master.centroid.position.change$Change.A2)[1:6]))
(u2 <- (6*6) + (6*(6+1))/2 - sum(rank(Master.centroid.position.change$Change.A2)[7:12]))

# PCoA 3
wilcox.test(Master.centroid.position.change$Change.A3 ~ Master.centroid.position.change$Method)
wilcox_test(Master.centroid.position.change$Change.A3 ~ Master.centroid.position.change$Method)
tapply(Master.centroid.position.change$Change.A3, Master.centroid.position.change$Method, median)

(u1 <- (6*6) + (6*(6+1))/2 - sum(rank(Master.centroid.position.change$Change.A3)[1:6]))
(u2 <- (6*6) + (6*(6+1))/2 - sum(rank(Master.centroid.position.change$Change.A3)[7:12]))


# PCoA 4
wilcox.test(Master.centroid.position.change$Change.A4 ~ Master.centroid.position.change$Method)
wilcox_test(Master.centroid.position.change$Change.A4 ~ Master.centroid.position.change$Method)
tapply(Master.centroid.position.change$Change.A4, Master.centroid.position.change$Method, median)

(u1 <- (6*6) + (6*(6+1))/2 - sum(rank(Master.centroid.position.change$Change.A4)[1:6]))
(u2 <- (6*6) + (6*(6+1))/2 - sum(rank(Master.centroid.position.change$Change.A4)[7:12]))

# ---- Figure 3: Magnitude of shifts in centroid ---- 

# Tidying
Master.Centroid.Shift$A <- as.factor(Master.Centroid.Shift$A)
Master.Centroid.Shift$A <- plyr::revalue(Master.Centroid.Shift$A, c("1" = "1 (39.5%)",
                                                                    "2" = "2 (19.2%)",
                                                                    "3" =  "3 (8.5%)",
                                                                    "4" =  "4 (5.8%)"))
Master.Centroid.Shift$Method <- factor(reorder(Master.Centroid.Shift$Method, levels = c("Geographic range","Presence/Absence")))


# Main plot
(p2 <- ggplot(Master.Centroid.Shift, aes(x=reorder(A, desc(A)), y = Change, group = Method, shape = Method)) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.5,
               alpha = 0.5, colour = "black") +
    geom_point(aes(fill = Continent),
               position = position_jitterdodge(dodge.width = 0.9, jitter.width = 0.2),
               size = 3) +
    scale_shape_manual(values = c(21,24)) +
    scale_fill_manual(values = good.colours) +
    geom_pointrange(mapping = aes(x = reorder(A, desc(A)), y = Change),
                    stat = "summary",
                    fun.ymin = function(z) {quantile(z,0.25)},
                    fun.ymax = function(z) {quantile(z,0.75)},
                    fun.y = median, 
                    position = position_jitterdodge(dodge.width = 0.3, jitter.width = 0.2),
                    colour = "black",
                    fill = "black",
                    size = .7) +
    theme_bw() +
    coord_flip() +
    theme(axis.text = element_text(size = 15),
          axis.text.y = element_text(size = 15), 
          axis.title = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          panel.grid = element_blank()) +
    geom_vline(xintercept = c(1.5,2.5,3.5), colour = "black") +
    xlab("Principal coordinate axis") + ylab("Centroid shift") +
    ylim(-0.15,0.08))

#Save without the centroid positional plot
pdf("./Results/Figures/Figure 3.pdf", height = 4.5, width = 7.5)
p2
dev.off()

#————————————————————————————————————————————————————————————————————————————————————————————####
#  ----- Sensitivity analyses -----
#————————————————————————————————————————————————————————————————————————————————————————————####

# ---- 1. BPMM on systematically reviewed species ----

# Identify the systematic species
systematic.species <- spec[!is.na(spec$systematic.species),]

# Prune tree to only the species list being used in this study that have been systematically checked
pruned.forest <- lapply(tree_samp, 
                        function(x) {
                          drop.tip(x, x$tip.label[!x$tip.label %in% levels(factor(systematic.species$Binomial.1.2))])
                        }
)
class(pruned.forest) <- "multiPhylo"

# Check a single tree
plot(pruned.forest[[1]])

# Number of tips in both
unlist(lapply(pruned.forest, Ntip)) # Should be 555 as those are species which have been systematically checked.

# The multinomial2 model requires predictor variables as integers
is.integer(systematic.species$CurrentRange) # TRUE
is.integer(systematic.species$Rangeloss) # FALSE - not sure why but change to be integer...

# Running the MCMCglmm on multiple trees - need to save it as mulTree object
# Creating the mulTree object - originally try on a subset of 2 trees
mulTree_data <- as.mulTree(data = systematic.species,
                           tree = pruned.forest[1:100],
                           taxa = "Binomial.1.2")
str(mulTree_data)
show(mulTree_data)

# The MCMC parameters (model iterations, thinning interval, burnin)
my_parameters <- c(200000, # Number of mcmcglmm iterations,
                   100,    # Thinning interval, default = 10
                   10000)  # Burnin period - discard first 3000

# The MCMCglmm priors
# Uninformative priors for:
# (R) the model residuals
# (G) the random effects [phylogenetic covariance matrix]
my_priors <- list(R = list(V = 1/2, nu = 0.002),# Covariance matrix of model residuals
                  G = list(G1 = list(V = 1/2, nu = 0.002))) # Random effects priors: V = scale matrix, nu = df

# This object is classified as "mulTree" and contains different elements to be
# passed to the mulTree function
class(mulTree_data) ; names(mulTree_data)


# ---- Hierarchical model running ----
# Create list to store the DIC outputs
DIC.list <- list()

# 1. Global model
mulTree(cbind(Rangeloss, CurrentRange) ~ scaled.mass * scaled.diet, 
        mulTree.data = mulTree_data,
        priors = my_priors,
        parameters = my_parameters, 
        output = "./Results/MCMCglmm chain outputs/systematic species/Model1/range_loss", # Added on an extra number to make sure it can be separated from the others
        family = "multinomial2",
        ESS = 1000, chains = 2) # Using the default assessment of effective sample size

# Get average DIC value and interquartile range from all of the models
DICs <- matrix(ncol = 1, nrow = length(pruned.forest[1:100]))
setwd("./Results/MCMCglmm chain outputs/systematic species/Model1")
for(i in 1:length(pruned.forest[1:100])) {
  x <- read.mulTree(paste0("range_loss-tree",i,"_chain1"),
                    model = TRUE)
  DICs[i,1] <- x$DIC
  print(i)
}
DIC.list[[1]] <- DICs

# re-set working directory
setwd("../../../../")


# 2. No interaction model
mulTree(cbind(Rangeloss, CurrentRange) ~ scaled.mass + scaled.diet, 
        mulTree.data = mulTree_data,
        priors = my_priors,
        parameters = my_parameters, 
        output = "./Results/MCMCglmm chain outputs/systematic species/Model2/range_loss",
        family = "multinomial2",
        ESS = 1000, chains = 2)

# Get average DIC value and interquartile range from all of the models
DICs <- matrix(ncol = 1, nrow = length(pruned.forest[1:100]))
setwd("./Results/MCMCglmm chain outputs/systematic species/Model2")
for(i in 1:length(pruned.forest[1:100])) {
  x <- read.mulTree(paste0("range_loss-tree",i,"_chain1"),
                    model = TRUE)
  DICs[i,1] <- x$DIC
  print(i)
}
DIC.list[[2]] <- DICs

# Re-set working directory
setwd("../../../../")


# 3. Only mass
mulTree(cbind(Rangeloss, CurrentRange) ~ scaled.mass, 
        mulTree.data = mulTree_data,
        priors = my_priors,
        parameters = my_parameters, 
        output = "./Results/MCMCglmm chain outputs/systematic species/Model3/range_loss",
        family = "multinomial2",
        ESS = 1000, chains = 2)

# Get average DIC value and interquartile range from all of the models
DICs <- matrix(ncol = 1, nrow = length(pruned.forest[1:100]))
setwd("./Results/MCMCglmm chain outputs/systematic species/Model3")
for(i in 1:length(pruned.forest[1:100])) {
  x <- read.mulTree(paste0("range_loss-tree",i,"_chain1"),
                    model = TRUE)
  DICs[i,1] <- x$DIC
  print(i)
}
DIC.list[[3]] <- DICs

# Re-set working directory
setwd("../../../../")


# 4. Only diet
mulTree(cbind(Rangeloss, CurrentRange) ~ scaled.diet, 
        mulTree.data = mulTree_data,
        priors = my_priors,
        parameters = my_parameters, 
        output = "./Results/MCMCglmm chain outputs/systematic species/Model4/range_loss",
        family = "multinomial2",
        ESS = 1000, chains = 2)

# Get average DIC value and interquartile range from all of the models
DICs <- matrix(ncol = 1, nrow = length(pruned.forest[1:100]))
setwd("./Results/MCMCglmm chain outputs/systematic species/Model4")
for(i in 1:length(pruned.forest[1:100])) {
  x <- read.mulTree(paste0("range_loss-tree",i,"_chain1"),
                    model = TRUE)
  DICs[i,1] <- x$DIC
  print(i)
}

DIC.list[[4]] <- DICs

# Identify which of the models is the best - lowest average DIC values
# The DIC value can be interpreted in a similar way to an AIC value
lapply(DIC.list, quantile)
lapply(DIC.list,mean)

# ---- Diagnostic plots ----
setwd("../../../../")

# Model 1
setwd("./Results/MCMCglmm chain outputs/systematic species/Model1")
chain.1 <- read.mulTree(paste0("range_loss-tree","20","_chain1"),
                        model = TRUE)
chain.2 <- read.mulTree(paste0("range_loss-tree","20","_chain2"),
                        model = TRUE)

### Diagnostic check for how well the chains converged on the same result for posterior
### distributions of the fixed effects
sol <- bind_rows(as.data.frame(chain.1$Sol[, 1:4]), 
                 as.data.frame(chain.2$Sol[, 1:4]))

sol["chain"] <- gl(2, 1900) # Number of chains, number of iterations(?)
sol["sample"] <- rep(1:1900, 2)
sol <- gather(sol, key = "variable", value = "value", -chain, -sample)

left <- ggplot(sol, aes(x = sample, y = value, col = chain)) + # Should look like a fuzzy caterpillar
  geom_line() + 
  facet_wrap(~ variable, scales = "free", nrow = 4) +
  theme(legend.position="none") + 
  ylab("")
right <- ggplot(sol, aes(x = value, col = chain)) + # Distribution of fixed effects parameter estimates
  geom_density() +
  geom_rug() +
  facet_wrap(~ variable, scales = "free", nrow = 4) +
  theme(legend.position="none") + 
  labs(x = "", y = "")
grid.arrange(left, right, nrow = 1)

# Need to also get the distribution of the random effects parameters
sol <- bind_rows(as.data.frame(chain.1$VCV[, 1:2]), 
                 as.data.frame(chain.2$VCV[, 1:2]))

sol["chain"] <- gl(2, 1900) # Number of chains, number of iterations(?)
sol["sample"] <- rep(1:1900, 2)
sol <- gather(sol, key = "variable", value = "value", -chain, -sample)

left <- ggplot(sol, aes(x = sample, y = value, col = chain)) + # Should look like a fuzzy caterpillar
  geom_line() + 
  facet_wrap(~ variable, scales = "free", nrow = 4) +
  theme(legend.position="none") + 
  ylab("")
right <- ggplot(sol, aes(x = value, col = chain)) + # Distribution of fixed effects parameter estimates
  geom_density() +
  geom_rug() +
  facet_wrap(~ variable, scales = "free", nrow = 4) +
  theme(legend.position="none") + 
  labs(x = "", y = "")
grid.arrange(left, right, nrow = 1)

# Checking convergence for our fixed factors (close to one is good)
class(chain.2$Sol)
gelman.diag(mcmc.list(as.mcmc(chain.1$Sol[, 1:4]),
                      as.mcmc(chain.2$Sol[, 1:4])),
            autoburnin = FALSE)

# Checking convergence for our random terms
gelman.diag(mcmc.list(chain.1$VCV,
                      chain.2$VCV),
            autoburnin = FALSE)

### Checking only the first chain
effectiveSize(chain.1$Sol[, 1:4]) 
effectiveSize(chain.1$VCV) 

# acf plot for the fixed estimates
acf(chain.1$Sol[, 1], lag.max = 20)
acf(chain.1$Sol[, 2], lag.max = 20)
acf(chain.1$Sol[, 3], lag.max = 20)
acf(chain.1$Sol[, 4], lag.max = 20)

# acf plot for the first random term in our model (the phyl term)
acf(chain.1$VCV[, 1], lag.max = 20)
acf(chain.1$VCV[, 2], lag.max = 20)


# ---- Model predictions on raw data ----

# Model 1 and 2 were equal models: keep interaction between body mass and
# vertebrate consumption.

# Diagnoses and summary of the best-supported model.
# Convergence diagnosis test to see if the two chains converged foreach tree 
read.mulTree("range_loss", convergence = TRUE) ## As indicated here, the chains converged for both chains!

# Reading all the models to perform the MCMCglmm analysis on multiple trees
all_models <- read.mulTree("range_loss") 
str(all_models) ## Contains 39,600 estimations of the Intercept and the terms!
class(all_models)

# Summarising the results by estimating the highest density regions 
# and their associated 95 and 50 confidence intervals (default) 
summarised_results <- summary(all_models, use.hdr = FALSE)
summarised_results

# Create new dataframe for predicting
colnames(systematic.species)
new.df <- systematic.species[c(1,7,8,10:11,13)]
colnames(new.df)
colnames(new.df)[1] <- "animal"
colnames(new.df)

# Create matrix to store model predictions:
model.predicts <- matrix(nrow = nrow(systematic.species), ncol = length(tree_samp))
model.predicts <- list()

# Get model estimates from each of the models
for(i in 1:length(tree_samp)) {
  one_model <- read.mulTree(paste0("range_loss-tree",i,"_chain1"),
                            model = TRUE) ## This model is a normal MCMCglmm object that has been ran on one single tree class(one_model) ; names(one_model)
  
  # Get model estimates from each of the models
  test <- predict(one_model, new.df,
                  type = "response", interval = "confidence")
  model.predicts[[i]] <- test
  print(i)
}

predict.matrix <- matrix(nrow = nrow(systematic.species), ncol = 100)
lower.matrix <- matrix(nrow = nrow(systematic.species), ncol = 100)
upper.matrix <- matrix(nrow = nrow(systematic.species), ncol = 100)

# Get average for the 100 models
for(i in 1:100) {
  x <- model.predicts[[i]]
  predict.matrix[,i] <- x[,1]
  lower.matrix[,i] <- x[,2]
  upper.matrix[,i] <- x[,3]
  print(i)
}

systematic.species$PredictedLoss <- apply(predict.matrix,1,mean)
systematic.species$PredictedLossLower <- apply(lower.matrix,1,mean)
systematic.species$PredictedLossUpper <- apply(upper.matrix,1,mean)
systematic.species$ProportionalPredictedLoss <- systematic.species$PredictedLoss/systematic.species$PresnatRange
systematic.species$ProportionalPredictedLossLower <- systematic.species$PredictedLossLower/systematic.species$PresnatRange
systematic.species$ProportionalPredictedLossUpper <- systematic.species$PredictedLossUpper/systematic.species$PresnatRange

# Add prediction to species dataframe
systematic.species$pred.rangeloss <- test
systematic.species$pred.proploss <- systematic.species$pred.rangeloss/systematic.species$PresnatRange

# ggplot of model predictions on raw data
# Create new factor of vertebrate consumption categories
systematic.species$size.cat <- NA
systematic.species[systematic.species$Diet.Vertebrate >= 70,]$size.cat <- "1"
systematic.species[systematic.species$Diet.Vertebrate < 70,]$size.cat <- "2"

# Creating two figures separately:
# 
# 1. Predictions and raw data for hypercarnivores (>=70%)
# 2. Predictions and raw data for non-hypercarnivores (5-69%)

# Note: Need to do a hack around getting the shaded ribbon lines around the estimate

# 1. Hypercarnivores (>=70%)
g1 <- ggplot(data = systematic.species[systematic.species$size.cat == "1",]) +
  geom_point(aes(x = log.bodymass, y = Prop_RangeContraction*-1),
             alpha = 0.3) +
  geom_smooth(aes(x = log.bodymass, y = ProportionalPredictedLossLower*-1),
              method = "loess", se = FALSE, linetype = 0, size = 0.5, colour = "black") +
  geom_smooth(aes(x = log.bodymass, y = ProportionalPredictedLossUpper*-1),
              method = "loess", se = FALSE, linetype = 0, alpha = 0.1, size = 0.5, colour = "black") +
  #scale_colour_brewer(palette = "Paired") + 
  scale_x_continuous(breaks = c(0:6),
                     labels = c("0","0.01", "0.1", "1", "10", "100", "1000"),
                     limits = c(0,6)) +
  xlab("Body mass (kg)") +
  ylab("Proportional range difference") +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
g1

# build plot object for rendering 
gg1 <- ggplot_build(g1)

# extract data for the loess lines from the 'data' slot
df2 <- data.frame(x = gg1$data[[2]]$x,
                  ymin = gg1$data[[2]]$y,
                  ymax = gg1$data[[3]]$y) 

# use the loess data to add the 'ribbon' to plot 
gg1 <- g1 +
  geom_smooth(aes(x = log.bodymass, y = ProportionalPredictedLoss*-1),
              method = "loess", se = FALSE, colour = "black") +
  geom_ribbon(data = df2, aes(x = x, ymin = ymin, ymax = ymax),
              fill = "grey", alpha = 0.4, linetype = 0)
gg1


# 2. Non-hypercarnivores (5-69%)
g2 <- ggplot(data = systematic.species[systematic.species$size.cat == "2",]) +
  geom_point(aes(x = log.bodymass, y = Prop_RangeContraction*-1),
             alpha = 0.3) +
  geom_smooth(aes(x = log.bodymass, y = ProportionalPredictedLossLower*-1),
              method = "loess", se = FALSE, linetype = 0, size = 0.5, colour = "black") +
  geom_smooth(aes(x = log.bodymass, y = ProportionalPredictedLossUpper*-1),
              method = "loess", se = FALSE, linetype = 0, alpha = 0.1, size = 0.5, colour = "black") +
  scale_x_continuous(breaks = c(0:6),
                     labels = c("0","0.01", "0.1", "1", "10", "100", "1000"),
                     limits = c(0,6)) +
  xlab("Body mass (kg)") +
  ylab("Proportional range difference") +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
g2

# build plot object for rendering 
gg2 <- ggplot_build(g2)

# extract data for the loess lines from the 'data' slot
df2 <- data.frame(x = gg2$data[[2]]$x,
                  ymin = gg2$data[[2]]$y,
                  ymax = gg2$data[[3]]$y) 

# use the loess data to add the 'ribbon' to plot 
gg2 <- g2 +
  geom_smooth(aes(x = log.bodymass, y = ProportionalPredictedLoss*-1),
              method = "loess", se = FALSE, colour = "black") +
  geom_ribbon(data = df2, aes(x = x, ymin = ymin, ymax = ymax),
              fill = "grey", alpha = 0.4, linetype = 0)


gg1
gg2

tiff("../../../Figures/BPMM_Results systematic species.tif",
     height = 1500, width = 3300, res = 400)
plot_grid(gg2, gg1,
          align = "v", labels = "auto")
dev.off()

pdf("../../../Figures/BPMM_Results systematic species.pdf",
    height = 4, width = 7)
plot_grid(gg2, gg1,
          align = "v", labels = "auto")
dev.off()

#————————————————————————————————————————————————————————————————————————————————————————————####
# ----- 2. Functional diversity with systematically-checked species ----

# Create a dataframe with just those species that were systematically checked:
systematic.species <- species2 %>% filter(systematic.species == 1)

# Remove species that are not in the systematically checked species list:
keep <- which(!colnames(current.output[c(8:1088)]) %in% levels(factor(systematic.species$Binomial.1.2)))
current.output <- cbind(current.output[c(1:7)], current.output[c(8:1088)][-keep])

keep <- which(!colnames(present.natural.output[c(8:1088)]) %in% levels(factor(systematic.species$Binomial.1.2)))
present.natural.output <- cbind(present.natural.output[c(1:7)], present.natural.output[c(8:1088)][-keep])

# Calculate carnivore continental range
# cbind current and present-natural outputs to the geographical information and coordinates, then rbind both scenarios
global.scenarios.merged <- rbind(current.output,present.natural.output)
global.scenarios.merged$Scenario <- as.factor(global.scenarios.merged$Scenario)
global.scenarios.merged$Continent <- as.factor(global.scenarios.merged$Continent)
global.scenarios.merged$realm <- as.factor(global.scenarios.merged$realm)
str(global.scenarios.merged)

# Sum all species per continent 
continents <- levels(global.scenarios.merged$Continent)
scenario <- levels(global.scenarios.merged$Scenario)

carnivore.continent.range <- data.frame()

for (i in 1:length(scenario)) {
  x <- global.scenarios.merged[global.scenarios.merged$Scenario == scenario[[i]],]
  continent.range <- matrix(ncol = length(colnames(current.output[-c(1:7)])), nrow = 6)
  
  for (j in 1:length(continents)) {
    x <- x[!is.na(x$Continent),]
    y <- x[x$Continent == continents[[j]],]
    y <- y[c(8:length(current.output))]
    z <- colSums(y)
    continent.range[j,] <- z
    print(j)
  }
  
  continent.range <- as.data.frame(continent.range)
  carnivore.continent.range <- rbind(continent.range, carnivore.continent.range)
  print(i)
}

colnames(carnivore.continent.range) <- colnames(current.output[c(8:length(current.output))])
carnivore.continent.range$Scenario <- c("Present-Natural","Present-Natural","Present-Natural","Present-Natural","Present-Natural","Present-Natural",
                                        "Current","Current","Current","Current","Current","Current")
carnivore.continent.range$Continent <- c(continents,continents)
carnivore.continent.range <- carnivore.continent.range[c((length(carnivore.continent.range)-1):(length(carnivore.continent.range)),
                                                         1:(length(carnivore.continent.range)-2))]
str(carnivore.continent.range)

#  ---- a) Functional diversity with range size ----

# 1. Community Matrix
row.names(carnivore.continent.range) <- paste(carnivore.continent.range$Scenario, carnivore.continent.range$Continent, sep = "_")
carnivore.continent.range <- carnivore.continent.range[-c(1:2)]
carnivore.not.on.continents <- subset(carnivore.continent.range, select = names(which(colSums(carnivore.continent.range) == 0))) # which species are we losing here
carnivore.not.on.continents <- as.data.frame(colnames(carnivore.not.on.continents))
colnames(carnivore.not.on.continents) <- "Binomial.1.2"
carnivore.not.on.continents <- merge(carnivore.not.on.continents,species,by = "Binomial.1.2") # Species being lost due to ranges not falling into continental shapefiles...

carnivore.continent.range <- subset(carnivore.continent.range, select = names(which(colSums(carnivore.continent.range) != 0))) 
str(carnivore.continent.range)

# 2. Traits
new.species <- as.data.frame(colnames(carnivore.continent.range))
colnames(new.species)[1] <- "Binomial.1.2"
traits.sorted <- merge(species, new.species)
row.names(traits.sorted) <- traits.sorted$Binomial.1.2

str(traits.sorted)
row.names(traits.sorted) <- traits.sorted$Binomial.1.2
traits.sorted <- traits.sorted[c(10,12:14)]

# Transform body mass 
traits.sorted$Mass.g <- log10(traits.sorted$Mass.g)

# Scale continuous traits
traits.sorted <- scale(traits.sorted)

# Put matrix and traits into a list to call dbFD
joined.world.function <- list(traits.sorted, carnivore.continent.range)
names(joined.world.function)[1] <- "traits"
names(joined.world.function)[2] <- "presence"
str(joined.world.function$traits)
str(joined.world.function$presence)

# Set weights
w <- c(1,1/3,1/3,1/3)

# Call dbFD 
global.FD.continents.GR <- dbFD(joined.world.function$traits,
                                joined.world.function$presence,
                                w,
                                m = 4,
                                stand.FRic = TRUE,
                                calc.CWM = TRUE,
                                print.pco = TRUE) # 0.74 R^2 value

#     Extract and look at results of differences between scenarios at the level of continents
FD.Continents.GR <- as.data.frame(cbind(global.FD.continents.GR$nbsp,
                                        global.FD.continents.GR$FRic,
                                        global.FD.continents.GR$FDis))

FD.Continents.GR$Scenario <- c("Present-Natural","Present-Natural","Present-Natural","Present-Natural","Present-Natural","Present-Natural",
                               "Current","Current","Current","Current","Current","Current")

FD.Continents.GR$Region.names <- c(continents,continents)
str(FD.Continents.GR)
FD.Continents.GR$Scenario <- as.factor(FD.Continents.GR$Scenario)
FD.Continents.GR$Scenario <- relevel(FD.Continents.GR$Scenario, "Present-Natural")
FD.Continents.GR$Region.names <- as.factor(FD.Continents.GR$Region.names)
FD.Continents.GR$Method <- "Geographic Range"
FD.Continents.GR$Region <- "Continent"
colnames(FD.Continents.GR)[1:3] <- c("SR","FRic","FDis")


# ----  b) Functional diversity without geographic range size (PRESENCE-ABSENCE) ----

# 1. change community matrix
str(carnivore.continent.range)
carnivore.continent.presence <- carnivore.continent.range
carnivore.continent.presence[carnivore.continent.presence > 1] <- 1
str(carnivore.continent.presence)

# 2. Traits
new.species <- as.data.frame(colnames(carnivore.continent.presence))
colnames(new.species)[1] <- "Binomial.1.2"
traits.sorted <- merge(species, new.species)
row.names(traits.sorted) <- traits.sorted$Binomial.1.2

str(traits.sorted)
row.names(traits.sorted) <- traits.sorted$Binomial.1.2
traits.sorted <- traits.sorted[c(10,12:14)]
str(traits.sorted)

# Transform body mass
traits.sorted$Mass.g <- log10(traits.sorted$Mass.g)

# Scale traits
traits.sorted <- scale(traits.sorted)

# Put matrix and traits into a list to call dbFD
joined.world.function <- list(traits.sorted, carnivore.continent.presence)
names(joined.world.function)[1] <- "traits"
names(joined.world.function)[2] <- "presence"
str(joined.world.function$traits)
str(joined.world.function$presence)

# Set weights
w <- c(1,1/3,1/3,1/3)

# Call dbFD
global.FD.continents.no.GR <- dbFD(joined.world.function$traits, 
                                   joined.world.function$presence,
                                   w,
                                   m = 4,
                                   stand.FRic = TRUE, 
                                   calc.CWM = TRUE, print.pco = TRUE) # only 0.74 R^2 value

#     Extract and look at results of differences between scenarios at the level of continents
FD.Continents.no.GR <- as.data.frame(cbind(global.FD.continents.no.GR$nbsp,
                                           global.FD.continents.no.GR$FRic,
                                           global.FD.continents.no.GR$FDis))

FD.Continents.no.GR$Scenario <- c("Present-Natural","Present-Natural","Present-Natural","Present-Natural","Present-Natural","Present-Natural",
                                  "Current","Current","Current","Current","Current","Current")

FD.Continents.no.GR$Region.names <- continents
str(FD.Continents.no.GR)
FD.Continents.no.GR$Scenario <- as.factor(FD.Continents.no.GR$Scenario)
FD.Continents.no.GR$Scenario <- relevel(FD.Continents.no.GR$Scenario, "Present-Natural")
FD.Continents.no.GR$Region.names <- as.factor(FD.Continents.no.GR$Region.names)
FD.Continents.no.GR$Method <- "P/A"
FD.Continents.no.GR$Region <- "Continent"
colnames(FD.Continents.no.GR)[1:3] <- c("SR","FRic","FDis")


#————————————————————————————————————————————————————————————————————————————————————————————####
#  ---- 3. Global trait space with systematically-checked species ----

# Get taxonomic information
traits <- read.csv("./Data/Phylacine_Traits/Trait_data.csv", header = TRUE)
str(traits)

# First what do the axes/PCs actually relate to 
functional.axes <- global.FD.continents.GR$x.axes[c(1:4)] # All of these will be the same, regardless of the scenario

functional.axes$Binomial.1.2 <- row.names(functional.axes)
functional.axes <- functional.axes[c(5,1:4)]
functional.axes.and.traits <- merge(functional.axes, traits, by.y = "Binomial.1.2")

#correlations of A1, A2, and A3 with diet traits
summary(lm(data = functional.axes.and.traits, A1~log10(Mass.g)))
cor.test(functional.axes.and.traits$A1,log10(functional.axes.and.traits$Mass.g))
A1.mass <- ggplot(functional.axes.and.traits,aes(x=A1,y=log10(Mass.g))) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("") + ylab(~Log[10]~'(body mass)')

summary(lm(data = functional.axes.and.traits, A1~Diet.Plant))
cor.test(functional.axes.and.traits$A1,functional.axes.and.traits$Diet.Plant)
A1.plant <- ggplot(functional.axes.and.traits,aes(x=A1,y=Diet.Plant)) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("") + ylab("Plant (%)")

summary(lm(data = functional.axes.and.traits, A1~Diet.Vertebrate))
cor.test(functional.axes.and.traits$A1,functional.axes.and.traits$Diet.Vertebrate)
A1.vert <- ggplot(functional.axes.and.traits,aes(x=A1,y=Diet.Vertebrate)) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("") + ylab("Vertebrate (%)")

summary(lm(data = functional.axes.and.traits, A1~Diet.Invertebrate))
cor.test(functional.axes.and.traits$A1,functional.axes.and.traits$Diet.Invertebrate)
A1.invert <- ggplot(functional.axes.and.traits,aes(x=A1,y=Diet.Invertebrate)) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("A1 (40%)") + ylab("Invertebrate (%)")

summary(lm(data = functional.axes.and.traits, A2~log10(Mass.g)))
cor.test(functional.axes.and.traits$A2,log10(functional.axes.and.traits$Mass.g))
A2.mass <- ggplot(functional.axes.and.traits,aes(x=A2,y=log10(Mass.g))) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("") + ylab("")

summary(lm(data = functional.axes.and.traits, A2~Diet.Plant))
cor.test(functional.axes.and.traits$A2,functional.axes.and.traits$Diet.Plant)
A2.plant <- ggplot(functional.axes.and.traits,aes(x=A2,y=Diet.Plant)) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("") + ylab("")

summary(lm(data = functional.axes.and.traits, A2~Diet.Vertebrate))
cor.test(functional.axes.and.traits$A2,functional.axes.and.traits$Diet.Vertebrate)
A2.vert <-ggplot(functional.axes.and.traits,aes(x=A2,y=Diet.Vertebrate)) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("") + ylab("") 

summary(lm(data = functional.axes.and.traits, A2~Diet.Invertebrate))
cor.test(functional.axes.and.traits$A2,functional.axes.and.traits$Diet.Invertebrate)
A2.invert <- ggplot(functional.axes.and.traits,aes(x=A2,y=Diet.Invertebrate)) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("A2 (19%)") + ylab("")

summary(lm(data = functional.axes.and.traits, A3~log10(Mass.g)))
cor.test(functional.axes.and.traits$A3, log10(functional.axes.and.traits$Mass.g))
A3.mass <- ggplot(functional.axes.and.traits,aes(x=A3,y=log10(Mass.g))) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("") + ylab("")

summary(lm(data = functional.axes.and.traits, A3~Diet.Plant))
cor.test(functional.axes.and.traits$A3, functional.axes.and.traits$Diet.Plant)
A3.plant <- ggplot(functional.axes.and.traits,aes(x=A3,y=Diet.Plant)) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("") + ylab("")

summary(lm(data = functional.axes.and.traits, A3~Diet.Vertebrate))
cor.test(functional.axes.and.traits$A3, functional.axes.and.traits$Diet.Vertebrate)
A3.vert <- ggplot(functional.axes.and.traits,aes(x=A3,y=Diet.Vertebrate)) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("") + ylab("")

summary(lm(data = functional.axes.and.traits, A3~Diet.Invertebrate))
cor.test(functional.axes.and.traits$A3, functional.axes.and.traits$Diet.Invertebrate)
A3.invert <-ggplot(functional.axes.and.traits,aes(x=A3,y=Diet.Invertebrate)) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("A3 (9%)") + ylab("")

summary(lm(data = functional.axes.and.traits, A4~log10(Mass.g)))
cor.test(functional.axes.and.traits$A4, log10(functional.axes.and.traits$Mass.g))
A4.mass<- ggplot(functional.axes.and.traits,aes(x=A4,y=log10(Mass.g))) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("") + ylab("")

summary(lm(data = functional.axes.and.traits, A4~Diet.Plant))
cor.test(functional.axes.and.traits$A4, functional.axes.and.traits$Diet.Plant)
A4.plant <- ggplot(functional.axes.and.traits,aes(x=A4,y=Diet.Plant)) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("") + ylab("")

summary(lm(data = functional.axes.and.traits, A4~Diet.Vertebrate))
cor.test(functional.axes.and.traits$A4, functional.axes.and.traits$Diet.Vertebrate)
A4.vert <- ggplot(functional.axes.and.traits,aes(x=A4,y=Diet.Vertebrate)) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("") + ylab("")

summary(lm(data = functional.axes.and.traits, A4~Diet.Invertebrate))
cor.test(functional.axes.and.traits$A4, functional.axes.and.traits$Diet.Invertebrate)
A4.invert <- ggplot(functional.axes.and.traits,aes(x=A4,y=Diet.Invertebrate)) +
  geom_point(alpha = .2) +
  geom_smooth(colour = "black", method = "lm") +
  theme_bw() + xlab("A4 (5%)") + ylab("")

tiff("./Results/Supplementary materials/Figures/PCs_Traits_Data2_Updated_Systematic.tif", height = 6700, width = 8700, res = 800)
grid.arrange(A1.mass, A2.mass, A3.mass, A4.mass,
             A1.plant,A2.plant, A3.plant, A4.plant,
             A1.vert, A2.vert, A3.vert, A4.vert,
             A1.invert, A2.invert, A3.invert, A4.invert,
             nrow = 4, ncol = 4)
dev.off()

# PC1 = Small, insectivorous species (+ve) or larger herbivorous (-ve)
# PC2 = herbivores (+ve) or carnivores (-ve)
# PC3 = larger insectivores (+ve) or medium/small carnivorous/herbivorous (-ve)?
# PC4 = More generalist

#————————————————————————————————————————————————————————————————————————————————————————————####
# ---- 4. Continental trait space with systematically-checked species -----

# ---- A. Continental with geographic range ----
# Get species range per continent
carnivore.continent.range$Scenario <- c("Present-Natural","Present-Natural","Present-Natural","Present-Natural","Present-Natural","Present-Natural",
                                        "Current","Current","Current","Current","Current","Current")
carnivore.continent.range$Continent <- c(continents,continents)
carnivore.continent.range <- carnivore.continent.range[c((length(carnivore.continent.range)-1):(length(carnivore.continent.range)),
                                                         1:(length(carnivore.continent.range)-2))]

# gather to go from wide to long
carnivore.continent.range_long <- gather(carnivore.continent.range, Binomial.1.2, Range.Size,3:length(carnivore.continent.range),factor_key = TRUE)
carnivore.range.pca.taxonomy <- merge(carnivore.continent.range_long,
                                      functional.axes.and.traits, by.y = "Binomial.1.2", all = TRUE)
str(carnivore.range.pca.taxonomy)

# Remove absences
continental.carnivore.trait.space <- carnivore.range.pca.taxonomy[carnivore.range.pca.taxonomy$Range.Size > 0,]
continental.carnivore.trait.space$Continent <- factor(continental.carnivore.trait.space$Continent, levels = c("N. America","Europe","Asia",
                                                                                                              "S. America","Africa","Australia"))
# Calculate the shifted centroid which calculates FDis in FD so the thing which is
# actually representative of the metric

continental.carnivore.trait.space.test <- continental.carnivore.trait.space

continental.carnivore.trait.space.test$centroid.A1 <- continental.carnivore.trait.space.test$Range.Size*continental.carnivore.trait.space.test$A1
continental.carnivore.trait.space.test$centroid.A2 <- continental.carnivore.trait.space.test$Range.Size*continental.carnivore.trait.space.test$A2
continental.carnivore.trait.space.test$centroid.A3 <- continental.carnivore.trait.space.test$Range.Size*continental.carnivore.trait.space.test$A3
continental.carnivore.trait.space.test$centroid.A4 <- continental.carnivore.trait.space.test$Range.Size*continental.carnivore.trait.space.test$A4

centroid.location.A1 <- continental.carnivore.trait.space.test %>% group_by(Scenario,Continent) %>% summarise(centroid.A1 = sum(centroid.A1)/sum(Range.Size))
centroid.location.A2 <- continental.carnivore.trait.space.test %>% group_by(Scenario,Continent) %>% summarise(centroid.A2 = sum(centroid.A2)/sum(Range.Size))
centroid.location.A3 <- continental.carnivore.trait.space.test %>% group_by(Scenario,Continent) %>% summarise(centroid.A3 = sum(centroid.A3)/sum(Range.Size))
centroid.location.A4 <- continental.carnivore.trait.space.test %>% group_by(Scenario,Continent) %>% summarise(centroid.A4 = sum(centroid.A4)/sum(Range.Size))

centroids <- cbind(centroid.location.A1, 
                   centroid.location.A2,
                   centroid.location.A3,
                   centroid.location.A4)

centroids <- centroids[c(1:3,6,9,12)]
colnames(centroids)[3:6] <- c("A1","A2","A3","A4")

# Need to add functional richness - get polygons 
find_hull <- function(df) df[chull(df$A1, df$A2), ]

master.hulls <- list()

for (i in 1:length(scenario)) {
  continental.hulls <- data.frame()
  hull.data <- data.frame()
  scenarionew <- scenario[[i]]
  list.number <- scenario[[i]]
  x <- continental.carnivore.trait.space.test[continental.carnivore.trait.space.test$Scenario == scenario[[i]],]
  
  
  for (j in 1:length(continents)) {
    y <- x[x$Continent == continents[[j]],]
    y.hulls <- plyr::ddply(y,"Continent", find_hull)
    hull.data <- rbind(y.hulls, hull.data)
    print(j)
  }
  
  hull.data$scenarios <- scenarionew
  master.hulls[[list.number]] <- hull.data
}

# Created muster hulls for all scenarios and for each continent
master.hullage <- rbind(master.hulls$Current,
                        master.hulls$`Present-Natural`)
master.hullage$Scenario <- as.factor(master.hullage$Scenario)
str(master.hullage$Scenario)
master.hullage$Scenario <- relevel(master.hullage$Scenario, "Present-Natural")

# ---- Figure S10a: Contienntal trait space with systematically-reviewed species ----
tiff("./Results/Supplementary materials/Figures/Continental_Trait Space.GR_Systematic species.tif", height = 3700, width = 6830, res = 800)
p <- ggplot(continental.carnivore.trait.space.test, aes(x=A1, y = A2, colour = Scenario, size = Range.Size)) +
  geom_point(alpha = 0.15) +
  geom_density_2d(size = .1, alpha = .5) +
  scale_colour_manual(values = c("#339900","#660000")) +
  #scale_colour_manual(values = c("#999999","#000000")) +
  geom_point(data = centroids, size = 2, alpha = 1, shape = 3, stroke = 2) +
  geom_polygon(data = master.hullage, aes(x = A1, y = A2, fill = scenarios), 
               alpha = 0, size = 0.6, linetype = 2) +
  scale_fill_discrete(guide = FALSE) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid=element_blank(),
        strip.text = element_text(face="bold", size=9,lineheight=5.0),
        strip.background = element_rect(fill= c("black","white"), colour="black",
                                        size=1)) +
  facet_wrap(.~Continent, ncol = 3) +
  xlab("A1 (40%)") + ylab("A2 (19%)")

g <- ggplot_gtable(ggplot_build(p))
strip_t <- which(grepl('strip-t', g$layout$name))
fills <- c("#FFCC00","#CCCC99","#FFCC99","#66CC99","#CC99CC","#99CCFF")
k <- 1
for (i in strip_t) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)
dev.off()

tiff("./Results/Figures/Updated continent trait space_Data2_Updated.tif", height = 3000, width = 5000, res = 800)
ggarrange(g,
          labels = c("a"))
dev.off()


# ---- B. Continental without geographic range with systematic-checked species -----

#   Repeat steps: data processing using presence-absence community matrix
str(carnivore.continent.presence)
carnivore.continent.presence$Scenario <- c("Present-Natural","Present-Natural","Present-Natural","Present-Natural","Present-Natural","Present-Natural",
                                           "Current","Current","Current","Current","Current","Current")
continents <- c("Africa","Asia","Australia","Europe","N. America","S. America")
carnivore.continent.presence$Continent <- c(continents,continents)
carnivore.continent.presence <- carnivore.continent.presence[c((length(carnivore.continent.presence)-1):(length(carnivore.continent.presence)),
                                                               1:(length(carnivore.continent.presence)-2))]

# gather to go from wide to long
carnivore.continent.presence.long <- gather(carnivore.continent.presence,
                                            Binomial.1.2, Range.Size,3:length(carnivore.continent.presence), factor_key = TRUE)
carnivore.continent.presence.long.taxonomy <- merge(carnivore.continent.presence.long,
                                                    functional.axes.and.traits, by.y = "Binomial.1.2", all = TRUE)
carnivore.continent.presence.long.taxonomy <- carnivore.continent.presence.long.taxonomy[carnivore.continent.presence.long.taxonomy$Range.Size > 0,]
carnivore.continent.presence.long.taxonomy$Continent <- factor(carnivore.continent.presence.long.taxonomy$Continent, levels = c("N. America","Europe","Asia",
                                                                                                                                "S. America","Africa","Australia"))

# Calculate shifted centroid to look at how FDis was calculated
centroid.location.A1 <- carnivore.continent.presence.long.taxonomy %>% group_by(Scenario,Continent) %>% summarise(centroid.A1 = sum(A1)/sum(Range.Size))
centroid.location.A2 <- carnivore.continent.presence.long.taxonomy %>% group_by(Scenario,Continent) %>% summarise(centroid.A2 = sum(A2)/sum(Range.Size))
centroid.location.A3 <- carnivore.continent.presence.long.taxonomy %>% group_by(Scenario,Continent) %>% summarise(centroid.A3 = sum(A3)/sum(Range.Size))
centroid.location.A4 <- carnivore.continent.presence.long.taxonomy %>% group_by(Scenario,Continent) %>% summarise(centroid.A4 = sum(A4)/sum(Range.Size))

centroids.No.GR <- cbind(centroid.location.A1,
                         centroid.location.A2,
                         centroid.location.A3,
                         centroid.location.A4)

centroids.No.GR <- centroids.No.GR[c(1:3, 6, 9,12)]
colnames(centroids.No.GR)[3:6] <- c("A1","A2","A3","A4")


# plot centroids.No.GR on each plot
tiff("./Results/Supplementary materials/Figures/Continental_Trait Space.No.GR_Systematic species.tif", height = 3700, width = 6830, res = 800)
p <- ggplot(carnivore.continent.presence.long.taxonomy, aes(x=A1, y = A2, colour = Scenario, size = Range.Size)) +
  geom_point(alpha = 0.15, size = 3) +
  geom_density_2d(size = .1, alpha = .5) +
  scale_colour_manual(values = c("#339900","#660000")) +
  geom_point(data = centroids.No.GR, size = 2, alpha = 1, shape = 3, stroke = 2) +
  geom_polygon(data = master.hullage, aes(x = A1, y = A2, fill = scenarios), 
               alpha = 0, size = 0.6, linetype = 2) +
  scale_fill_discrete(guide = FALSE) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        strip.text = element_text(face="bold", size=9,lineheight=5.0),
        strip.background = element_rect(fill="white", colour="black",
                                        size=1)) +
  facet_wrap(.~Continent, ncol = 3) +
  xlab("A1 (40%)") + ylab("A2 (19%)")

g <- ggplot_gtable(ggplot_build(p))
strip_t <- which(grepl('strip-t', g$layout$name))
fills <- c("#FFCC00","#CCCC99","#FFCC99","#66CC99","#CC99CC","#99CCFF")
k <- 1
for (i in strip_t) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)
dev.off()

#————————————————————————————————————————————————————————————————————————————————————————————####
#  ----- 5. Functional diversity differences with systematic species ----

# First step, tidy data which I apparently have not done
Master.Global.FD.Regions <- rbind(FD.Continents.GR,
                                  FD.Continents.no.GR)

# ---- 1. Functional richness - not influenced by weighting ----
fric.change <- Master.Global.FD.Regions[Master.Global.FD.Regions$Method == "P/A",]
fric.change <- fric.change[fric.change$Scenario == "Current",]$FRic - fric.change[fric.change$Scenario == "Present-Natural",]$FRic
median(fric.change)


# Wilcox test:
fric.change <- Master.Global.FD.Regions[Master.Global.FD.Regions$Method == "P/A",]

test <- data.frame(scenario = c("PN","PN","PN","PN","PN","PN", "C","C","C","C","C","C"),
                   value = c(fric.change[fric.change$Scenario == "Present-Natural",]$FRic,
                             fric.change[fric.change$Scenario == "Current",]$FRic))
wilcoxsign_test(test$value ~ test$scenario, alternative = "less")


# Get the relative change in functional richness (% lost from PN)
fric.change <- Master.Global.FD.Regions[Master.Global.FD.Regions$Method == "P/A",]
fric.change.rel <- (fric.change[fric.change$Scenario == "Current",]$FRic - fric.change[fric.change$Scenario == "Present-Natural",]$FRic)/fric.change[fric.change$Scenario == "Present-Natural",]$FRic
median(fric.change.rel)
quantile(fric.change.rel)


# ---- 2a. Functional dispersion with geographic range ----
fdis.change.gr <- Master.Global.FD.Regions[Master.Global.FD.Regions$Method == "Geographic Range",]
fdis.change.gr <- fdis.change.gr[fdis.change.gr$Scenario == "Current",]$FDis - fdis.change.gr[fdis.change.gr$Scenario == "Present-Natural",]$FDis
median(fdis.change.gr)

# Wilcoxon test:
fdis.change.gr <- Master.Global.FD.Regions[Master.Global.FD.Regions$Method == "Geographic Range",]

test <- data.frame(scenario = c("PN","PN","PN","PN","PN","PN", "C","C","C","C","C","C"),
                   value = c(fdis.change.gr[fdis.change.gr$Scenario == "Present-Natural",]$FDis,
                             fdis.change.gr[fdis.change.gr$Scenario == "Current",]$FDis))
wilcoxsign_test(test$value ~ test$scenario, alternative = "less")

# Relative change
fdis.change.gr <- Master.Global.FD.Regions[Master.Global.FD.Regions$Method == "Geographic Range",]
fdis.change.gr <- (fdis.change.gr[fdis.change.gr$Scenario == "Current",]$FDis - fdis.change.gr[fdis.change.gr$Scenario == "Present-Natural",]$FDis)/fdis.change.gr[fdis.change.gr$Scenario == "Present-Natural",]$FDis
median(fdis.change.gr)


# ---- 2b. Functional dispersion without geographic range ----
fdis.change.pa <- Master.Global.FD.Regions[Master.Global.FD.Regions$Method == "P/A",]
fdis.change.pa <- fdis.change.pa[fdis.change.pa$Scenario == "Current",]$FDis - fdis.change.pa[fdis.change.pa$Scenario == "Present-Natural",]$FDis

# Wilcoxon tests:
fdis.change.pa <- Master.Global.FD.Regions[Master.Global.FD.Regions$Method == "P/A",]

test <- data.frame(scenario = c("PN","PN","PN","PN","PN","PN", "C","C","C","C","C","C"),
                   value = c(fdis.change.pa[fdis.change.pa$Scenario == "Present-Natural",]$FDis,
                             fdis.change.pa[fdis.change.pa$Scenario == "Current",]$FDis))
wilcoxsign_test(test$value ~ test$scenario, alternative = "less")


fdis.change.pa <- Master.Global.FD.Regions[Master.Global.FD.Regions$Method == "P/A",]
fdis.change.pa <- (fdis.change.pa[fdis.change.pa$Scenario == "Current",]$FDis - fdis.change.pa[fdis.change.pa$Scenario == "Present-Natural",]$FDis)/fdis.change.pa[fdis.change.pa$Scenario == "Present-Natural",]$FDis
median(fdis.change.pa)
quantile(fdis.change.pa)

# ---- Figure S10(b,c): Functional diversity differences -----
FRic.Change.Continents <- Master.Global.FD.Regions[Master.Global.FD.Regions$Region == "Continent",]
FRic.Change.Continents.GR <- Master.Global.FD.Regions[Master.Global.FD.Regions$Method == "Geographic Range",]
FRic.Change.Continents.PA <- Master.Global.FD.Regions[Master.Global.FD.Regions$Method == "P/A",]

FRic.Change.Continents.GR$Region.names <- factor(FRic.Change.Continents.GR$Region.names, 
                                                 levels = c("N. America","S. America",
                                                            "Europe","Africa",
                                                            "Asia","Australia"))

good.colours <- c("#FFCC99","#99CCFF","#CCCC99","#CC99CC","#FFCC00","#66CC99")
good.colours <- c("#66CC99","#FFCC00","#CC99CC","#CCCC99","#99CCFF","#FFCC99")


levels(FRic.Change.Continents.GR$Scenario)[levels(FRic.Change.Continents.GR$Scenario) == "Present-Natural"] <- "PN"
levels(FRic.Change.Continents.GR$Scenario)[levels(FRic.Change.Continents.GR$Scenario) == "Current"] <- "C"

# Functional richness plot
p1 <- ggplot(FRic.Change.Continents.GR, aes(y = FRic, x = Scenario, colour = Region.names, group = Scenario)) +
  geom_point(size = 3, aes(colour = Region.names)) +
  scale_colour_manual(values = good.colours, name = "Continent") +
  geom_line(aes(group = Region.names), colour = "black", linetype = "dashed") +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())

FDis.Change.Continents <- Master.Global.FD.Regions[Master.Global.FD.Regions$Region == "Continent",]

FDis.Change.Continents$Scenario <- factor(FDis.Change.Continents$Scenario, levels = c("Present-Natural","Current"))

levels(FDis.Change.Continents$Scenario)[levels(FDis.Change.Continents$Scenario) == "Present-Natural"] <- "PN"
levels(FDis.Change.Continents$Scenario)[levels(FDis.Change.Continents$Scenario) == "Current"] <- "C"

FDis.Change.Continents$Region.names <- factor(FDis.Change.Continents$Region.names, 
                                              levels = c("N. America","S. America",
                                                         "Europe","Africa",
                                                         "Asia","Australia"))
# Functional dispersion plot
p2 <- ggplot(FDis.Change.Continents, aes(y = FDis, x = Scenario, colour = Region.names, group = Scenario)) +
  geom_point(size = 3, aes(colour = Region.names, shape = Method)) +
  scale_colour_manual(values = good.colours, name = "Continent") +
  geom_line(aes(group = Region.names), colour = "black", linetype = "dashed") +
  facet_wrap(. ~ Method) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())

tiff("./Results/Supplementary materials/Figures/Continental.FRic.FDis.Change_Systematic species.tif", height = 2000, width = 5000, res = 800)
ggarrange(p1 + theme(legend.position="none",
                     axis.title.x = element_blank()),
          # axis.title = element_text(size = 12),
          #axis.text = element_text(size = 10)),
          p2 + theme(axis.title.x = element_blank(),
                     #axis.title = element_text(size = 12),
                     #axis.text = element_text(size = 10),
                     legend.position = "none"),
          labels = c("b","c"),
          ncol = 2, nrow = 1)
dev.off()



#————————————————————————————————————————————————————————————————————————————————————————————####
# ----- 6. Magnitude of change of centroid from the present-natural to the current: -----

# ---- With geographic range ----
str(centroids)

continents <- c("N. America","Europe","Asia","S. America","Africa","Australia")

centroids.2 <- cbind(centroids[1:6,], centroids[7:12,])
centroids.2 <- centroids.2[c(1:6,9:12)]
centroids.2$Change.A1 <- centroids.2$A1-centroids.2$A11
centroids.2$Change.A2 <- centroids.2$A2-centroids.2$A21
centroids.2$Change.A3 <- centroids.2$A3-centroids.2$A31
centroids.2$Change.A4 <- centroids.2$A4-centroids.2$A41

centroid.position.change <- centroids.2[c(2,11:14)]
centroid.position.change$Method <- "Geographical Range"

centroid.change.continent.GR <- as.data.frame(cbind(centroids.2$Continent, centroids.2$Change.A1))
colnames(centroid.change.continent.GR) <- c("Continent","Change")
centroid.change.continent.GR$Continent <- continents
centroid.change.continent.GR$A <- "1"

centroid.change.continent.GR.2 <- as.data.frame(cbind(centroids.2$Continent, centroids.2$Change.A2))
colnames(centroid.change.continent.GR.2) <- c("Continent","Change")
centroid.change.continent.GR.2$Continent <- continents
centroid.change.continent.GR.2$A <- "2"

centroid.change.continent.GR.3 <- as.data.frame(cbind(centroids.2$Continent, centroids.2$Change.A3))
colnames(centroid.change.continent.GR.3) <- c("Continent","Change")
centroid.change.continent.GR.3$Continent <- continents
centroid.change.continent.GR.3$A <- "3"

centroid.change.continent.GR.4 <- as.data.frame(cbind(centroids.2$Continent, centroids.2$Change.A4))
colnames(centroid.change.continent.GR.4) <- c("Continent","Change")
centroid.change.continent.GR.4$Continent <- continents
centroid.change.continent.GR.4$A <- "4"

centroid.change.continent.GR <- rbind(centroid.change.continent.GR, 
                                      centroid.change.continent.GR.2,
                                      centroid.change.continent.GR.3,
                                      centroid.change.continent.GR.4)

centroid.change.continent.GR$Method <- "Geographic range"

# ---- Without geographic range ----
str(centroids)

centroids.2 <- cbind(centroids.No.GR[1:6,], centroids.No.GR[7:12,])
centroids.2 <- centroids.2[c(1:6,9:12)]
centroids.2$Change.A1 <- centroids.2$A1-centroids.2$A11
centroids.2$Change.A2 <- centroids.2$A2-centroids.2$A21
centroids.2$Change.A3 <- centroids.2$A3-centroids.2$A31
centroids.2$Change.A4 <- centroids.2$A4-centroids.2$A41

centroid.position.change.2 <- centroids.2[c(2,11:14)]
centroid.position.change.2$Method <- "Presence/Absence"
Master.centroid.position.change <- rbind(centroid.position.change,
                                         centroid.position.change.2)

centroid.change.continent.no.GR <- as.data.frame(cbind(centroids.2$Continent, centroids.2$Change.A1))
colnames(centroid.change.continent.no.GR) <- c("Continent","Change")
centroid.change.continent.no.GR$Continent <- continents
centroid.change.continent.no.GR$A <- "1"

centroid.change.continent.no.GR.2 <- as.data.frame(cbind(centroids.2$Continent, centroids.2$Change.A2))
colnames(centroid.change.continent.no.GR.2) <- c("Continent","Change")
centroid.change.continent.no.GR.2$Continent <- continents
centroid.change.continent.no.GR.2$A <- "2"

centroid.change.continent.no.GR.3 <- as.data.frame(cbind(centroids.2$Continent,centroids.2$Change.A3))
colnames(centroid.change.continent.no.GR.3) <- c("Continent","Change")
centroid.change.continent.no.GR.3$Continent <- continents
centroid.change.continent.no.GR.3$A <- "3"

centroid.change.continent.no.GR.4 <- as.data.frame(cbind(centroids.2$Continent,centroids.2$Change.A4))
colnames(centroid.change.continent.no.GR.4) <- c("Continent","Change")
centroid.change.continent.no.GR.4$Continent <- continents
centroid.change.continent.no.GR.4$A <- "4"

centroid.change.continent.no.GR <- rbind(centroid.change.continent.no.GR,
                                         centroid.change.continent.no.GR.2,
                                         centroid.change.continent.no.GR.3,
                                         centroid.change.continent.no.GR.4)

centroid.change.continent.no.GR$Method <- "Presence/Absence"


Master.Centroid.Shift <- rbind(centroid.change.continent.GR, centroid.change.continent.no.GR)

Master.Centroid.Shift$Continent <- factor(Master.Centroid.Shift$Continent, levels = c("N. America","S. America","Europe","Africa","Asia","Australia"))
Master.Centroid.Shift$Method <- factor(Master.Centroid.Shift$Method, levels = c("Presence/Absence","Geographic range"))

Master.centroid.position.change$Method <- factor(Master.centroid.position.change$Method, levels = c("Presence/Absence","Geographical Range"))

# ---- Statistics ----

# Kruskal-Wallis
# 1. Change between principle components
Master.Centroid.Shift.PA <- Master.Centroid.Shift[Master.Centroid.Shift$Method == "Presence/Absence",]
Master.Centroid.Shift.PA$A <- as.factor(Master.Centroid.Shift.PA$A)
kruskal.test(Master.Centroid.Shift.PA$Change, Master.Centroid.Shift.PA$A)
dunn.test(Master.Centroid.Shift.PA$Change,Master.Centroid.Shift.PA$A)

Master.Centroid.Shift.GR <- Master.Centroid.Shift[Master.Centroid.Shift$Method == "Geographic range",]
Master.Centroid.Shift.GR$A <- as.factor(Master.Centroid.Shift.GR$A)
kruskal.test(Master.Centroid.Shift.GR$Change, Master.Centroid.Shift.GR$A)
dunn.test(Master.Centroid.Shift.GR$Change,Master.Centroid.Shift.GR$A)


# Mann-Whitney U-tests
# 2. Change within principle components

# PCoA1
wilcox_test(Master.centroid.position.change$Change.A1 ~ Master.centroid.position.change$Method)
tapply(Master.centroid.position.change$Change.A1, Master.centroid.position.change$Method, median)

(u1 <- (6*6) + (6*(6+1))/2 - sum(rank(Master.centroid.position.change$Change.A1)[1:6]))
(u2 <- (6*6) + (6*(6+1))/2 - sum(rank(Master.centroid.position.change$Change.A1)[7:12]))

# PCoA2
wilcox_test(Master.centroid.position.change$Change.A2 ~ Master.centroid.position.change$Method)
tapply(Master.centroid.position.change$Change.A2, Master.centroid.position.change$Method, median)

(u1 <- (6*6) + (6*(6+1))/2 - sum(rank(Master.centroid.position.change$Change.A2)[1:6]))
(u2 <- (6*6) + (6*(6+1))/2 - sum(rank(Master.centroid.position.change$Change.A2)[7:12]))

# PCoA3
wilcox_test(Master.centroid.position.change$Change.A3 ~ Master.centroid.position.change$Method)
tapply(Master.centroid.position.change$Change.A3, Master.centroid.position.change$Method, median)

(u1 <- (6*6) + (6*(6+1))/2 - sum(rank(Master.centroid.position.change$Change.A3)[1:6]))
(u2 <- (6*6) + (6*(6+1))/2 - sum(rank(Master.centroid.position.change$Change.A3)[7:12]))

# PCoA4
wilcox_test(Master.centroid.position.change$Change.A4 ~ Master.centroid.position.change$Method)
tapply(Master.centroid.position.change$Change.A4, Master.centroid.position.change$Method, median)

(u1 <- (6*6) + (6*(6+1))/2 - sum(rank(Master.centroid.position.change$Change.A4)[1:6]))
(u2 <- (6*6) + (6*(6+1))/2 - sum(rank(Master.centroid.position.change$Change.A4)[7:12]))

# Plots
good.colours <- c("#66CC99","#CC99CC","#99CCFF","#FFCC00","#CCCC99","#FFCC99")
Master.Centroid.Shift$Method <- factor(Master.Centroid.Shift$Method, levels = c("Geographic range", "Presence/Absence"))


# ---- Figure S11: Centroid directional shifts for systematically reviewed species ----
p2 <- ggplot(Master.Centroid.Shift,
             aes(x=reorder(A, desc(A)), y = Change, group = Method, shape = Method)) +
  # stat_summary(fun.data=median_hilow,
  # stat_summary(fun.y = "mean", 
  #              fun.ymin = function(x) mean(x) - sd(x), 
  #             fun.ymax = function(x) mean(x) + sd(x), 
  #              geom = "pointrange", position = position_dodge(width = 0.9), size = 1) +
  geom_point(aes(colour = Continent),
             position = position_jitterdodge(dodge.width = 0.9,
                                             jitter.width = 0.2), size = 2) +
  geom_pointrange(mapping = aes(x = reorder(A, desc(A)), y = Change),
                  stat = "summary",
                  fun.ymin = function(z) {quantile(z,0.25)},
                  fun.ymax = function(z) {quantile(z,0.75)},
                  fun.y = median, 
                  position = position_jitterdodge(dodge.width = 0.5, jitter.width = 0.2),
                  colour = "black",
                  size = .7) +
  scale_colour_manual(values = good.colours) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1, colour = "grey") +
  theme_bw() +
  coord_flip() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15)) +
  geom_vline(xintercept = c(1.5,2.5,3.5), colour = "black") +
  geom_vline(xintercept = c(1,2,3,4), colour = "grey", linetype = 2) +
  xlab("Principal coordinate axis") + ylab("Centroid shift") +
  ylim(-0.08,0.15)

#Save without the centroid positional plot
tiff("./Results/Supplementary materials/Figures/Continental.Centroid.Magnitude.Shift.Figures_Systematic.tif", height = 4000, width = 6500, res = 800)
p2
dev.off()


#————————————————————————————————————————————————————————————————————————————————————————————####
#  ----- Additional supplementary materials -----
#————————————————————————————————————————————————————————————————————————————————————————————####

# ---- Bias in range fragmentation/size difference between the PN and current ----
master.range.stats <- read.csv("./Data/SpeciesRanges_Metrics and Traits.csv")
colnames(master.range.stats)

#subset to include data needed
range.bias.data <- master.range.stats[c(2:6,33,8,20,18,30:32)]
range.bias.data <- merge(range.bias.data, traits[c(1,2)], by = "Binomial.1.2")
colnames(range.bias.data)


# ---- Figure S12: Trait correlates with geographic range metrics ----
(p1 <- ggplot(range.bias.data, aes(x = Diet.Vertebrate, y = Current_Clumps, 
                                              colour = Order.1.2, group = Order.1.2)) +
  geom_point()+
  theme_bw() +
  scale_colour_manual(values = colours, name = "Order") +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  xlab("Vertebrate in diet (%)") + ylab(""))

(p2 <- ggplot(range.bias.data, aes(x = log(Mass.g), y = Current_Clumps,
                                              colour = Order.1.2, group = Order.1.2)) +
  geom_point() +
  theme_bw() +
  scale_colour_manual(values = colours, name = "Order") +
  xlab("Body mass (log(g))") + ylab(""))


(p3 <- ggplot(range.bias.data, aes(x = References, y = Current_Clumps,
                                  colour = Order.1.2, group = Order.1.2)) +
  geom_point() +
  theme_bw() +
  scale_colour_manual(values = colours, name = "Order") +
  xlab("References") + ylab(""))


(master.p <- ggarrange(p1 + theme(legend.position="none"),
                      p2 + theme(legend.position="none"),
                      p3+ theme(legend.position="none"), ncol = 1, nrow = 3, 
                      common.legend = TRUE, 
                      legend = "right",
                      labels = c("a","b","c")))

tiff("./Results/Supplementary materials/Figures/Range Fragmentation and Traits.tif", height = 5500, width = 4500, res = 800)
annotate_figure(master.p,
                left = text_grob("Range fragmentation (Number of isloated patched)", rot = 90))
dev.off()

# ---- Principal component analyses on range change metrics ----
colnames(range.bias.data)

# Get current range metrics
pca.current.ranges <- range.bias.data[c(8,10,12)]
row.names(pca.current.ranges) <- range.bias.data$Binomial.1.2
pca.current.ranges <- na.omit(pca.current.ranges)

# Get present-natural range metrics
pca.presnat.ranges <- range.bias.data[c(7,9,11)]
row.names(pca.presnat.ranges) <- range.bias.data$Binomial.1.2
pca.presnat.ranges <- na.omit(pca.presnat.ranges)

#Look at these in a PCA
pca.change.ranges <- data.frame(range.bias.data$Binomial.1.2)
pca.change.ranges$Binomial.1.2 <- range.bias.data$Binomial.1.2
pca.change.ranges$size_change <- range.bias.data$n.cell-range.bias.data$presnat_n.cell
pca.change.ranges$CAI_change <- range.bias.data$core.area.index-range.bias.data$presnat_core.area.index
pca.change.ranges$clumps_change <- range.bias.data$Current_Clumps-range.bias.data$PN_Clumps

row.names(pca.change.ranges) <- range.bias.data$Binomial.1.2
pca.change.ranges <- pca.change.ranges[,-c(1,2)]
pca.change.ranges <- na.omit(pca.change.ranges)

pca.3 <- prcomp(scale(pca.change.ranges))
summary(pca.3)
print(pca.3)
pca.axes <- as.data.frame(pca.3$x)
pca.axes$Binomial.1.2 <- row.names(pca.axes)
autoplot(prcomp(scale(pca.change.ranges)),
         loadings = TRUE,
         loadings.label = TRUE,loadings.label.size = 5)

# This figure potentially shows how species ranges can change in different ways, and how 
# add vertebrate and mass to this graph
pca.axes <- merge(pca.axes,species, by.y = "Binomial.1.2")
(p4 <- ggplot(pca.axes, aes(x = PC1, y = PC2, size = Mass.g, colour = Diet.Vertebrate,
                           label = Binomial.1.2)) +
  geom_point() +
  #geom_label() +
  xlab("PC1 (69.4%)") + ylab("PC2 (18.6%)") +
  theme_bw())
ggplotly(p4)

ggsave("./Results/Supplementary materials/Figures/Change in georgaphic range in PCA.pdf", p4,
       width = 4.5, height = 4.5, units = "in")

# ---- Species' geographic range differences ----

# Create function to extract raster values from data frame and convert back into a dataframe.
create.raster <- function(x, y) {
  if(y == "current") {
    which(colnames(current.range.full) %in% x)
    test <- current.range.full[c(which(colnames(current.range.full) %in% x))]
    m <- as.matrix(matrix(test[,1], ncol = 360, nrow = 142, byrow = TRUE))
    r <<- raster(m)
    crs(r) <<- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    extent(r) <<- c(-17367530, 17367530, -6301393, 7296404)
  } else {
    which(colnames(presnat.range.full) %in% x)
    test <- presnat.range.full[c(which(colnames(presnat.range.full) %in% x))]
    m <- as.matrix(matrix(test[,1], ncol = 360, nrow = 142, byrow = TRUE))
    r <<- raster(m)
    crs(r) <<- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    extent(r) <<- c(-17367530, 17367530, -6301393, 7296404)
  }
}


# Example species: Jaguar
create.raster("Panthera_onca", "pn")
plot(r)


#  1. Leopard
par(mfrow = c(1,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

create.raster("Panthera_pardus", "current")
plot(r)
plot(continents, add = TRUE)

create.raster("Panthera_pardus", "pn")
leopard.pn <- r
plot(leopard.pn)

create.raster("Panthera_pardus", "current")
leopard.c <- r
str(leopard.c)

# both ranges on the same plot
master.leopard <- cbind(as.data.frame(leopard.pn), as.data.frame(leopard.c))
master.leopard <- master.leopard[1] + master.leopard[2]

master.leopard <- matrix(master.leopard$layer,
                         ncol = 360, nrow = 142, byrow = TRUE)
master.leopard <- raster(master.leopard)
str(master.leopard)
projection(master.leopard) <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
extent(master.leopard) <- c(-17367530, 17367530, -6301393, 7296404)

plot(master.leopard, axes = F, box = F, legend = F,
     col = c("#FFFFFF","#996633","#3366CC"))
plot(continents, add = TRUE)


# 2. Spotted hyaena
par(mfrow = c(1,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

create.raster("Crocuta_crocuta", "current")
plot(r)
plot(continents, add = TRUE)

create.raster("Crocuta_crocuta", "pn")
leopard.pn <- r
plot(leopard.pn)

create.raster("Crocuta_crocuta", "current")
leopard.c <- r
str(leopard.c)

# both ranges on the same plot
master.hyena <- cbind(as.data.frame(leopard.pn), as.data.frame(leopard.c))
master.hyena <- master.hyena[1] + master.hyena[2]

master.hyena <- matrix(master.hyena$layer,
                         ncol = 360, nrow = 142, byrow = TRUE)
master.hyena <- raster(master.hyena)
str(master.hyena)
projection(master.hyena) <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
extent(master.hyena) <- c(-17367530, 17367530, -6301393, 7296404)

plot(master.hyena, axes = F, box = F, legend = F,
     col = c("#FFFFFF","#996633","#3366CC"))
plot(continents, add = TRUE)


# 3. Wolverine
par(mfrow = c(1,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

create.raster("Gulo_gulo", "current")
plot(r)
plot(continents, add = TRUE)

create.raster("Gulo_gulo", "pn")
leopard.pn <- r
plot(leopard.pn)

create.raster("Gulo_gulo", "current")
leopard.c <- r
str(leopard.c)

# both ranges on the same plot
master.wolverine <- cbind(as.data.frame(leopard.pn), as.data.frame(leopard.c))
master.wolverine <- master.wolverine[1] + master.wolverine[2]

master.wolverine <- matrix(master.wolverine$layer,
                         ncol = 360, nrow = 142, byrow = TRUE)
master.wolverine <- raster(master.wolverine)
str(master.wolverine)
projection(master.wolverine) <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
extent(master.wolverine) <- c(-17367530, 17367530, -6301393, 7296404)

plot(master.wolverine, axes = F, box = F, legend = F,
     col = c("#FFFFFF","#996633","#3366CC"))
plot(continents, add = TRUE)


# 4. Coyote
par(mfrow = c(1,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

create.raster("Canis_latrans", "current")
plot(r)
plot(continents, add = TRUE)

create.raster("Canis_latrans", "pn")
leopard.pn <- r
plot(leopard.pn)

create.raster("Canis_latrans", "current")
leopard.c <- r
str(leopard.c)

# both ranges on the same plot
master.coyote <- cbind(as.data.frame(leopard.pn), as.data.frame(leopard.c))
master.coyote <- master.coyote[1] + master.coyote[2]

master.coyote <- matrix(master.coyote$layer,
                           ncol = 360, nrow = 142, byrow = TRUE)
master.coyote <- raster(master.coyote)
str(master.coyote)
projection(master.coyote) <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
extent(master.coyote) <- c(-17367530, 17367530, -6301393, 7296404)

plot(master.coyote, axes = F, box = F, legend = F,
     col = c("#FFFFFF","#CCCCCC","#000000"))
plot(continents, add = TRUE)


# 5. Orangutan
par(mfrow = c(1,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

create.raster("Pongo_pygmaeus", "current")
plot(r)
plot(continents, add = TRUE)

create.raster("Pongo_pygmaeus", "pn")
leopard.pn <- r
plot(leopard.pn)

create.raster("Pongo_pygmaeus", "current")
leopard.c <- r
str(leopard.c)

# both ranges on the same plot
master.pongo <- cbind(as.data.frame(leopard.pn), as.data.frame(leopard.c))
master.pongo <- master.pongo[1] + master.pongo[2]

master.pongo <- matrix(master.pongo$layer,
                         ncol = 360, nrow = 142, byrow = TRUE)
master.pongo <- raster(master.pongo)
str(master.pongo)
projection(master.pongo) <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
extent(master.pongo) <- c(-17367530, 17367530, -6301393, 7296404)

plot(master.pongo, axes = F, box = F, legend = F,
     col = c("#FFFFFF","#996633","#3366CC"))
plot(continents, add = TRUE)


# 6. Weasel 
par(mfrow = c(1,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

create.raster("Mustela_nigripes", "current")
plot(r)
plot(continents, add = TRUE)

create.raster("Mustela_nigripes", "pn")
leopard.pn <- r
plot(leopard.pn)

create.raster("Mustela_nigripes", "current")
leopard.c <- r
str(leopard.c)

# both ranges on the same plot
master.weasel <- cbind(as.data.frame(leopard.pn), as.data.frame(leopard.c))
master.weasel <- master.weasel[1] + master.weasel[2]

master.weasel <- matrix(master.weasel$layer,
                         ncol = 360, nrow = 142, byrow = TRUE)
master.weasel <- raster(master.weasel)
str(master.weasel)
projection(master.weasel) <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
extent(master.weasel) <- c(-17367530, 17367530, -6301393, 7296404)

plot(master.weasel, axes = F, box = F, legend = F,
     col = c("#FFFFFF","#996633","#3366CC"))
plot(continents, add = TRUE)
 

# 7. Catopuma
par(mfrow = c(1,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

create.raster("Catopuma_temminckii", "current")
plot(r)
plot(continents, add = TRUE)

create.raster("Catopuma_temminckii", "pn")
leopard.pn <- r
plot(leopard.pn)

create.raster("Catopuma_temminckii", "current")
leopard.c <- r
str(leopard.c)

# both ranges on the same plot
master.catopuma <- cbind(as.data.frame(leopard.pn), as.data.frame(leopard.c))
master.catopuma <- master.catopuma[1] + master.catopuma[2]

master.catopuma <- matrix(master.catopuma$layer,
                         ncol = 360, nrow = 142, byrow = TRUE)
master.catopuma <- raster(master.catopuma)
str(master.catopuma)
projection(master.catopuma) <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
extent(master.catopuma) <- c(-17367530, 17367530, -6301393, 7296404)

plot(master.catopuma, axes = F, box = F, legend = F,
     col = c("#FFFFFF","#996633","#3366CC"))
plot(continents, add = TRUE)


# 8. Lion
par(mfrow = c(1,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

create.raster("Panthera_leo", "current")
plot(r)
plot(continents, add = TRUE)

create.raster("Panthera_leo", "pn")
leopard.pn <- r
plot(leopard.pn)

create.raster("Panthera_leo", "current")
leopard.c <- r
str(leopard.c)

# both ranges on the same plot
master.lion <- cbind(as.data.frame(leopard.pn), as.data.frame(leopard.c))
master.lion <- master.lion[1] + master.lion[2]

master.lion <- matrix(master.lion$layer,
                         ncol = 360, nrow = 142, byrow = TRUE)
master.lion <- raster(master.lion)
str(master.lion)
projection(master.lion) <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
extent(master.lion) <- c(-17367530, 17367530, -6301393, 7296404)

plot(master.lion, axes = F, box = F, legend = F,
     col = c("#FFFFFF","#996633","#3366CC"))
plot(continents, add = TRUE)


# 9. Shrew
par(mfrow = c(1,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

create.raster("Sorex_minutus", "current")
plot(r)
plot(continents, add = TRUE)

create.raster("Sorex_minutus", "pn")
leopard.pn <- r
plot(leopard.pn)

create.raster("Sorex_minutus", "current")
leopard.c <- r
str(leopard.c)

# both ranges on the same plot
master.shrew <- cbind(as.data.frame(leopard.pn), as.data.frame(leopard.c))
master.shrew <- master.shrew[1] + master.shrew[2]

master.shrew <- matrix(master.shrew$layer,
                         ncol = 360, nrow = 142, byrow = TRUE)
master.shrew <- raster(master.shrew)
str(master.shrew)
projection(master.shrew) <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
extent(master.shrew) <- c(-17367530, 17367530, -6301393, 7296404)

plot(master.shrew, axes = F, box = F, legend = F,
     col = c("#FFFFFF","#996633","#3366CC"))
plot(continents, add = TRUE)


# ---- Figure S4: Species range differences ----
par(mfrow = c(3,1))
plot(master.leopard, axes = F, box = F, legend = F, col = c("#FFFFFF","#996633","#3366CC"))
plot(continents, add = TRUE)
plot(master.lion, axes = F, box = F, legend = F, col = c("#FFFFFF","#996633","#3366CC"))
plot(continents, add = TRUE)
plot(master.weasel, axes = F, box = F, legend = F, col = c("#FFFFFF","#996633","#3366CC"))
plot(continents, add = TRUE)

plot(master.hyena, axes = F, box = F, legend = F, col = c("#FFFFFF","#996633","#3366CC"))
plot(continents, add = TRUE)
plot(master.pongo, axes = F, box = F, legend = F, col = c("#FFFFFF","#996633","#3366CC"))
plot(continents, add = TRUE)
plot(master.shrew, axes = F, box = F, legend = F, col = c("#FFFFFF","#996633","#3366CC"))
plot(continents, add = TRUE)

plot(master.wolverine, axes = F, box = F, legend = F, col = c("#FFFFFF","#996633","#3366CC"))
plot(continents, add = TRUE)
plot(master.catopuma, axes = F, box = F, legend = F, col = c("#FFFFFF","#996633","#3366CC"))
plot(continents, add = TRUE)
plot(master.coyote, axes = F, box = F, legend = F, col = c("#FFFFFF","#CCCCCC","#000000"))
plot(continents, add = TRUE)



# ---- Figure S7: Theoretial figure showing calculation of functional diversity metrics ----

made.up.data <- data.frame(c(1,1.5,2.4,4.6,5.6,7.4,7.6,7.8,1,1.5,2.4,4.6,5.6,7.4,7.6, 7.8),
                           c(7.8,4.3,2.4,6.4,1.8,7.1,8.1,3.2,7.8,4.3,2.4,6.4,1.8,7.1,8.1,3.2),
                           c(1,1,1,1,1,1,1,1,5.5,7,6,8,3,0.3,1.5,1.2),
                           c("Unweighted","Unweighted","Unweighted","Unweighted","Unweighted","Unweighted","Unweighted","Unweighted",
                             "Weighted","Weighted","Weighted","Weighted","Weighted","Weighted","Weighted","Weighted"))
colnames(made.up.data) <- c("Trait.1","Trait.2","Weight.value","Weighted")

made.up.data.centroids.unweighted <- made.up.data %>% group_by(Weighted) %>%  summarise(Trait.1 = sum(Trait.1)/sum(Weight.value),
                                                                                        Trait.2 = sum(Trait.2)/sum(Weight.value))

made.up.data.centroids.weighted <- made.up.data[made.up.data$Weighted == "Weighted",]

made.up.data.centroids.weighted$weighted.Trait.1 <- made.up.data.centroids.weighted$Trait.1*made.up.data.centroids.weighted$Weight.value
made.up.data.centroids.weighted$weighted.Trait.2 <- made.up.data.centroids.weighted$Trait.2*made.up.data.centroids.weighted$Weight.value

made.up.data.centroids.weighted <- made.up.data.centroids.weighted %>% summarise(Trait.1 = sum(weighted.Trait.1)/sum(Weight.value),
                                                                                 Trait.2 = sum(weighted.Trait.2)/sum(Weight.value))

str(made.up.data.centroids.weighted)
made.up.data.centroids.weighted$Weighted <- "Weighted"
made.up.data.centroids.weighted$Weight.value <- 5


str(made.up.data.centroids.unweighted)
made.up.data.centroids.unweighted <- made.up.data.centroids.unweighted[made.up.data.centroids.unweighted$Weighted == "Unweighted",]
made.up.data.centroids.unweighted$Weighted <- "Unweighted"
made.up.data.centroids.unweighted$Weight.value <- 5

made.up.data$Weighted <- relevel(made.up.data$Weighted, "Weighted")

find_hull <- function(df) df[chull(df$Trait.1, df$Trait.2),]
y.hulls <- plyr::ddply(made.up.data, "Weighted", find_hull)

df <- data.frame(c(3.291692,3.291692,3.291692,3.291692,3.291692,3.291692,3.291692,3.291692,
                   4.7375,4.7375,4.7375,4.7375,4.7375,4.7375,4.7375,4.7375),
                 c(4.988308,4.988308,4.988308,4.988308,4.988308,4.988308,4.988308,4.988308,
                   5.1375,5.1375,5.1375,5.1375,5.1375,5.1375,5.1375,5.1375),
                 c(1,1.4,2.5,4.6,5.6,7.4,7.6,7.8,1,1.4,2.5,4.6,5.6,7.4,7.6,7.8),
                 c(7.8,4.3,2.4,6.4,1.8,7.1,8.1,3.2,7.8,4.3,2.4,6.4,1.8,7.1,8.1,3.2),
                 c("Weighted","Weighted","Weighted","Weighted","Weighted","Weighted","Weighted","Weighted",
                   "Unweighted","Unweighted","Unweighted","Unweighted","Unweighted","Unweighted","Unweighted","Unweighted"),
                 c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2)) 

colnames(df) <- c("Centroid.x","Centroid.y","Trait.1","Trait.2","Weighted","Weight.value") 


tiff("./Results/Supplementary materials/Figures/theoretical figure.tif", height = 3300, width = 3300, res = 800)
ggplot(made.up.data, aes(x = Trait.1, y = Trait.2, size = Weight.value, colour = Weighted)) +
  geom_point(alpha = 0.5) +
  geom_point(data = made.up.data.centroids.weighted, size = 10) +
  geom_point(data = made.up.data.centroids.unweighted, size = 10) +
  geom_polygon(data = y.hulls, aes(x = Trait.1, y = Trait.2), 
               alpha = 0, size = 0.2, linetype = 2, colour = "Black") +
  geom_segment(aes(x = Centroid.x, y = Centroid.y, xend = Trait.1, yend = Trait.2, size = 0.1), data = df) +
  scale_colour_manual(values = c("#999999","#000000")) +
  theme_classic() +
  xlab("Trait 1") + ylab("Trait 2") +
  theme(legend.position = "none") 
dev.off()                    

#————————————————————————————————————————————————————————————————————————————————————————————####
#  ----- End of script -----
#————————————————————————————————————————————————————————————————————————————————————————————####