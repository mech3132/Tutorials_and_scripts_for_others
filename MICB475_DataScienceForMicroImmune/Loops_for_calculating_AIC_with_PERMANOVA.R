library(phyloseq)
library(tidyverse)
library(vegan)

# Here, you need to "source" code from the file, AIC_for_adonis2.R
# This includes custom scripts that calculate AIC for permanovas
source("AIC_for_adonis2.R")

# Let's say you have a phyloseq object called "phyloseqObject"
# Here, I'm going to use the atacama dataset as an example. I haven't included 
# the code here for it, so you should replace atacama_final with whatever 
# your phyloseq object is called.
load("atacam_final.RData")
phyloseqObject <- atacama_final

## First, extract sample data of interest
sampdat <- data.frame(sample_data(phyloseqObject))
# It is important that we filter out any NAs in the dataset
sampdat_filt <- sampdat %>% select(elevation
                                   , transectname
                                   , depth
                                   , ph
                                   , toc
                                   , averagesoilrelativehumidity
                                   , averagesoiltemperature
                                   , vegetation
) %>%
  drop_na()

# We also need to filter the phyloseq object to remove those dropped samples
phyloseqObject_filt <- prune_samples(rownames(sampdat_filt), phyloseqObject)


# Our goal here is to iterate through all combinations of predictors
# How do we do that?

####### Before the loop #######
# First, let's make a matrix that descrbies all the combination of predictors we want.
# This code below makes a data frame of TRUE's and FALSE's that describe
# ALL combinations of predictors (each row is a combination of predictors, where 1 means the predictor is present and 0 means it is absent)
combinations <- expand.grid(rep(list(c(FALSE,TRUE)),ncol(sampdat_filt))) 
# We actually want to get rid of the row with all FALSE because that model doesn't have any predictors
combinations <- combinations[-1,]

# Now, let's make a distance matrix response variable. This will be the same between all models.
dm <- distance(phyloseqObject_filt, method="bray")

# Finally, make a vector to save all our results (AIC values) 
# We want it the same lengths as there are rows in combinations, since that is what we are looping through
aic_results <- vector(length=nrow(combinations))

##### Start writing our loop ####
# Use the draft code below to make a loop
# We want to go through ROW and use that to subset the metadata. 
# For example, let's say we are on row 10
r <- 10 # setting arbitrary variable to 10
# We extract row 3 from "combinations" and use "which" to tell us 
toinclude <- unlist(combinations[r,])
# If we look at this object, it is a vector with two TRUEs

# Next, we can select the COLUMNS of the sample data with this vector
newSampDat <- data.frame(sampdat_filt[,toinclude])

# Then, we can use all the columns from this new sample data to run a PERMANOVA
permanova_results <- adonis2(dm ~ ., data=newSampDat, permutations = 10000)

# Use the custom sourced AICc.PERMANOVA2 function to calculate AIC value from PERMANOVA
aic <- AICc.PERMANOVA2(permanova_results)

# Save the AIC value into our vector
# We can conveniently use the 'r' variable (the row we're on) to index the aic_results vector
aic_results[r] <- aic$AIC



########## Analyze results ###########
# After you make aic_results

# Find out which is the smallest AIC
minRow <- which.min(aic_results)
# Find the combination of TRUEs and FALSEs result in the smallest AIC
bestCombo <- unlist(combinations[minRow,]) # minRow is an index-- it is getting the minRow-th row of combinations
# Find out the column names where you have TRUE
colnames(data.frame(sampdat_filt[,bestCombo])) 
# The above list of variables is your "best" model!
# 