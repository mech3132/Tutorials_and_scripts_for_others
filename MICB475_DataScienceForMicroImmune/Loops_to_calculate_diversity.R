library(MASS)
library(phyloseq)
library(tidyverse)
library(vegan)
# NOTE: above, it's important that you load MASS before tidyverse.
# The both have functions called 'select', and if you load tidyverse first,
# MASS::select() will mask tidyverse::select(), which will cause bugs below.


# Let's say you have a phyloseq object called "phyloseqObject"
# Here, I'm going to use the atacama dataset as an example. I haven't included 
# the code here for it, so you should replace atacama_final with whatever 
# your phyloseq object is called.
load("atacam_final.RData")
phyloseqObject <- atacama_final

######## Setting up richness loop ########
# Calculate richness metric
estrich <- estimate_richness(phyloseqObject)
# Extract sample data
sampdat <- data.frame(sample_data(phyloseqObject))


# Here, I am filtering metadata so that it is only the metadata of interest
sampdat_filt <- sampdat %>% select(elevation
                                   , transectname
                                   , depth
                                   , ph
                                   , toc
                                   , averagesoilrelativehumidity
                                   , averagesoiltemperature
                                   , vegetation
                                   )


######## Alpha diversity loop to test each individual predictor ###########
# To run linear models for each predictor, 
# you can loop through each column.

# Make vector of column names
allPredictors <- colnames(sampdat_filt)
# Make a vector that is just shannon
shannon_vec <- estrich$Shannon

# For each loop, you are going to run a linear model.
# You need somewhere to save the output, so we are going to make a 'list', which we will populate with the results
resultsList <- list()
for ( x in allPredictors ) {
  # Make linear model
  # Here, the response variable stays the same every time. It is the shannon diversity vector.
  # Then, the function get() asks lm() to treat 'x' as a variable rather than as 'x' itself. 
  # Using get() allows you to 'get' the vaues stored WITHIN x
  model <- lm(shannon_vec ~ get(x), data=sampdat_filt)
  # Save the model summary
  modelsummary <- summary(model)
  # Finally, we are going to save the modelsummary in our list
  # We want to name each list item with the variable 'x', so we do this:
  resultsList[[x]] <- modelsummary
  
}

# Now, if we call our list, "resultsList", it will return ALL model summary results
resultsList

# We can save this as a text file for easier reading
# Because a "list" is not a table (and therefore not easy to "write" for a computer)
# we are going to use a function called sink() to "capture" the output of resultsList and save it. Like this:
sink(file="resultsList_shannon_all_predictors.txt") # Opens file up to save
resultsList # Captures this output
sink() # closes the file and saves it


######## See if you can figure out beta diversity yourself! ########
# reminders:
# for beta diversity, the response variable will be your distance matrix
# You will need to create the distance matrix with your phyloseq object like this:
dm <- distance(phyloseqObject, method="bray")
# The distance matrix will be the same for all tests done by the loop, so this shoud go on the outside of the loop

# Then, you will use adonis2() (from the vegan package) instead of lm() 
permanova_results <- adonis2(dm ~ elevation, data=sampdat_filt, permutations = 10000) 
permanova_results
# You will want to save these outputs, just as you saved alpha diversity
# Note: I suggest increasing your "permutations" to 10,000 to allow for lower p-values.
# The default permutation size is 1000, which only allows for pvalues of 0.001 or more (because, by definition, your sample can only be compared to 1000 permutations (1/1000))
# The run time is slightly longer, but your results will be easier to interpret. 

######## stepAIC for alpha diversity ############

# Make a model with ALL predictors
# Here, the '.' is a "wildcard" symbol that means "use all columns in my data"
model_full <- lm(shannon_vec ~ ., data=sampdat_filt)
# Now, use stepAIC from the MASS package:
stepaic_results <- stepAIC(model_full)
stepaic_results # This is the best fit model

# You can compare AIC values like this:
AIC(stepaic_results) # This is your best fit model AIC
AIC(model_full) # This was your full model (all predictors)
AIC(lm(shannon_vec ~ averagesoilrelativehumidity, data=sampdat_filt )) 
# In analysis above, averagesoilrelativehumidity was the "strongest" predictor so let's compare just averagesoilrelativehumidity to our final model

# An AIC difference of approximately -2 means it is "meaningfully different"
# Our final model is MUCH better than a model with only our best predictor, averagesoilrelativehumidity
# AIC difference is 17.50188-40.97226 = -23

#### Our final model keeps:
# elevation
# transectname
# toc
# averavesoilrealtivehumidity
# averagesoiltemperature
# vegetation
#### But drops:
# depth
# ph

