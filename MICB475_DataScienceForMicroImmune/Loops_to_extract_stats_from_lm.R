library(tidyverse)
library(phyloseq)

#### Code for extracting statistics #######


# Let's say you have a phyloseq object called "phyloseqObject"
# Here, I'm going to use the atacama dataset as an example. I haven't included 
# the code here for it, so you should replace atacama_final with whatever 
# your phyloseq object is called.
load("atacam_final.RData")
phyloseqObject <- atacama_final

######## Setting up loop ########
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

### This is the code for the previous loop I wrote for you:
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


######### Stats extraction ###########
# Now, instead of saving each model result into a list, 
# let's say we want to make a table
# To do this, we need to make a DATA FRAME
results_in_table <- data.frame()

## Now, within a loop (you can actually integrate the code below 
# into the original loop above)

# Let's say you run a linear model with the predictor 'x'
# We're going to arbitrarily assign 'x' as elevation
 x<- "elevation"
lm_example <- lm(shannon_vec ~ get(x), data=sampdat)
# To extract values, you can use summary to produce a table of coefficients
summary(lm_example)$coefficients
# Then, you can "insert" these values into the data frame. 
# Note! Make sure you make an extra column that tells you which predictor you're using
# Since you used get(x), you won't be able to tell which lines are from
# which predictors after the table is all merged
new_results <- summary(lm_example)$coefficients %>%
  as.data.frame() %>% rownames_to_column(var="CoefficientName") %>% # this makes the rownames a column
  mutate(predictor = x) # This tells you which predictor you are using

# Now, you can "rbind" (row-bind) your new results with your data frame
# and save the results BACK into the 'results_in_table' object 
results_in_table <- rbind(results_in_table, new_results)

### Note: the next time your code goes through the loop, you'll add
# on the additional info. For example:
x <- "ph" # re-set x to ph this time
lm_example <- lm(shannon_vec ~ get(x), data=sampdat)
new_results <- summary(lm_example)$coefficients %>%
  as.data.frame() %>% rownames_to_column(var="CoefficientName") %>% # this makes the rownames a column
  mutate(predictor = x) # This tells you which predictor you are using
results_in_table <- rbind(results_in_table, new_results)
results_in_table # View results here; now there are TWO sets of results


######### Plotting results ##########
# After you loop through all variables, you want to:
results_in_table_filt <- results_in_table %>%
  filter(CoefficientName != "(Intercept)") %>% # Filters out intercept coefficient
  rowwise() %>% mutate(CoefficientName = gsub("get(x)", predictor, CoefficientName, fixed=TRUE)) %>% # substitutes (gsub) get(x) with the predictor name, by row
  ungroup() # This just gets rid of the "rowwise" designation, which groups data by rows

# You can now either print the table as-is:
results_in_table_filt
# Save it like this:
write.table(results_in_table_filt, file="results_table.txt", quote = FALSE, sep="\t", row.names = FALSE)
# and open in excel


# Or make a plot like this:
ggplot(results_in_table_filt, aes(x=CoefficientName, y=Estimate, col=`Pr(>|t|)`<0.05)) +
  geom_point() +
  geom_errorbar(aes(ymin=Estimate-`Std. Error`, ymax= Estimate+`Std. Error`))+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) 


######## PERMANOVA stat extraction ###########
#### For extracting stats from permanovas, you would do this:
dm <- distance(phyloseqObject, method="bray")
x <- "elevation"
permanova_results <- adonis2(dm ~ get(x), data=sampdat_filt, permutations = 10000) 
new_results <- permanova_results %>% as.data.frame() %>%
  rownames_to_column(var="CoefficientName") %>%
  mutate(predictor=x)
new_results


# And after merging all tables, you should filter out "Residual" and "Total"
# from column 'CoefficientName'. Adjust the code above, where you remove '(Intercept)'
