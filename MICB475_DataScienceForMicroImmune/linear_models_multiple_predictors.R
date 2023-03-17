library(tidyverse)
library(phyloseq)

### Let's use the fecal dataset as practise

#### Read in files ####
sampdatFP <- "fecal_sample-metadata.txt"
sampdat <- read.delim(sampdatFP)
otuFP <- "fecal_feature-table.txt"
otu <- read.delim(otuFP, skip=1)
taxaFP <- "fecal_taxonomy.txt"
taxa <- read.delim(taxaFP)

#### Format files to be a phyloseq object ####
### Sample data ###
sampdat_phylo <- sampdat[,-1]
rownames(sampdat_phylo) <- sampdat[,1]
SAMP <- sample_data(sampdat_phylo)
### OTU table ###
otu_phylo <- otu[,-1]
rownames(otu_phylo) <- otu[,1]
OTU <- otu_table(otu_phylo, taxa_are_rows = TRUE)
### Taxonomy ###
taxa_phylo <- taxa %>% 
  separate(Taxon, sep="; ", into=c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  select(-Consensus, -Feature.ID) %>% as.matrix()
rownames(taxa_phylo) <- taxa$Feature.ID
TAXONOMY <- tax_table(taxa_phylo)

### Make phyloseq object ###
fecal <- phyloseq(SAMP, OTU, TAXONOMY)
# Export sample metadata as a data frame
sampdat <- data.frame(sample_data(fecal))
# Export alpha diversity
alphadiv <- estimate_richness(fecal)

sampdat_wdiv <- data.frame(sampdat, alphadiv)


### Let's look at column names
head(sampdat_wdiv)

# Let's say you are interested  in changes in diversity through time ('week')
# You also want to know differences in treatment.group
# And finally, account for effects of sex.


# In this dataset, let's filter out 'donor' because we are only interested in knowing
# treatment vs control
sampdat_filt <- sampdat_wdiv %>%
  filter(treatment.group !="donor")

######## Treatment ###########
# Let's say we're interested in knowing whether treatment.group affects diversity. 
# We can plot this relationship with a boxplot like so:
ggplot(sampdat_filt) +
  geom_boxplot(aes(x=treatment.group, y=Shannon))

# Now, you may want to conduct a statistical test to see whether this difference
# is significant
# To run a simple linear model with treatment.group as predictor, we do this:
lm_treatmentonly <- lm(Shannon ~ treatment.group, data=sampdat_filt)
summary(lm_treatmentonly)
# We see there is no effect ot treatment (relative to control, which is the intercept)

######### Treatment and time ###########
# We may also be interested in knowing whether diversity changes with time.
# You can plot it like this-- here, I am colouring in dots by treatment group.
ggplot(sampdat_filt, aes(x=week, y=Shannon, col=treatment.group)) +
  geom_point()
# We can also add a "linear model" to the plot to see the overall trend:
ggplot(sampdat_filt, aes(x=week, y=Shannon, col=treatment.group)) +
  geom_point()+
  geom_smooth(method="lm")
# (note, we can also facet instead of using colours if you prefer)
ggplot(sampdat_filt, aes(x=week, y=Shannon)) +
  geom_point()+
  geom_smooth(method="lm") +
  facet_grid(.~treatment.group)
# Interesting!! The relationship between time and diversity is different between the two treatment groups

# Let's add time as a CONTINUOUS predictor in our model
lm_treatmenttime <- lm(Shannon ~ treatment.group*week, data=sampdat_filt)
summary(lm_treatmenttime)
# In this model, we use treatment.group*time, which is equivalent to:
# treatment.group + time + treatment.group:time
#### treatment.group will tell you whether there is a difference in the INTERCEPTS of control and treatment
#### time will tell you whether slope for control is significantlly different from zero
#### treatment.group:time will tell you whether the slope (changes to diversity through time) is different between control and treatment

# What we see is that diversity changes with time- and that it changes differently between control and treatment


##### Treatment, time, and gender ##########
# Finally, let's say you want to account for possible effects of gender.
# Let's try a simple plot first:
ggplot(sampdat_filt, aes(x=treatment.group, y=Shannon, col=gender)) +
  geom_boxplot()
# It seems that men, on average, have more diversity than women.
# If we add gender to our previous plot, what does it look like?
ggplot(sampdat_filt, aes(x=week, y=Shannon, col=gender)) +
  geom_point()+
  geom_smooth(method="lm") +
  facet_grid(.~treatment.group)
# The signal is very faint and difficult to tell if the effect is significant. 
# So let's run a model. Note: here, I'm going to include only the effects of 
# gender. This is because I do not expect men and women to differ in
# how their gut microbiome changes through time, but it does look like men
# might have, on average, richer microbiota than women.
# This is an arbitrary decision-- there is no 'right' or 'wrong' answer
# in how you build a model, but you should be able to justify why you included
# the terms that you include. 
lm_treatmenttimegender <- lm(Shannon ~ treatment.group*week +gender, data=sampdat_filt)
summary(lm_treatmenttimegender)
### Interpretation is as follows:
# Intercept-- Control y-intercept is different from zero
# treatment.grouptreatment-- The y-intercept for the treatmentline is smaller (-0.79625) than the y-intercept for control
# week-- Diversity does not change significantly through time for CONTROL individuals (p=0.06)
# genderm-- Diversity does not differ significantly between men and women (after accounting for all other factors)
# treatment.grouptreatment:week-- There is a significant difference in the way that treatment individuals change through time compared to control
# Note that treatment.grouptreatment:week has a positive estimate (0.10229). This means the slope for treatment is 0.10229 steeper than control.

