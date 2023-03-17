#!bin/bash

# Combine all Data
library(tidyverse)

#### Pathways ####
bd_PW <- "1_process_counts/downstream/dat_count_agg.txt"
meta_PW <- "0_raw_data/metadata_simplified.txt"

dir.create("2_combine_data")
### Load ####
bd_dat <- read.delim(bd_PW)
meta_dat <- read.delim(meta_PW)

# Check it loaded correctly
head(bd_dat) # head() looks at first few lines
head(meta_dat)

# Combine using tidyverse's "full_join"
# You can use "left_join" if you want to use the first data frame as the "template"; 
# "right_join" if you want to use the second data framea s the "template", 
# or "inner_join" if you only want to keep rows and columns that are common between both datasets.

all_dat <- full_join(bd_dat, meta_dat) %>%
  unite(Inhibitory, RichLevel, Rep, col = SampleID, sep="_", remove = FALSE) # set remove=FALSE so we don't lose those columns

dir.create("2_combine_data/downstream")
write.table(all_dat, quote=FALSE, row.names = FALSE, sep="\t", file = "2_combine_data/downstream/full_metadata.txt")
