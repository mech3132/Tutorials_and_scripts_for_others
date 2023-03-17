#!bin/bash

### Processing Bd microscopy counts and measurements ####
library(tidyverse)

#### Pathways ####
bd_count_dat <- "0_raw_data/Bd_count_data_AW.txt"

#### Load ####
dat_count <- read.delim(bd_count_dat)

#### Output folder ####
dir.create("1_process_counts")

## Count data processing ##
colnames(dat_count) # Look at column names, just to make sure we uploaded data correctly
# Look at problematic ones first
dat_count[which(dat_count$too_blurry_to_count==TRUE),c("Rep","Jar","Picture","Count","Notes")] # %>% View() # View makes it pop up in a nice window
# After reviewing all of them, I think we should just omit these counts. They are too difficult to count.
# I filter these images out below.


### QUADRANTS #### 
# Basically, we need to record which quadrants we measured each time so we have Zoosporangia/area.
# This is a very long, confusing section because we need to adjust all entries to either be TRUE or FALSE
# "NA" is treated very differently in R-- they are neither TRUE nor FALSE, so filtering by TRUE/FALSE won't work

# Key for new quadrant designation:
# 1=TopLeft, 2=TopRight, 3=BottomLeft, 4=BottomRIght, 5=Top, 6=Bottom 7=Left, 8=Right, 9=Diagonla, 10=FULL

# Notice how below, each line is a "new" command, separated by the pipe (%>%) at the end.
# This makes it easier to read (for future you, and for collaborators), but also bulks up your code.
dat_count_cleaned <- dat_count %>% 
  mutate(too_blurry_to_count=ifelse(is.na(too_blurry_to_count), FALSE, too_blurry_to_count)) %>% #### Blurry counts ####
  filter(!too_blurry_to_count) %>% # Filter out all "NOT" "too_blurry_to_count"
  mutate(pinpricks_present = ifelse(pinpricks_present%in% c("TRUE","TRUE?") , TRUE, FALSE)) %>%  #### Pinpricks ####
  mutate(clouds_present = ifelse(is.na(clouds_present), FALSE, clouds_present))%>%   #### Clouds ####
  mutate(very_small = ifelse(very_small %in% c("TRUE","Possibly","TRUE?"), TRUE, FALSE)) %>%   #### very small ####
  mutate(Count = ifelse(Count=="??", NA, as.numeric(Count))) %>%   #### very small ####
  rowwise() %>%  # Make calculations by row 
  #### Fixing quadrants to be T/F, no NAs ####
  mutate(Full = !any(!is.na(c(TopLeft, TopRight, BottomLeft, BottomRight)))
         , Top = sum(c(TopLeft, TopRight, is.na(c(BottomLeft, BottomRight))))==4
         , Bottom = sum(c(BottomLeft, BottomRight, is.na(c(TopLeft, TopRight ))))==4
         , Left = sum(c(BottomLeft, TopLeft, is.na(c(BottomRight, TopRight ))))==4
         , Right = sum(c(BottomRight, TopRight, is.na(c(BottomLeft, TopLeft ))))==4
         , OnlyTopLeft = sum(c(TopLeft, is.na(c(BottomLeft, TopRight, BottomRight ))))==4
         , OnlyTopRight = sum(c(TopRight, is.na(c(BottomLeft, TopLeft, BottomRight ))))==4
         , OnlyBottomLeft = sum(c(BottomLeft, is.na(c(TopLeft, TopRight, BottomRight ))))==4
         , OnlyBottomRight = sum(c(BottomRight, is.na(c(BottomLeft, TopRight, TopLeft ))))==4
         , Diagonal = (sum(c(BottomRight, TopLeft, is.na(c(TopRight, BottomLeft))))==4 |sum(c(TopRight, BottomLeft, is.na(c(BottomRight, TopLeft))))==4)
         )  %>%
  mutate(Top = ifelse(is.na(Top), FALSE, Top)
         , Bottom = ifelse(is.na(Bottom), FALSE, Bottom)
         , Left = ifelse(is.na(Left), FALSE, Left)
         , Right = ifelse(is.na(Right), FALSE, Right)
         , OnlyTopLeft = ifelse(is.na(OnlyTopLeft), FALSE, OnlyTopLeft)
         , OnlyTopRight = ifelse(is.na(OnlyTopRight), FALSE, OnlyTopRight)
         , OnlyBottomLeft = ifelse(is.na(OnlyBottomLeft), FALSE, OnlyBottomLeft)
         , OnlyBottomRight = ifelse(is.na(OnlyBottomRight), FALSE, OnlyBottomRight)
         , Diagonal = ifelse(is.na(Diagonal), FALSE, Diagonal)
  ) %>%
  mutate(image_quadrant = ifelse(Full, 10, ifelse(OnlyTopLeft,1,ifelse(OnlyTopRight,2,ifelse(OnlyBottomLeft,3,ifelse(OnlyBottomRight,4,ifelse(Top,5,ifelse(Bottom,6,ifelse(Left,7,ifelse(Right,8,ifelse(Diagonal,9,NA))))))))))) %>%   
  #### Manual area in put ####
  mutate(Area_whole = 0.16651143, Area_unit = "mm2"
         , Area_counted = ifelse(image_quadrant %in% c(1,2,3,4), Area_whole/4, ifelse(image_quadrant %in% c(5,6,7,8,9), Area_whole/2, Area_whole))) %>%
  select(Rep, Jar, Picture, Count, very_small, pinpricks_present, clouds_present, image_quadrant, Area_counted, Area_unit, Notes) %>% #### remove unwanted cols####
  filter(!is.na(Rep)) %>%
  mutate(zoosp_dens_perimage = Count/Area_counted
    , Rep = paste0("Rep",Rep))%>%
  rename(JarID = Jar)

# Aggregate counts by Sample
dat_count_agg <- dat_count_cleaned %>%
  group_by(Rep, JarID) %>% # Summarize by unique combinations of these columns
  summarize(Bd_density = sum(Count, na.rm=TRUE)/sum(Area_counted)
            , Bd_sd = sd(zoosp_dens_perimage, na.rm=TRUE)
            , pinpricks_present = any(pinpricks_present)
            , clouds_present = any(clouds_present))
### After doing a bunch of data manipulation, always plot it out to see if you've made a grave mistake.
## Plotting variability of counts
dat_count_cleaned %>%
  ggplot()+geom_point(aes(x=as.factor(JarID), y=zoosp_dens_perimage, col=pinpricks_present)) +
  facet_wrap(.~Rep, nrow=3)


### Write out for downstream
dir.create("1_process_counts/downstream")
write.table(dat_count_cleaned, "1_process_counts/downstream/dat_count_clean.txt", quote=FALSE, row.names = FALSE,sep="\t")
write.table(dat_count_agg, "1_process_counts/downstream/dat_count_agg.txt", quote=FALSE, row.names = FALSE,sep="\t")

# Note: you can store output in text files or .RData files (using save())-- I like text files because they're easy to read outside of R, but if it's a super large intermediate file I'd just go with RData

