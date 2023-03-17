#!bin/bash Rscript

#### Load data ####

dat <- read.csv(file.choose())

#### Select one isolate ####


#### Plot out ####

plot(dat$time_h, dat$OD)
