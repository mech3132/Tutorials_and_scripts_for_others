#!bin/bash
library("tidyverse")
library("growthcurver")

dat <- read.csv(file="01_combined_dat/01_sample_dat.csv")

dat$time_original

dat_adj <- dat %>% separate(col = time_original, into=c("dd","rest"), sep = "[.]", fill="left") %>% 
  separate(col=rest,into=c("hh","mm","ss"),sep="[:]", fill="left") %>%
  mutate(dd = ifelse(is.na(dd), 0, dd)) %>% 
  mutate(dd = as.numeric(dd), hh = as.numeric(hh), mm = as.numeric(mm), ss = as.numeric(ss))%>% 
  mutate(time_h = dd*24 + hh + mm/60 + ss/360) 

#subtract the blank then fit growth curve

str(dat)

dat$isolate

dat[c(1,10,20),]
#c is for combining variables

dat_a1 <- dat_adj[dat_adj$well == "A1",]
dat_a1 <- dat_adj[dat_adj$isolate == "ATD 4B4",]

# To save plot, surround plotting function with png() and dev.off()
png(filename = "Downloads/2022-05-11_growthcurve_aileen/01_combined_dat/a1.png")
plot(x=dat_a1$time_h, y=dat_a1$OD, type='l')
dev.off()

##### ---------------- Some other plotting tricks ---------------- ######

#### If you want to change the x and y axis labels, you can do it like so:
plot(x=dat_a1$time_h, y=dat_a1$OD, type='l', xlab = "Time (h)", ylab="OD")

##### If you want to plot dots AND lines, use type='b' (both) instead of 'l' (line)
plot(x=dat_a1$time_h, y=dat_a1$OD, type='b')


#dev.off closes the file so that it saves
#####

list_isolate <- unique(dat_adj$isolate)
list_isolate

for ( i in list_isolate ) {
  
  # Get subset of data where 'isolate' is equal to 'i', save as dat_i
  dat_i <- dat_adj[dat_adj$isolate == i,]
  
  # Create the file path you want to save your graph under.
  # The "paste0()" function takes character strings and "pastes" them together into ONE string.
  # Here, we are pasting together the file path you want, plus a file name that will change with every loop
  # *Change to YOUR file path on your computer here!*
  filepath_tosave <- paste0("Downloads/2022-05-11_growthcurve_aileen/01_combined_dat/", i, ".png") # See how the file name will change every loop? Since i changes every loop?
  
  #   png(filename = '01_combined_dat/dat_a1.png') # This was the old version-- see how we've replaced it with our string?
  png(filename = filepath_tosave)
  plot(x=dat_i$time_h, y=dat_i$OD, type='l',ylab="OD", xlab="Time (h)", col="blue")
  # Can you play around with this "plot" function to get the colour, shape, and axis labels you want?
  dev.off()
  
}


dat_c8 <- dat_adj[dat_adj$isolate == "ATT 4B4",]
unique(dat_c8$well.1)
# To save plot, surround plotting function with png() and dev.off()
png(filename = "Downloads/2022-05-11_growthcurve_aileen/01_combined_dat/ATT 4B4.png")
plot(x=dat_c8$time_h, y=dat_c8$OD, type='l', xlab = "Time (h)", ylab="OD", col="blue")
dev.off()

dat_f8 <- dat_adj[dat_adj$isolate == "ATL 4A1",]
unique(dat_f8$well.1)
# To save plot, surround plotting function with png() and dev.off()
png(filename = "Downloads/2022-05-11_growthcurve_aileen/01_combined_dat/ATL 4A1.png")
plot(x=dat_f8$time_h, y=dat_f8$OD, type='l', xlab = "Time (h)", ylab="OD", col="blue")
dev.off()

dat_g9 <- dat_adj[dat_adj$isolate == "ATT 3B5",]

# To save plot, surround plotting function with png() and dev.off()
png(filename = "Downloads/2022-05-11_growthcurve_aileen/01_combined_dat/ATT 3B5.png")
plot(x=dat_g9$time_h, y=dat_g9$OD, type='l', xlab = "Time (h)", ylab="OD", col="blue")
dev.off()


# =========================

unique(dat_adj$isolate)
c("BLANK B4","BLANK G6", "BLANK H12")

dat_blank <- dat_adj[dat_adj$isolate %in% c("BLANK B4","BLANK G6", "BLANK H12"),]

dat_blank_median <- dat_blank %>%
  group_by(time_h) %>%
  summarise(blankOD = median(OD))

plot(dat_blank_median$time_h, dat_blank_median$blankOD, type="l")


full_join(dat_adj, dat_blank_median) %>%
  mutate()

