dat <- read.delim("0_raw_data/metadata_raw.txt")

CV <- dat %>%
  select(Type, RichLevel, DateStart, Rep, Inhibitory,  ave_log10Copies_qPCR, sd_withinsample_qPCR) %>%
  filter(Type=="CV") %>%
  select(-Type) %>% rename(BacterialCopyNumber_mean=ave_log10Copies_qPCR, BacterialCopyNumber_sd=sd_withinsample_qPCR)
BD <- dat %>%
  select(Type, RichLevel, DateStart, Rep, Inhibitory, JarID, ave_log10Copies_qPCR, sd_withinsample_qPCR) %>%
  filter(Type=="Bd") %>%
  select(-Type) %>% rename(BdCopyNumber_mean=ave_log10Copies_qPCR, BdCopyNumber_sd=sd_withinsample_qPCR)

joined <- full_join(CV, BD)

write.table(joined, quote = FALSE, row.names = FALSE, sep="\t", file = "0_raw_data/metadata_simplified.txt")
