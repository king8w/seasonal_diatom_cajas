## Perform richness calculation of phytoplankton data

phytoplankton <- read.csv("data/phytoplankton.csv", sep=";") %>%
  filter(division!="Bacillariophyta") %>%
  select(-c("division", "class", "order", "family")) 

phytoplankton <- as.data.frame(t(phytoplankton))
str(phytoplankton)

phytoplankton_raw <- read.csv("data/phytoplankton.csv", sep=";") %>%
  filter(division!="Bacillariophyta") %>%
  select(-c("division", "class", "order", "family")) 

df <- mutate_all(phytoplankton, function(x) as.numeric(as.character(x)))
colnames(df) <- phytoplankton_raw$genera
df <- df[-1,]

write.csv(df, "data/phytoplankton_nondiatoms.csv")

# calculate cel/L mean across lake depths
n <-3 # there are 3 depth observations
df_mean <- aggregate(df, list(rep(1:(nrow(df) %/% n + 1), each = n, len = nrow(df))), mean)[-1]
row.names(df_mean) <- c("CJ-001", "CJ-003", "CJ-020", "CJ-023", "CJ-024", "CJ-039", "CJ-043", "CJ-051", "CJ-097", "CJ-121")

# calculate phytoplankton (non-diatoms) richness per lake
richness <- apply(df_mean>0,1,sum)
richness_seq <- rep(richness, 4)
nms <- names(richness_seq)
richness_seq_df <- as.data.frame(richness_seq) # this is the ordered vector of phytoplankton richness
write.csv(richness_seq_df, "data/phyto_richness.csv", row.names = TRUE)


