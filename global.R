Lang <- read.csv("changes.csv")$x
#if (is.na(Lang)) Lang = "EN"
print(paste0("layout ", Lang, ".R"))
source(paste0("layout ", Lang, ".R"))
