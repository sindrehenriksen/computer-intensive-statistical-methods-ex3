#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#rm(list = ls())
## ---- libsData3
library(tidyverse)

data.df <- data.frame(
  z = read.table("../data/z.txt",header=F,col.names = "z"),
  u = read.table("../data/u.txt",header=F,col.names = "u")
  )


