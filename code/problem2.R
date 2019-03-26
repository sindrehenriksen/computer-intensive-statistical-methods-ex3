setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

## --- labsdata2
library(tidyverse)
library(plotly)
bili.df <- as.data.frame(read.table("../data/bilirubin.txt",header=T))

## --- boxplot.bili
boxplot.bili <- ggplot(bili.df, aes(x= pers, y = log(meas),fill=pers)) + 
  geom_boxplot()+
  ggtitle("boxplot of bilirubin data")

## --- break
ggsave("../figures/boxplot_bili.pdf", plot = boxplot.bili, device = NULL, 
       path = NULL,scale = 1, width = 5.5, height = 2*4, 
       units = "in",dpi = 300, limitsize = TRUE)

## --- linearReg2
log.y <- lm(log(meas)~pers,data=bili.df)
sum.log.y <- summary(log.y)
Fval <- sum.log.y$fstatistic
