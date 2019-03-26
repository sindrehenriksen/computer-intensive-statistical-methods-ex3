setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list = ls())

## ---- labsdata2
library(tidyverse)
bili.df <- as.data.frame(read.table("../data/bilirubin.txt",header=T))

## ---- boxplot.bili
# boxplot of data
boxplot.bili <- ggplot(bili.df, aes(x= pers, y = log(meas),fill=pers)) + 
  geom_boxplot()+
  ggtitle("boxplot of bilirubin data")

boxplot.bili
## ---- break
ggsave("../figures/boxplot_bili.pdf", plot = boxplot.bili, device = NULL, 
       path = NULL,scale = 1, width = 5.5, height = 2*4, 
       units = "in",dpi = 300, limitsize = TRUE)

## ---- linearReg2
# linear regression of the data
log.y <- lm(log(meas)~pers,data=bili.df)
# getting the summary
sum.log.y <- summary(log.y)
# fetching f.statistics
Fval <- as.numeric(sum.log.y$fstatistic[1])

## ---- permTest
permTest <- function(bili.df){
  # temporary dataframe
  bili.p.df = bili.df
  # permutation on the grouping column
  bili.p.df$pers = bili.df$pers[sample(1:nrow(bili.p.df))]
  # return of summary.lm and fetching the fstat
  return(as.numeric(summary(lm(log(meas)~pers,
                               data = bili.p.df))$fstatistic[1]))
  }

## ---- permutationTest
N <- 999
Fvals <- vector()
# running the 999 permutions in permTest
for (i in seq(1,N)){
  Fvals = c(Fvals, permTest(bili.df))
}
# dataframe for ggplot histogram
fvals.df <- enframe(Fvals)

# histogram of fstat of permutations
fstat.plot<-ggplot(fvals.df)+
  geom_histogram(aes(x = value, y = ..density..),bins= 25)+
  geom_vline(xintercept = Fval, color = "red")+
  ggtitle("F-statistics of permutations")

fstat.plot

# calculation of p-value
p.F <- mean(Fvals >= Fval)
cat("p-value of F-statistics:", p.F)
## ---- break
ggsave("../figures/his_fstat.pdf", plot = fstat.plot, device = NULL, 
       path = NULL,scale = 1, width = 5.5, height = 5.5, 
       units = "in",dpi = 300, limitsize = TRUE)
save(file = "../data/variables/p_val_fstat.Rdata", p.F)

## ---- printPval
load(file = "../data/variables/p_val_fstat.Rdata")
cat("p-value of F-statistics:", p.F)
