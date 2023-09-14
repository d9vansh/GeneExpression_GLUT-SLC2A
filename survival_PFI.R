#/mnt/swift/devansh/SLC2

library(survminer)
library(survival)
library(tidyverse)
library(dplyr)

data <- working.data
cdata <- data

#finding median
data$SLC2A2 <- as.numeric(data$SLC2A2)
med <- median(data$SLC2A2)

cdata$SLC2A2 <- as.numeric(cdata$SLC2A2)

counter <- 1
for (i in cdata$SLC2A2){
  if (i < med){
    i = '1'
  }
  else if (i >= med){
    i = '2'
  }
  cdata$SLC2A2[counter] <- i
  counter <- counter + 1
}
cdata$SLC2A2 <- as.integer(cdata$SLC2A2)
typeof(cdata$SLC2A2)
sum(cdata$SLC2A2 == 1)

fit <- survfit(Surv(PFI.time, PFI) ~ SLC2A2, data = cdata)
ggsurvplot(fit,
           pval = TRUE)