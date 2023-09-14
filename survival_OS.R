#/mnt/swift/devansh/SLC2

library(survminer)
library(survival)
library(tidyverse)
library(dplyr)

data <- working.data
cdata <- data

#finding median
data$SLC2A14 <- as.numeric(data$SLC2A14)
med <- median(data$SLC2A14)

cdata$SLC2A14 <- as.numeric(cdata$SLC2A14)

counter <- 1
for (i in cdata$SLC2A14){
  if (i < med){
    i = '1'
  }
  else if (i >= med){
    i = '2'
  }
  cdata$SLC2A14[counter] <- i
  counter <- counter + 1
}
cdata$SLC2A14 <- as.integer(cdata$SLC2A14)
typeof(cdata$SLC2A14)
sum(cdata$SLC2A14 == 1)

fit <- survfit(Surv(OS.time, OS) ~ SLC2A14, data = cdata)
ggsurvplot(fit,
           pval = TRUE)