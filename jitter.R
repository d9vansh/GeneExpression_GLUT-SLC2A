#/mnt/swift/devansh/SLC2

library(tidyverse)
library(ggplot2)

is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

calc_boxplot_stat <- function(x) {
  coef <- 1.5
  n <- sum(!is.na(x))
  # calculate quantiles
  stats <- quantile(x, probs = c(0.0, 0.25, 0.5, 0.75, 1.0))
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
  iqr <- diff(stats[c(2, 4)])
  # set whiskers
  outliers <- x < (stats[2] - coef * iqr) | x > (stats[4] + coef * iqr)
  if (any(outliers)) {
    stats[c(1, 5)] <- range(c(stats[2:4], x[!outliers]), na.rm = TRUE)
  }
  return(stats)
}

data1 <- working_tcga_data
data <- data1[,!names(data1) %in% c("SLC2A2", "SLC2A7")]

gathered <- data %>% gather(key = genename, value = expression.level, SLC2A1:SLC2A14)
typeof(gathered$expression.level)
gathered$expression.level <- as.double(gathered$expression.level)

gathered %>% mutate(outlier = is_outlier(expression.level)) %>% filter(outlier == FALSE) -> ooutfilgathered

typeof(ooutfilgathered$genename)
ooutfilgathered$genename <- as.factor(ooutfilgathered$genename)
ooutfilgathered$genename <- factor(ooutfilgathered$genename, c("SLC2A1", "SLC2A3", "SLC2A4", "SLC2A5", "SLC2A6", "SLC2A8", "SLC2A9", "SLC2A10", "SLC2A11", "SLC2A12", "SLC2A13", "SLC2A14"))

plot1 <- ggplot(ooutfilgathered, aes(x=genename, y=expression.level, color=genename)) + 
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", width=0.5) + geom_jitter(alpha=0.6, width = 0.3) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base=10))  +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, colour = "black")) + 
  labs(title="HNSCC", x ="Genes", y = "Expression Level")
plot1

ggsave('jitter.png',plot = plot1, width = 317, height = 165, units = "mm", dpi = 300)