#/mnt/swift/devansh/SLC2

library(tidyverse)

data <- working_tcga_data
dup.data <- data[,!names(data) %in% c("SLC2A2", "SLC2A7")]
dup.data %>% gather(key = 'gene.name', value = "expression.level", SLC2A1:SLC2A14) -> gathered
gathered$gene.name <- as.character(gathered$gene.name)
gathered$gene.name <- factor(gathered$gene.name, levels = c("SLC2A1", "SLC2A2", "SLC2A3","SLC2A4", "SLC2A5", "SLC2A6", "SLC2A7", "SLC2A8", "SLC2A9", "SLC2A10", "SLC2A11", "SLC2A12", "SLC2A13", "SLC2A14"))

typeof(gathered$expression.level)
gathered$expression.level <- as.double(gathered$expression.level)

po <- position_dodge(width = 0.8)

plot1 <-ggplot(gathered, aes(x=gene.name, y=expression.level, fill = Group)) + 
  stat_boxplot(geom = "errorbar", width = 0.4, position = position_dodge(width = 0.75)) + 
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim =  c(-1.5, 3.5)) +
  scale_fill_manual(values = c("#A1A1A1","#FFFFFF","#555555")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, colour = "black"), axis.text.y = element_text(colour = "black"), axis.title.x = element_text(colour = "black", size = 12), axis.title.y = element_text(colour = "black", size=12), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  labs(title="HNSCC", x ="Genes", y = "Expression Level")
plot1

stars.pval <- function(p.value)
{
  unclass(
    symnum(p.value, corr = FALSE, na = FALSE,
           cutpoints = c(0, 0.001, 0.01, 0.05, 1),
           symbols = c("***", "**", "*", "NS"))
  )
}

#p-values for Group
anova.pval <- c()
for (i in c(1, 3:6, 8:14)) {
  res.aov <- aov(expression.level ~ Group, data = gathered %>% filter(gene.name == paste('SLC2A',i,sep = "")))
  anova.pval <- c(anova.pval, round((summary(res.aov)[[1]][["Pr(>F)"]][[1]]),10))
}
stars.pval(anova.pval) -> anova.dot

#adding asterisks to the plot
plot1 + annotate("text", x=c(1), y=c(3.5), label=anova.dot[1], size=4.5) -> g2
for (i in c(2:12)) {
  g2 + annotate("text", x=c(i), y= 3.5, label=anova.dot[i], size=4.5) -> g2
}

#adding ___ to the plot
g2 + annotate("text", x=c(1), y=c(3.5), label="____", size=4.5) -> g2
for (i in c(2:12)) {
  g2 + annotate("text", x=c(i), y=c(3.5), label="____", size=4.5) -> g2
}
g2

ggsave('boxplot_group_FINALd.png',plot = g2, width = 345, height = 189, units = "mm", dpi = 300)