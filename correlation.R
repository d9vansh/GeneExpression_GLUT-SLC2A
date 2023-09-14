#/mnt/swift/devansh/SLC2

library(tidyverse)
library(reshape2)
library(Hmisc)

data <- working_tcga_data

wdata <- data %>% select(SLC2A1, SLC2A3, SLC2A4, SLC2A5, SLC2A6, SLC2A8, SLC2A9, SLC2A10, SLC2A11, SLC2A12, SLC2A13, SLC2A14)

wdata <- apply(wdata, 2, as.numeric)


#plotting correlation matrix--------------------------
cormat <- round(cor(wdata),2)
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(cormat)

melted_cormat <- melt(upper_tri, na.rm = TRUE)

#plotting crude heatmap
heatmap_plot <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 10, hjust = 1),
        axis.text.y = element_text(size = 10))+
  coord_fixed()

heatmap_plot

#calculate p values for heatmap----------------------------
p = rcorr(as.matrix(wdata[sapply(wdata,is.numeric)]))$P

get_upper_tri(p) -> p_upper

melted_p <- melt(p_upper, na.rm = T)
melted_p$value <- round(melted_p$value,3)
#calculating asterisks
stars.pval <- function(p.value)
{
  unclass(
    symnum(p.value, corr = FALSE, na = FALSE,
           cutpoints = c(0, 0.001, 0.01, 0.05, 1),
           symbols = c("***", "**", "*", "NS"))
  )
}
stars.pval(melted_p$value) -> melted_p$asterisk

# adding correlation scores and p value asterisks
final_heatmap <- heatmap_plot + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 3, vjust = 0) +
  geom_text(data = melted_p, aes(Var2, Var1, label = asterisk), size=3, vjust = 1.5) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
plot(final_heatmap)

ggsave('correlationplotd27.png',plot = final_heatmap, width = 345, height = 189, units = "mm", dpi = 300)