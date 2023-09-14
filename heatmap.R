#/mnt/swift/devansh/SLC2

library(dplyr)
library(pheatmap)
library(gridExtra)
library(scales)
library(patchwork)
library(tidyverse)
library(tidyr)
library(ggsci)

wdata <- fdata
rownames(wdata) <- wdata[,1]

#colnames(wdata) <- wdata[1,]

wdata %>% dplyr::select(Group, hpv.status, SLC2A1, SLC2A3, SLC2A4, SLC2A5, SLC2A6, SLC2A8, SLC2A9, SLC2A10, SLC2A11, SLC2A12, SLC2A13, SLC2A14, tobacco_smoking_history, gender) -> hdata
hdata %>% arrange(Group) -> hdataarr

#BA barcodes
hdataarr %>% filter(Group == "BA") -> hdataba
FOOba <- as.matrix(hdataba)
FOOba <- dist(FOOba)
FOOba <- hclust(FOOba)
crap_vec <- c()
for (variable in FOOba$order) {
  crap_vec <- c(crap_vec, rownames(hdataba)[variable])
}
crap_vecba <- crap_vec

#CL barcodes
hdataarr %>% filter(Group == "CL") -> hdatacl
FOOcl <- as.matrix(hdatacl)
FOOcl <- dist(FOOcl)
FOOcl <- hclust(FOOcl)
crap_vec <- c()
for (variable in FOOcl$order) {
  crap_vec <- c(crap_vec, rownames(hdatacl)[variable])
}
crap_veccl <- crap_vec

#MS barcodes
hdataarr %>% filter(Group == "MS") -> hdatams
FOOms <- as.matrix(hdatams)
FOOms <- dist(FOOms)
FOOms <- hclust(FOOms)
crap_vec <- c()
for (variable in FOOms$order) {
  crap_vec <- c(crap_vec, rownames(hdatams)[variable])
}
crap_vecms <- crap_vec

sample_order <- c(crap_vecba, crap_veccl, crap_vecms)

#clustering Y-axis genes
#FOOgene <- as.matrix(hdataarr)
#FOOgene <- dist(FOOgene)
#FOOgene <- hclust(FOOgene)
#crap_vec <- c()
#for (variable in FOOgene$order) {
#  crap_vec <- c(crap_vec, colnames(hdataarr)[variable])
#}
#crap_vecgene <- crap_vec

hdata$barcode <- rownames(hdata)

hdata[match(sample_order, hdata$barcode, hdata$hpv.status),] -> hdatasorted
hdatasorted %>% gather(key = 'Gene', value = 'expression', SLC2A1:SLC2A14) -> ghdata

#ghdata$Gene <- factor(ghdata$Gene, levels = c("SLC2A5", "SLC2A7", "SLC2A14", "SLC2A11", "SLC2A6", "SLC2A13", "SLC2A2", "SLC2A4", "SLC2A9", "SLC2A12", "SLC2A3", "SLC2A1", "SLC2A10", "SLC2A8"))
ghdata$Gene <- factor(ghdata$Gene, levels = c("SLC2A1", "SLC2A3", "SLC2A4", "SLC2A5", "SLC2A6", "SLC2A8", "SLC2A9", "SLC2A10", "SLC2A11", "SLC2A12", "SLC2A13", "SLC2A14"))
ghdata$barcode <- factor(ghdata$barcode, levels = sample_order)
ghdata$hpv.status <- factor(ghdata$hpv.status, levels = hpv.status)

typeof(ghdata$expression)
ghdata$expression <- as.double(ghdata$expression)

ghdata$Expression <- ">5.0"
ghdata$Expression[ghdata$expression <= 5.0] <- "NO"
  
slc.plot <- ggplot() + 
  geom_tile(data = subset(ghdata, expression <= 5.0), 
            aes(x=barcode, y=Gene, fill=expression)) +
  scale_fill_gradient2(high = "red",mid = "white", low = "blue") +
  geom_tile(data = subset(ghdata, expression > 5.0), 
            aes(x=barcode, y=Gene, color=Expression)) +
  scale_colour_manual(name = "Expression", values = c('black')) +
  theme(text = element_text(size = 10), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

slc.plot

ghdata$dummy <- 'value'

group.plot <- ggplot() +
  geom_tile(data = ghdata, aes(x=barcode, y=dummy, fill=Group)) +
  scale_fill_uchicago() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

hpv.plot <- ggplot() +
  geom_tile(data = ghdata, aes(x=barcode, y=dummy, fill=hpv.status)) +
  scale_fill_jco() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

smoking.plot <- ggplot() +
  geom_tile(data = ghdata, aes(x=barcode, y=dummy, fill=tobacco_smoking_history)) +
  scale_fill_tron() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

gender.plot <- ggplot() +
  geom_tile(data = ghdata, aes(x=barcode, y=dummy, fill=gender)) +
  scale_fill_simpsons() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

combined.plot <- group.plot+hpv.plot+smoking.plot+gender.plot+slc.plot+plot_layout(ncol = 1, heights = c(0.25, 0.25, 0.25, 0.25, 6))
combined.plot

ggsave('group.plot.png',plot = group.plot, width = 317, height = (165/6), units = "mm", dpi = 300)
ggsave('hpv.status.plot.png',plot = hpv.plot, width = 317, height = (165/6), units = "mm", dpi = 300)
ggsave('slc.plotd.png',plot = slc.plot,width = 317, height = 165, units = "mm", dpi = 300)
ggsave('combined.plot2.png',plot = combined.plot, width = 317, height = 165, units = "mm", dpi = 300)
ggsave('smoking.plot.png',plot = smoking.plot, width = 317, height = (165/6), units = "mm", dpi = 300)
ggsave('gender.plot.png',plot = gender.plot, width = 317, height = (165/6), units = "mm", dpi = 300)

ggsave('heatmap.png',plot = slc.plot, width = 317, height = 165, units = "mm", dpi = 300)