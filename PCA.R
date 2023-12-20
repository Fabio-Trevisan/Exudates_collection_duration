library(devtools)
library(ggbiplot)
library(ggpubr)

library(ggplot2)
library(remotes)

table <- read.csv("Exudates_ISTD_Area_.csv", sep=",",
                  header=T)
Caption <- "Area"



#Groups definition or combination ####
nr. = ncol(table)
table$Time <- as.factor(table$Time)

# Time ####
H2 <- table[table$`Time` == "2",]
x2 <- H2[1,3]
Title2 <- paste("Exudates", x2, sep="_")
H4 <- table[table$`Time` == "4",]
x4 <- H4[1,3]
Title4 <- paste("Exudates", x4, sep="_")
H6 <- table[table$`Time` == "6",]
x6 <- H6[1,3]
Title6 <- paste("Exudates", x6, sep="_")
Title <- "Exudates_Time"

# Principal components
Exudates.2 <- prcomp(H2[,c(4:nr.)], center = TRUE, scale=TRUE, retx = T) #(change number of columns -> compounds identified)
Exudates.4 <- prcomp(H4[,c(4:nr.)], scale = TRUE) #(change number of columns -> compounds identified)
Exudates.6 <- prcomp(H6[,c(4:nr.)], scale = TRUE) #(change number of columns -> compounds identified)

# PCA
PCA_2 <- ggbiplot(Exudates.2, var.axes=FALSE, ellipse=TRUE, labels=H2$Sample, groups=H2$Tr, scale = 1) +
  theme_classic() + 
  ggtitle(Title2) + 
  theme(plot.title = element_text(size=15, face="bold", hjust=0.5)) + 
  scale_color_manual(values=c("grey77", "darkorange2", "slateblue3", "skyblue3"))
PCA_4 <- ggbiplot(Exudates.4, var.axes=FALSE, ellipse=TRUE, labels=H4$Sample, groups=H4$Tr, scale = 1) +
  theme_classic() + 
  ggtitle(Title4) + 
  theme(plot.title = element_text(size=15, face="bold", hjust=0.5)) + 
  scale_color_manual(values=c("grey77", "darkorange2", "slateblue3", "skyblue3"))
PCA_6 <-ggbiplot(Exudates.6, var.axes=FALSE, ellipse=TRUE, labels=H6$Sample, groups=H6$Tr, scale = 1) +
  theme_classic() + 
  ggtitle(Title6) + 
  theme(plot.title = element_text(size=15, face="bold", hjust=0.5)) + 
  scale_color_manual(values=c("grey77", "darkorange2", "slateblue3", "skyblue3"))

# PCA loadings
PCAloadings.2 <- data.frame(Variables=rownames(Exudates.2$rotation), Exudates.2$rotation)
PCAloadings.4 <- data.frame(Variables=rownames(Exudates.4$rotation), Exudates.4$rotation)
PCAloadings.6 <- data.frame(Variables=rownames(Exudates.6$rotation), Exudates.6$rotation)

# add PCA scores to the dataset
H2[, c('PC1', 'PC2')] = Exudates.2$x[, 1:2]
H4[, c('PC1', 'PC2')] = Exudates.4$x[, 1:2]
H6[, c('PC1', 'PC2')] = Exudates.6$x[, 1:2]

# save variable loadings in a separate dataset
rot.2 = as.data.frame(Exudates.2$rotation[, 1:2])
rot.2$var = rownames(Exudates.2$rotation)
rot.4 = as.data.frame(Exudates.4$rotation[, 1:2])
rot.4$var = rownames(Exudates.4$rotation)
rot.6 = as.data.frame(Exudates.6$rotation[, 1:2])
rot.6$var = rownames(Exudates.6$rotation)

# rescale the loadings to fit nicely within the scatterplot of our data
mult.2 = max(abs(H2[, c('PC1', 'PC2')])) / max(abs(rot.2[, 1:2])) / 2
rot.2[, 1:2] = rot.2[, 1:2] * mult.2
mult.4 = max(abs(H4[, c('PC1', 'PC2')])) / max(abs(rot.4[, 1:2])) / 2
rot.4[, 1:2] = rot.4[, 1:2] * mult.4
mult.6 = max(abs(H6[, c('PC1', 'PC2')])) / max(abs(rot.6[, 1:2])) / 2
rot.6[, 1:2] = rot.6[, 1:2] * mult.6

# ggplot the scatterplot and rotation taken from separate data.frames
# if there are many variables to plot, you can play with ggrepel 
library(ggrepel) 
Loading_2 <- ggplot(data = rot.2, aes(x = 0, y = 0, xend = PC1, yend = PC2, label = var)) +
  geom_segment(color = 'black', arrow = arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(aes(PC1 * 1, PC2 * 1)) +
  theme_bw() +
  theme(panel.grid = element_blank())
Loading_4 <- ggplot(data = rot.4, aes(x = 0, y = 0, xend = PC1, yend = PC2, label = var)) +
  geom_segment(color = 'black', arrow = arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(aes(PC1 * 1, PC2 * 1)) +
  theme_bw() +
  theme(panel.grid = element_blank())
Loading_6 <- ggplot(data = rot.6, aes(x = 0, y = 0, xend = PC1, yend = PC2, label = var)) +
  geom_segment(color = 'black', arrow = arrow(length = unit(0.02, "npc"))) +
  geom_label_repel(aes(PC1 * 1, PC2 * 1)) +
  theme_bw() +
  theme(panel.grid = element_blank())

#arrange
x <- ggarrange(PCA_2, PCA_4, PCA_6, Loading_2, Loading_4, Loading_6, 
          labels = c("A", "B", "C", "D", "E", "F"),
          align = "hv",
          ncol = 3, nrow = 2)
x

x1 <- ggarrange(Loading_2, Loading_4, Loading_6,
               align = "hv",
               ncol = 3, nrow = 1)
x1



# Save
ggsave(filename = paste("PCA", Title, Caption, ".pdf", sep="_"), 
       plot = last_plot(), 
       dpi = 600, units = "cm", 
       width = 50, height = 50, 
       scale = 0.4)

ggsave(filename = paste("PCA_Exudates_Loading", ".pdf", sep=""), 
       plot = x1, 
       dpi = 600, units = "cm", 
       width = 150, height = 50, 
       scale = 0.4)



# Treatment ####
C <- table[table$`Tr` == "C",]
xC <- C[1,2]
TitleC <- paste("Exudates", xC, sep="_")
P <- table[table$`Tr` == "P",]
xP <- P[1,2]
TitleP <- paste("Exudates", xP, sep="_")
Fe <- table[table$`Tr` == "F",]
xF <- Fe[1,2]
TitleF <- paste("Exudates", xF, sep="_")
M <- table[table$`Tr` == "M",]
xM <- M[1,2]
TitleM <- paste("Exudates", xM, sep="_")
Title <- "Exudates_Treatment"

# Principal components
Exudates.C <- prcomp(C[,c(4:nr.)], scale = TRUE) #(change number of columns -> compounds identified)
Exudates.P <- prcomp(P[,c(4:nr.)], scale = TRUE) #(change number of columns -> compounds identified)
Exudates.F <- prcomp(Fe[,c(4:nr.)], scale = TRUE) #(change number of columns -> compounds identified)
Exudates.M <- prcomp(M[,c(4:nr.)], scale = TRUE) #(change number of columns -> compounds identified)

# PCA
PAC_C <- ggbiplot(Exudates.C, var.axes=FALSE, ellipse=TRUE,  labels=C$Sample, groups=C$Time, scale = 1) +
  theme_classic() + 
  ggtitle(TitleC) + 
  theme(plot.title = element_text(size=20, face="bold", hjust=0.5)) +
  scale_color_manual(values=c("lightblue","steelblue","navy"))
PAC_P <- ggbiplot(Exudates.P, var.axes=FALSE, ellipse=TRUE,  labels=P$Sample, groups=P$Time, scale = 1) +
  theme_classic() + 
  ggtitle(TitleP) + 
  theme(plot.title = element_text(size=20, face="bold", hjust=0.5)) +
  scale_color_manual(values=c("lightblue","steelblue","navy"))
PAC_F <- ggbiplot(Exudates.F, var.axes=FALSE, ellipse=TRUE,  labels=Fe$Sample, groups=Fe$Time, scale = 1) +
  theme_classic() + 
  ggtitle(TitleF) + 
  theme(plot.title = element_text(size=20, face="bold", hjust=0.5)) +
  scale_color_manual(values=c("lightblue","steelblue","navy"))
PAC_M <- ggbiplot(Exudates.M, var.axes=FALSE, ellipse=TRUE,  labels=M$Sample, groups=M$Time, scale = 1) +
  theme_classic() + 
  ggtitle(TitleM) + 
  theme(plot.title = element_text(size=20, face="bold", hjust=0.5)) +
  scale_color_manual(values=c("lightblue","steelblue","navy"))

y <- ggarrange(PAC_C, PAC_P, PAC_F, PAC_M, 
          labels = c("A", "B", "C", "D"),
          align = "hv",
          ncol = 2, nrow = 2)
y

# Save
ggsave(filename = paste("PCA", Title, Caption, ".pdf", sep="_"), 
       plot = last_plot(), 
       dpi = 600, units = "cm", 
       width = 50, height = 50, 
       scale = 0.4)



# Principal components (complete dataset) ####
Exudates.pr <- prcomp(table[,c(4:nr.)], scale = TRUE) #(change number of columns -> compounds identified)
Title <- "Exudates"


# PCA plot
ggbiplot(Exudates.pr, var.axes=FALSE, ellipse=TRUE,  labels=table$Sample, groups=table$Tr, scale = 1) +
  theme_classic() + 
  ggtitle(Title) + 
  theme(plot.title = element_text(size=20, face="bold", hjust=0.5)) +
  scale_color_manual(values=c("grey77", "darkorange2", "slateblue3", "skyblue3"))

# Save
ggsave(filename = paste("PCA", Title, Caption, ".pdf", sep="_"), 
       plot = last_plot(), 
       dpi = 600, units = "cm", 
       width = 50, height = 50, 
       scale = 0.4)

# PCA loading plot ####
PCA_loadings <- Exudates.pr$rotation

plot(PCA_loadings, #PCA loading plot W dots
     pch=21,
     bg="black",
     cex=1,
     main="Root Exudates") +
  text(PCA_loadings, labels=rownames(PCA_loadings), 
       cex = 0.5,
       adj=0)

plot(PCA_loadings, #PCA loading plot W/O dots
     type="n",
     main= "Root Exudates (ISTD)") +
  text(PCA_loadings, labels=rownames(PCA_loadings), 
       cex = 0.5,
       adj=0)

# Save
ggsave(filename = paste("PCA_Exudates_Loading2.0", Caption, ".pdf", sep="_"), 
       plot = last_plot(), 
       dpi = 600, units = "cm", 
       width = 50, height = 50, 
       scale = 0.4)