library(devtools)
library(ggbiplot)
library(ggpubr)

table <- read.csv("Exudates_ISTD_Area_.csv", sep=",",
                  header=T)
Caption <- "Area"

table <- read.csv("Exudates_ISTD_Height_.csv", sep=",",
                    header=T)
Caption <- "Hight"



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
Exudates.2 <- prcomp(H2[,c(4:nr.)], scale = TRUE) #(change number of columns -> compounds identified)
Exudates.4 <- prcomp(H4[,c(4:nr.)], scale = TRUE) #(change number of columns -> compounds identified)
Exudates.6 <- prcomp(H6[,c(4:nr.)], scale = TRUE) #(change number of columns -> compounds identified)

# PCA
PCA_2 <- ggbiplot(Exudates.2, var.axes=FALSE, ellipse=TRUE, labels=H2$Sample, groups=H2$Tr, scale = 1) +
  theme_classic() + 
  ggtitle(Title2) + 
  theme(plot.title = element_text(size=20, face="bold", hjust=0.5)) + 
  scale_color_manual(values=c("grey77", "darkorange2", "slateblue3", "skyblue3"))
PCA_4 <- ggbiplot(Exudates.4, var.axes=FALSE, ellipse=TRUE, labels=H4$Sample, groups=H4$Tr, scale = 1) +
  theme_classic() + 
  ggtitle(Title4) + 
  theme(plot.title = element_text(size=20, face="bold", hjust=0.5)) + 
  scale_color_manual(values=c("grey77", "darkorange2", "slateblue3", "skyblue3"))
PCA_6 <-ggbiplot(Exudates.6, var.axes=FALSE, ellipse=TRUE, labels=H6$Sample, groups=H6$Tr, scale = 1) +
  theme_classic() + 
  ggtitle(Title6) + 
  theme(plot.title = element_text(size=20, face="bold", hjust=0.5)) + 
  scale_color_manual(values=c("grey77", "darkorange2", "slateblue3", "skyblue3"))

ggarrange(PCA_2, PCA_4, PCA_6, 
          labels = c("A", "B", "C"),
          align = "hv",
          ncol = 2, nrow = 2)

# Save
ggsave(filename = paste("PCA", Title, Caption, ".pdf", sep="_"), 
       plot = last_plot(), 
       dpi = 600, units = "cm", 
       width = 50, height = 50, 
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

ggarrange(PAC_C, PAC_P, PAC_F, PAC_M, 
          labels = c("A", "B", "C", "D"),
          align = "hv",
          ncol = 2, nrow = 2)

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
ggsave(filename = paste("PCA_Exudates_Loading", Caption, ".pdf", sep="_"), 
       plot = last_plot(), 
       dpi = 600, units = "cm", 
       width = 50, height = 50, 
       scale = 0.4)