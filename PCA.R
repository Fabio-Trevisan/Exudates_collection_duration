library(devtools)
library(ggbiplot)



table <- read.csv("DATA_GC_Exudates_ISTD-correction_Hight.csv", sep=",",
                    header=T)
Caption <- "(ISTD_Hight)"
Title <- ("Exudates")
  
table <- read.csv("DATA_GC_Tissues_ISTD-correction_Area.csv", sep=",",
                  header=T)
Caption <- "(ISTD_Area)"
Title <- ("GC_Root")
Title <- ("GC_Shoot")

nr. = ncol(table)



#Groups definition or combination ####
H2 <- table[table$`Time` == 2,] #S/N ratio filtering = 2
H4 <- table[table$`Time` == 4,] #S/N ratio filtering = 4
H6 <- table[table$`Time` == 6,] #S/N ratio filtering = 6

R <- table[table$`Tissue` == "R",] #Selecting Root
S <- table[table$`Tissue` == "S",] #Selecting Shoot



# Principal components ####

# GC
Root.pr <- prcomp(R[,c(4:nr.)], scale = TRUE) #(change number of columns -> compounds identified)
Shoot.pr <- prcomp(S[,c(4:nr.)], scale = TRUE) #(change number of columns -> compounds identified)

# Exudates
Exudates.pr <- prcomp(table[,c(4:nr.)], scale = TRUE) #(change number of columns -> compounds identified)
Exudates.pr <- prcomp(H2[,c(4:nr.)], scale = TRUE) #(change number of columns -> compounds identified)
Exudates.pr <- prcomp(H4[,c(4:nr.)], scale = TRUE) #(change number of columns -> compounds identified)
Exudates.pr <- prcomp(H6[,c(4:nr.)], scale = TRUE) #(change number of columns -> compounds identified)

Summary <- summary(Exudates.pr) #summary results PCA (variance explained by each Principal Component)



# PCA plot ####
# GC Root
ggbiplot(Root.pr, var.axes=FALSE, ellipse=TRUE,  labels=R$Sample, groups=R$Tr, scale = 1) +
  theme_classic() + 
  ggtitle(Title) + 
  theme(plot.title = element_text(size=20, face="bold", hjust=0.5)) + 
  labs(caption = Caption) 

# GC Shoot
ggbiplot(Shoot.pr, var.axes=FALSE, ellipse=TRUE,  labels=S$Sample, groups=S$Tr, scale = 1) +
  theme_classic() + 
  ggtitle(Title) + 
  theme(plot.title = element_text(size=20, face="bold", hjust=0.5)) + 
  labs(caption = Caption) 

#Save
ggsave(filename = paste(Title, Caption, ".pdf", sep="_"), 
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







