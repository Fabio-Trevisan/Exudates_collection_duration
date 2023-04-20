library(devtools)
library(ggbiplot)

table <- read.csv("DATA_GC_Exudates_ISTD-correction_Area.csv", sep=",",
                  header=T)
Caption <- "(ISTD_Area)"

table <- read.csv("DATA_GC_Exudates_ISTD-correction_Hight.csv", sep=",",
                    header=T)
Caption <- "(ISTD_Hight)"

nr. = ncol(table)

table$Time <- as.factor(table$Time)



#Groups definition or combination ####
#Time
a <- table[table$`Time` == "2",]
a <- table[table$`Time` == "4",]
a <- table[table$`Time` == "6",]
x <- a[1,3]

#Tr
a <- table[table$`Tr` == "C",]
a <- table[table$`Tr` == "P",]
a <- table[table$`Tr` == "F",]
a <- table[table$`Tr` == "M",]
x <- a[1,2]

Title <- paste("Exudates", x, sep="_")



# Principal components ####
Exudates.pr <- prcomp(table[,c(4:nr.)], scale = TRUE) #(change number of columns -> compounds identified)

Exudates.a <- prcomp(a[,c(4:nr.)], scale = TRUE) #(change number of columns -> compounds identified)



# PCA plot ####
ggbiplot(Exudates.pr, var.axes=FALSE, ellipse=TRUE,  labels=table$Sample, groups=table$Tr, scale = 1) +
  theme_classic() + 
  ggtitle(Title) + 
  theme(plot.title = element_text(size=20, face="bold", hjust=0.5)) + 
  labs(caption = Caption) 

# Time
ggbiplot(Exudates.a, var.axes=FALSE, ellipse=TRUE,  labels=a$Sample, groups=a$Tr, scale = 1) +
  theme_classic() + 
  ggtitle(Title) + 
  theme(plot.title = element_text(size=20, face="bold", hjust=0.5)) + 
  labs(caption = Caption) 

# Tr
ggbiplot(Exudates.a, var.axes=FALSE, ellipse=TRUE,  labels=a$Sample, groups=a$Time, scale = 1) +
  theme_classic() + 
  ggtitle(Title) + 
  theme(plot.title = element_text(size=20, face="bold", hjust=0.5)) + 
  labs(caption = Caption) 



# Save ####
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







