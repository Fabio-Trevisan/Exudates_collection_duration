library(readr)
library(ggpubr)
library(tidyverse)
library(plyr)
library(dplyr)
library(rstatix)
library(purrr)
library(agricolae)
library(reshape2)


#Statistics C vs P vs F vs M ####

#Read CSV ####
#Exudates (Tr, Time, Sample_Number)
df <- read.csv("Exudates_ISTD_Clean.csv", sep=",",
                  header=T)
Tissue <- "Exudates"
df <- df[!str_detect(names(df), "Time")]
df <- df[!str_detect(names(df), "Sample")]



#Dataset Uniformation ####
names(df) <- substring(names(df), 3)
colnames(df)[colnames(df) == ""] <- "Tr"
table <- melt(df, na.rm = FALSE, value.name = "Relative_abbundance", id = 'Tr', variable = "RI_mz")

vector_RImz <- dput(colnames(df[,-1]))

table$RI_mz <- as.character(table$RI_mz)


#Assumptions ####
## 1. Homogeneity of variances
##Treatment
Levene_test <- lapply(split(table, table$RI_mz), function(i){
  levene_test(Relative_abbundance ~ Tr, data = i)
})
sink("Exudates_Levene_test_(Homogeneity).csv")
Levene_test 
sink(NULL)

##2. Normality
##Shapiro-Wilk test for all single treatments
SW_test <- table %>%
  group_by(Tr, RI_mz) %>%
  shapiro_test(Relative_abbundance)
View(SW_test)
write.table(SW_test, file = "Exudates_ShapiroWilk_test_(Normality).csv", quote = FALSE, sep = ";")

##3. Indipendency
#Data are indepent by experimental design!
  
  

#1way ANOVA ####
## for PostHoc Test
OneWay_Anova <- lapply(split(table, table$RI_mz), function(i){ 
    aov(Relative_abbundance ~ Tr, data = i)
})

## for Print
OneWay_Anova_Print <- lapply(split(table, table$RI_mz), function(i){ 
  anova(lm(Relative_abbundance ~ Tr, data = i))
})

##OneWayAnova save
sink(paste(Tissue, "OneWayAnova.csv", sep="_"))
OneWay_Anova_Print 
sink(NULL)



# P-Value print in original table ####
Anova_PValue <- lapply(vector_RImz, function(i){
  as.data.frame(OneWay_Anova_Print[[i]][["Pr(>F)"]])
})
names(Anova_PValue) <- vector_RImz #extract P_Value as list

for(i in vector_RImz) { #rename sub-title list to P-Value
  names(Anova_PValue[[i]]) <- "P_Value"
} 

P_Value <- data.frame(Anova_PValue) #transform list in dataframe
names(P_Value) <- vector_RImz #rename Columns
P_Value <- P_Value[-2,] #remove Row

df2 <- as.data.frame(t(df)) #reverse dataframe
names(df2) <-  unlist(df2[1,]) #rename Columns
df2 <- df2[-1,] #remove Row

df2$P_Value <- t(P_Value) #add column with P_Value to original table

z <- 0.05 #filtering ratio
Filtered <- df2[df2$`P_Value`<z,] #P_Value filtering >z
Filtered <- Filtered[,-19]

Clean_df <- data.frame(t(Filtered)) #invert row/columns
Clean_df <- data.frame(Names = row.names(Clean_df), Clean_df) #extract rownames into column
rownames(Clean_df) <- NULL #clean row names
Clean_df <- extract(Clean_df, Names, into = c("Tr"), "(.{1})", remove=T) #separate Treatment 
write.csv(Clean_df, file = paste(Tissue, "Significant.csv", sep="_"), row.names=FALSE)




#Tukey as post hoc test ####
HSD <- lapply(vector_RImz, function(i){ 
    HSD.test(OneWay_Anova[[i]], "Tr")
})
names(HSD) <- vector_RImz


##HSD_test save
HSD_groups <- lapply(vector_RImz, function(i){
    as.data.frame(HSD[[i]][["groups"]])
})
names(HSD_groups) <- vector_RImz

sink(paste(Tissue, "HSD.csv", sep="_"))
HSD_groups 
sink(NULL)
