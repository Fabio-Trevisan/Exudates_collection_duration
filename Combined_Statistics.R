library(readr)
library(ggpubr)
library(tidyverse)
library(plyr)
library(dplyr)
library(rstatix)
library(purrr)
library(agricolae)
library(reshape2)

#Exudates times statistics ####

#Read CSV ####
df <- read.csv("Exudates_ISTD_Area_.csv", sep=",",
                  header=T)
df <- df[!str_detect(names(df), "Sample")]
Caption <- "All"

df <- read.csv("Exudates_Significant.csv", sep=",",
               header=T)
Caption <- "Significant"

#Dataset Uniformation ####
names(df) <- substring(names(df), 3)
colnames(df)[colnames(df) == ""] <- "Tr"
colnames(df)[colnames(df) == "me"] <- "Time"
c <- c("Tr","Time")
table <- melt(df, na.rm = FALSE, value.name = "Value", id = c, variable = "RI_mz")

vector_RImz <- dput(colnames(df[,-(1:2)]))

table$RI_mz <- factor(table$RI_mz)
table$Time <- factor(table$Time)
table$Tr <- factor(table$Tr)

#create Subsets according to RI_mz ####
Subsets <- lapply(vector_RImz, function(i){ 
  i <- subset(table, RI_mz == i)
})

names(Subsets) <- vector_RImz



#Assumptions ####
## 1. Homogeneity of variances
##Tr*Time
Levene_test <- lapply(split(table, table$RI_mz), function(i){
  levene_test(Value ~ Tr * Time, data = i)
})

##2. Normality
##Shapiro-Wilk test for all single Trs
SW_test <- table %>%
  group_by(Time, Tr, RI_mz) %>%
  shapiro_test(Value)
View(SW_test)
write.table(SW_test, file = "13C_ShapiroWilk_test_results_2.0.csv", quote = FALSE, sep = ";")

##3. Indipendency
#Data are indepent by experimental design!
  
  
  
#2way ANOVA ####
##Multiple RI_mz
TwoWay_Anova <- lapply(split(table, table$RI_mz), function(i){
  anova(lm(Value ~ Tr * Time, data = i))
})
write.table(TwoWay_Anova, file = paste("TwoWay_Anova_", Caption, ".csv"), quote = FALSE, sep = ";")

sink(paste("TwoWay_Anova_2.0_", Caption, ".csv"))
TwoWay_Anova
sink(NULL)


#1way ANOVA (TIME)####
##Time.for tukey
OneWay_Anova_Ti <- lapply(vector_RImz, function(m){
  lapply(split(Subsets[[m]], Subsets[[m]][["Tr"]]), function(i){ 
    aov(Value ~ Time, data = i)
  })
})
names(OneWay_Anova_Ti) <- vector_RImz

##Time.for print
OneWay_Anova_Ti2 <- lapply(vector_RImz, function(m){
  lapply(split(Subsets[[m]], Subsets[[m]][["Tr"]]), function(i){ 
    anova(lm(Value ~ Time, data = i))
  })
})
names(OneWay_Anova_Ti2) <- vector_RImz

##OneWayAnova save
sink(paste("OneWayAnova_Ti_", Caption, ".csv" ))
OneWay_Anova_Ti2 
sink(NULL)



#Tukey as post hoc test (TIME)####
HSD_Ti <- lapply(vector_RImz, function(m){
  lapply(names(OneWay_Anova_Ti[[m]]), function(i){ 
    HSD.test(OneWay_Anova_Ti[[m]][[i]], "Time")
  })
})
names(HSD_Ti) <- vector_RImz
for(i in vector_RImz) {
  list <- names(OneWay_Anova_Ti[[i]]) 
  names(HSD_Ti[[i]]) <- list
}

##HSD_test save
HSD_Ti_groups <- lapply(vector_RImz, function(i){
  lapply(names(OneWay_Anova_Ti[[i]]), function(m){
    as.data.frame(HSD_Ti[[i]][[m]][["groups"]])
  })
})
names(HSD_Ti_groups) <- vector_RImz
for(i in vector_RImz) {
  list <- names(OneWay_Anova_Ti[[i]]) 
  names(HSD_Ti_groups[[i]]) <- list
}
sink(paste("HSD_Ti_", Caption, ".csv" ))
HSD_Ti_groups 
sink(NULL)



#1way ANOVA (TREATMENT)####
##Tr.for tukey
OneWay_Anova_Tr <- lapply(vector_RImz, function(m){
  lapply(split(Subsets[[m]], Subsets[[m]][["Time"]]), function(i){ 
    aov(Value ~ Tr, data = i)
  })
})
names(OneWay_Anova_Tr) <- vector_RImz


##Tr.for print
OneWay_Anova_Tr2 <- lapply(vector_RImz, function(m){
  lapply(split(Subsets[[m]], Subsets[[m]][["Time"]]), function(i){ 
    anova(lm(Value ~ Tr, data = i))
  })
})
names(OneWay_Anova_Tr2) <- vector_RImz


##OneWayAnova save
sink(paste("OneWayAnova_Tr_", Caption, ".csv" ))
OneWay_Anova_Tr2 
sink(NULL)




#Tukey as post hoc test (TREATMENT)####
HSD_Tr <- lapply(vector_RImz, function(m){
  lapply(names(OneWay_Anova_Tr[[m]]), function(i){ 
    HSD.test(OneWay_Anova_Tr[[m]][[i]], "Tr")
  })
})
names(HSD_Tr) <- vector_RImz
for(i in vector_RImz) {
  list <- names(OneWay_Anova_Tr[[i]]) 
  names(HSD_Tr[[i]]) <- list
}


##HSD_test save
HSD_Tr_groups <- lapply(vector_RImz, function(i){
  lapply(names(OneWay_Anova_Tr[[i]]), function(m){
    as.data.frame(HSD_Tr[[i]][[m]][["groups"]])
  })
})
names(HSD_Tr_groups) <- vector_RImz
for(i in vector_RImz) {
  list <- names(OneWay_Anova_Tr[[i]]) 
  names(HSD_Tr_groups[[i]]) <- list
}
sink(paste("HSD_Tr_", Caption, ".csv" ))
HSD_Tr_groups 
sink(NULL)
