library(readr)
library(ggpubr)
library(tidyverse)
library(plyr)
library(dplyr)
library(rstatix)
library(purrr)
library(agricolae)


#Fresh weight statistics (extended)####

#for Time ####
#Read CSV ####
table <- read.csv("DATA_weight.csv", sep=";",
                  header=T)

table$Tr <- factor(table$Tr)
table$Time <- factor(table$Time)

#Summary Table ####
Summary_table <- ddply(table, c("Tr", "Time"), summarise,
                       N    = sum(!is.na(FW)),
                       mean = mean(FW, na.rm=TRUE),
                       sd   = sd(FW, na.rm=TRUE),
                       se   = sd / sqrt(N))
Summary_table
write.table(Summary_table, file = "Weight_Summary_table.csv", quote = FALSE, sep = ";")


#Assumptions ####
## 1. Homogeneity of variances
##Tr*Time
Levene_test <- levene_test(FW ~ Tr * Time, data = table)
Levene_test_Time <- levene_test(FW ~ Time, data = table)
Levene_test_Tr <- levene_test(FW ~ Tr, data = table)

sink("Weight_Levene_test_Homogeneity.csv")
"Interaction Tr*Time"
Levene_test
"Time"
Levene_test_Time
"Tr"
Levene_test_Tr
sink(NULL)

##2. Normality
##Shapiro-Wilk test for all single Trs
SW_test <- table %>%
  group_by(Time, Tr) %>%
  shapiro_test(FW)
View(SW_test)
write.table(SW_test, file = "ShapiroWilk_test_Normality.csv", quote = FALSE, sep = ";")

##3. Indipendency
#Data are indepent by experimental design!



#2way ANOVA ####
##Multiple RI_mz
TwoWay_Anova <- anova(lm(FW ~ Tr * Time, data = table))

write.table(TwoWay_Anova, file = "Weight_TwoWay_Anova.csv", quote = FALSE, sep = ";")

sink("Weight_TwoWay_Anova.csv")
TwoWay_Anova
sink(NULL)


#1way ANOVA (TIME)####
vector_Tr <- c("C", "P", "Fe", "M")
vector_Time <- c("2", "4", "6")
  
##Time.for tukey
OneWay_Anova_Ti <- lapply(split(table, table$Tr), function(i){ 
    aov(FW ~ Time, data = i)
})

##Time.for print
OneWay_Anova_Ti2 <- lapply(split(table, table$Tr), function(i){ 
  anova(lm(FW ~ Time, data = i))
})

##OneWayAnova save
sink("Weight_OneWayAnova_Ti.csv")
OneWay_Anova_Ti2 
sink(NULL)



#Tukey as post hoc test (TIME)####
HSD_Ti <- lapply(vector_Tr, function(m){
   HSD.test(OneWay_Anova_Ti[[m]], "Time")
})
names(HSD_Ti) <- vector_Tr


##HSD_test save
HSD_Ti_groups <- lapply(vector_Tr, function(i){
  as.data.frame(HSD_Ti[[i]][["groups"]])
})
names(HSD_Ti_groups) <- vector_Tr

sink("Weight_HSD_Ti.csv")
HSD_Ti_groups 
sink(NULL)



#1way ANOVA (TREATMENT)####
##Tr.for tukey
OneWay_Anova_Tr <- lapply(split(table, table$Time), function(i){ 
    aov(FW ~ Tr, data = i)
})


##Tr.for print
OneWay_Anova_Tr2 <- lapply(split(table, table$Time), function(i){ 
    anova(lm(FW ~ Tr, data = i))
})


##OneWayAnova save
sink("Weight_OneWayAnova_Tr.csv")
OneWay_Anova_Tr2 
sink(NULL)




#Tukey as post hoc test (TREATMENT)####
HSD_Tr <- lapply(vector_Time, function(m){
   HSD.test(OneWay_Anova_Tr[[m]], "Tr")
})
names(HSD_Tr) <- vector_Time


##HSD_test save
HSD_Tr_groups <- lapply(vector_Time, function(i){
   as.data.frame(HSD_Tr[[i]][["groups"]])
})
names(HSD_Tr_groups) <- vector_Time


sink("Weight_HSD_Tr.csv")
HSD_Tr_groups 
sink(NULL)