library(ggpubr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(dplyr)

#Box plots

df <- read.csv("Exudates_ISTD_Area_.csv", sep=",",
               header=T)
df <- df[!str_detect(names(df), "Sample")]
Caption <- "All"

df <- read.csv("Exudates_Significant.csv", sep=",",
               header=T)
Caption <- "Significant"

#table<-table %>% select(1,2,3:to end)
table <- table[,-1]
table1<-melt(table, id=c("Tr","Time"))
melt(table, id=c("Tr","Time")) 
table_1_sum<-table1 %>% 
        group_by(Tr, Time, variable)


#comparisons
my_comparisons <-list(c("C","P"), 
                      c("C","F"),
                      c("C","M"),
                      c("P","F"),
                      c("P","M"),
                      c("F","M"))

my_comparisons <-list(c("C","P"), 
                      c("C","F"),
                      c("C","M"),
                      c("P","F"),
                      c("P","M"),
                      c("F","M"),
                      c("2","4"),
                      c("2","6"),
                      c("4","6"))



#Plots
table_1_sum$Time <- factor(table_1_sum$Time)

#Time
p2 <- ggplot(table_1_sum, aes(x= Time, value, fill=Tr)) + 
  theme_bw() + geom_boxplot(width=0.5) +
  facet_wrap(~variable, scales="free") + 
  ylab("Relative abundance") + xlab("Time") + 
  scale_fill_manual(values=c("grey77", "darkorange2", "slateblue3", "skyblue3"))
p2
t <- "Time"

#Treatment
p1 <- ggplot(table_1_sum, aes(x= Tr, value, fill=Time)) + 
  theme_bw() + geom_boxplot(width=0.5) +
  facet_wrap(~variable, scales="free") + 
  ylab("Relative abundance") + xlab("Treatments") + 
  scale_fill_manual(values=c("darkorange","red","brown"))
p1
t <- "Tr"
Stat <- "no-stat"

#not working
p1 + stat_compare_means(comparisons = my_comparisons, label = "p.ajust",  
                        method = "anova", 
                        p.adjust.method = "holm",
                        hide.ns = TRUE) +
  theme_bw() + ylab("Relative abundance")
Stat <- "stat"



# Save
ggsave(filename = paste("BoxPlot_Exudates", t, Stat, ".pdf", sep = "_"), plot = last_plot(), dpi = 600, units = "cm", width = 90, height = 90, scale = 1)

########
warnings()


