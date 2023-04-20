library(dplyr)
library(scales)
library(tidyr)
library(stringr)


#### ROOT EXUDATES ####
table <- read.table("Exudates_Height.txt", sep="\t",
                    header=F)

table <- read.table("Exudates_Area.txt", sep="\t",
                    header=F)

# Columns and lines cleaning ####
########## add cleaning of pool and pyr samples

table1 <- table[-c(1:4),-c(1,6:15,17:25,27:28)] 

names(table1) <-  unlist(table1[1,]) #rename dataset
table2 <- table1[-1,]

x <- 7 #number of categories created in MSDial
y <- c((ncol(table2)-x*2+1):ncol(table2)) #vector with columns to erase
table3 <- table2[,-y]

table4 <- table3[table3$Comment!="",] #removal unchecked peaks
table4 <- table4[table4$Comment!="no",] #removal checked peaks labelled with "no"

table4[,c(1:3,6:88)] <- sapply(table4[,c(1:3,6:88)], as.numeric) #columns from character to numeric



# S/N ratio filtering ####
z <- 3 #filtering ratio
SN_filtered <- table4[table4$`S/N average`>z,] #S/N ratio filtering >z
SN_filtered <- SN_filtered[,-6]



# Sample/Blank ratio filtering ####
########## Implement loop / apply to avoid repetition of commands
Blank_filtering <- SN_filtered

#caluclate means --> enter number of means to calculate according to the treatments
ExtrBlank <- Blank_filtering %>% select(matches("ExtrBlank"))  
C2 <- Blank_filtering %>% select(matches("C2")) 
C4 <- Blank_filtering %>% select(matches("C4")) 
C6 <- Blank_filtering %>% select(matches("C6"))
F2 <- Blank_filtering %>% select(matches("F2")) 
F4 <- Blank_filtering %>% select(matches("F4")) 
F6 <- Blank_filtering %>% select(matches("F6")) 
P2 <- Blank_filtering %>% select(matches("P2")) 
P4 <- Blank_filtering %>% select(matches("P4")) 
P6 <- Blank_filtering %>% select(matches("P6")) 
M2 <- Blank_filtering %>% select(matches("M2")) 
M4 <- Blank_filtering %>% select(matches("M4")) 
M6 <- Blank_filtering %>% select(matches("M6")) 

Blank_filtering$ExtrBlank_Mean <- rowMeans(ExtrBlank[], na.rm=TRUE)
Blank_filtering$C2_Mean <- rowMeans(C2[], na.rm=TRUE) 
Blank_filtering$C4_Mean <- rowMeans(C4[], na.rm=TRUE) 
Blank_filtering$C6_Mean <- rowMeans(C6[], na.rm=TRUE) 
Blank_filtering$F2_Mean <- rowMeans(F2[], na.rm=TRUE) 
Blank_filtering$F4_Mean <- rowMeans(F4[], na.rm=TRUE) 
Blank_filtering$F6_Mean <- rowMeans(F6[], na.rm=TRUE) 
Blank_filtering$P2_Mean <- rowMeans(P2[], na.rm=TRUE) 
Blank_filtering$P4_Mean <- rowMeans(P4[], na.rm=TRUE) 
Blank_filtering$P6_Mean <- rowMeans(P6[], na.rm=TRUE) 
Blank_filtering$M2_Mean <- rowMeans(M2[], na.rm=TRUE) 
Blank_filtering$M4_Mean <- rowMeans(M4[], na.rm=TRUE) 
Blank_filtering$M6_Mean <- rowMeans(M6[], na.rm=TRUE) 

#caluclate Sample/bBank ratio --> enter number of means to calculate according to the treatments
Blank_filtering$C2_Ratio <- Blank_filtering$C2_Mean/Blank_filtering$ExtrBlank_Mean 
Blank_filtering$C4_Ratio <- Blank_filtering$C4_Mean/Blank_filtering$ExtrBlank_Mean  
Blank_filtering$C6_Ratio <- Blank_filtering$C6_Mean/Blank_filtering$ExtrBlank_Mean  
Blank_filtering$F2_Ratio <- Blank_filtering$F2_Mean/Blank_filtering$ExtrBlank_Mean  
Blank_filtering$F4_Ratio <- Blank_filtering$F4_Mean/Blank_filtering$ExtrBlank_Mean 
Blank_filtering$F6_Ratio <- Blank_filtering$F6_Mean/Blank_filtering$ExtrBlank_Mean 
Blank_filtering$P2_Ratio <- Blank_filtering$P2_Mean/Blank_filtering$ExtrBlank_Mean  
Blank_filtering$P4_Ratio <- Blank_filtering$P4_Mean/Blank_filtering$ExtrBlank_Mean  
Blank_filtering$P6_Ratio <- Blank_filtering$P6_Mean/Blank_filtering$ExtrBlank_Mean  
Blank_filtering$M2_Ratio <- Blank_filtering$M2_Mean/Blank_filtering$ExtrBlank_Mean  
Blank_filtering$M4_Ratio <- Blank_filtering$M4_Mean/Blank_filtering$ExtrBlank_Mean 
Blank_filtering$M6_Ratio <- Blank_filtering$M6_Mean/Blank_filtering$ExtrBlank_Mean 

#Sample/Blank ratio filtration
Blank_filtering$Max_Ratio <- pmax(Blank_filtering$C2_Ratio, 
                                  Blank_filtering$C4_Ratio, 
                                  Blank_filtering$C6_Ratio,
                                  Blank_filtering$F2_Ratio, 
                                  Blank_filtering$F4_Ratio, 
                                  Blank_filtering$F6_Ratio, 
                                  Blank_filtering$P2_Ratio, 
                                  Blank_filtering$P4_Ratio, 
                                  Blank_filtering$P6_Ratio, 
                                  Blank_filtering$M2_Ratio, 
                                  Blank_filtering$M4_Ratio, 
                                  Blank_filtering$M6_Ratio)

w <- 3 #filtering ratio
Blank_filtered <- Blank_filtering[Blank_filtering$`Max_Ratio`>w,] #Sample/Blank ratio filtration > w

#clening ratios and means from dataset
Blank_filtered <- Blank_filtered[,1:65] #depends on number of samples check



# ISTD normalization #### 
Metabolites_info <- Blank_filtered[,1:5] #to add at the end 
df <- Blank_filtered[,6:65] #select only numerical variables
row.names(df) <-  unlist(Blank_filtered[,4]) #rename rows

ISTD <- df["*13C Citric acid (4TMS)",] #select ISTD
ISTD <- ISTD[rep(1, 46), ] #create df of equal size
ISTD_Norm <- df/ISTD #dived the 2 df
#poco elegante ma funziona
ISTD_Norm <- ISTD_Norm[!(row.names(ISTD_Norm) %in% "*13C Citric acid (4TMS)"), ] #remove ISTD line



# FW or root weight normalization #### 
FW <- read.csv("DATA_GC_exudates_Root_weight.csv", sep=";", header=F) #get FW or root weight
row.names(FW) <-  unlist(FW[,1]) #rename rows
names(FW) <-  unlist(FW[1,]) #rename columns
FW <- FW[-1,-1] #remove headings columns/rows
FW[] <- sapply(FW[], as.numeric) #set numeric variables
FW <- FW[rep(1,45),] #create df of equal size

# ISTD
FW_norm <- ISTD_Norm/FW #dived the 2 df

# no-ISTD
df <- df[!(row.names(df) %in% "*13C Citric acid (4TMS)"), ] #remove ISTD line
FW_norm <- df/FW #dive the 2 df



# Removal of outlayers####
#root exudates (16_M2 & 48_P6)
FW_norm <- FW_norm[,!(names(FW_norm) %in% c("16_M2","48_P6"))]



# Relative abundance ####
Relative_abundance <- as.data.frame(t(apply(FW_norm, 1, function(x) x/max(x))))



# Dataset cleaning #### 
# DATA_GC_Exudates_ISTD-correction
Clean_df <- data.frame(t(Relative_abundance)) #invert row/columns
Clean_df <- data.frame(Names = row.names(Clean_df), Clean_df) #extract rownames into column
rownames(Clean_df) <- NULL

Clean_df <- Clean_df %>%
  separate(Names, c("Sample_Number", "Tr_Time")) #separate Sample Number from Treatment/Time
Clean_df <- extract(Clean_df, Tr_Time, into = c("Tr", "Time"), "(.{1})(.{1})", remove=T) #separate Treatment from Time

# DATA_GC_Exudates_noISTD-correction
Clean_df2 <- data.frame(t(Relative_abundance)) #invert row/columns
Clean_df2 <- data.frame(Names = row.names(Clean_df2), Clean_df2) #extract rownames into column
rownames(Clean_df2) <- NULL

Clean_df2 <- Clean_df2 %>%
  separate(Names, c("Sample_Number", "Tr_Time")) #separate Sample Number from Treatment/Time
Clean_df2 <- extract(Clean_df2, Tr_Time, into = c("Tr", "Time"), "(.{1})(.{1})", remove=T) #separate Treatment from Time

# DATA_GC_Exudates_noRelativeAbundance_ISTD-correction
Clean_df3 <- data.frame(t(FW_norm)) #invert row/columns
Clean_df3 <- data.frame(Names = row.names(Clean_df3), Clean_df3) #extract rownames into column
rownames(Clean_df3) <- NULL

Clean_df3 <- Clean_df3 %>%
  separate(Names, c("Sample_Number", "Tr_Time")) #separate Sample Number from Treatment/Time
Clean_df3 <- extract(Clean_df3, Tr_Time, into = c("Tr", "Time"), "(.{1})(.{1})", remove=T) #separate Treatment from Time





# Save table ####
write.csv(Clean_df, "DATA_GC_Exudates_ISTD-correction_Area.csv", row.names=FALSE)

write.csv(Clean_df, "DATA_GC_Exudates_noISTD-correction_Hight.csv", row.names=FALSE)

write.csv(Clean_df3, "DATA_GC_Exudates_noRelativeAbundance_ISTD-correction.csv", row.names=FALSE)


