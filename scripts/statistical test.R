library(tidyverse)
library(car)        
library(rstatix)    
library(ggpubr)
library(ggplot2)
dir.create("statistical test")
dir.create("statistical test/adult")
#Downloading the data 
adult_data<-read.csv("enrichement analysis/adult/STRING_analysis.csv",header = T,sep = ",")
# exploring the data statistically 
summary(adult_data)
sum(is.na(adult_data$disorder_percentage))
#visualization
ggplot(adult_data, aes(x = disorder_percentage)) +
  geom_histogram(
    bins = 30,
    color = "black"
  ) +
  theme_classic() +
  labs(
    title = "Distribution of IDR percentage in adult DEGs",
    x = "IDR percentage",
    y = "Number of genes"
  )

#QQplot
ggqqplot(adult_data$disorder_percentage,title = "Q-Q plot of IDR percentage")
#Normality test
shapiro.test(adult_data$disorder_percentage) 
#correlation between logFC and IDR_percentage--------------------------------------------------------
cor.test(
  adult_data$logFC,
  adult_data$disorder_percentage,
  method = "spearman"
)
#scatter plot
ggplot(adult_data, aes(x = logFC, y = disorder_percentage)) +
  geom_point(size=5,alpha = 0.7) +
  labs(
    title = "Relationship between logFC and Disorder Percentage",
    x = "log2 Fold Change",
    y = "Disorder Percentage (%)"
  ) +
  theme_classic()
#Normality test for the degree--------------------------------------------------------------------------
shapiro.test(adult_data$Degree) 
#correlation between logFC and Degree
cor.test(
  adult_data$logFC,
  adult_data$Degree,
  method = "spearman",
  alternative = "greater"
)
ggplot(adult_data, aes(x = logFC, y = Degree)) +
  geom_point(size = 4) +
  theme_classic() +
  labs(
    x = "log2 Fold Change",
    y = "STRING Degree",
    title = "logFC vs STRING Degree"
  )
#Normality test for the length------------------------------------------------------------------------------
shapiro.test(adult_data$length) 
#correlation between logFC and Degree
cor.test(
  adult_data$length,
  adult_data$disorder_percentage,
  method = "spearman")
ggplot(adult_data, aes(x = length, y = disorder_percentage)) +
  geom_point(size = 4) +
  theme_classic() +
  labs(
    x = "Protein length",
    y = "IDR percentage",
    title = "Protein length vs IDR percentage"
  )
#regression model between logFC and degree, length and IDR percentage
model<-lm(logFC~disorder_percentage+length+Degree,data = adult_data)
summary(model)
coeff_table <- as.data.frame(summary(model)$coefficients)
coeff_table
write.csv(coeff_table,"statistical test/adult/coeff_table.csv",row.names = F)
#Pediatric part--------------------------------------------------------------------------------------
#Downloading the data 
pediatric_data<-read.csv("enrichement analysis/pediatric/STRING_analysis.csv",header = T,sep = ",")
# exploring the data statistically 
summary(pediatric_data)
sum(is.na(pediatric_data$disorder_percentage))
#omit genes with NA IDR_percentage
pediatric_IDR_clean <- pediatric_data %>%
  filter(!is.na(disorder_percentage))
dim(pediatric_IDR_clean)
#visualization
ggplot(pediatric_IDR_clean, aes(x = disorder_percentage)) +
  geom_histogram(
    bins = 30,
    color = "black"
  ) +
  theme_classic() +
  labs(
    title = "Distribution of IDR percentage in pediatric DEGs",
    x = "IDR percentage",
    y = "Number of genes"
  )

#QQplot
ggqqplot(pediatric_IDR_clean$disorder_percentage,title = "Q-Q plot of IDR percentage")
#Normality test-----------------------------------------------------------------------------------------
shapiro.test(pediatric_IDR_clean$disorder_percentage) # IDR_percentage non normale avec pvalue de 2.329e-06
#correlation between avg_log2FC and IDR_percentage
cor.test(
  pediatric_IDR_clean$logFC,
  pediatric_IDR_clean$disorder_percentage,
  method = "spearman"
)
ggplot(pediatric_IDR_clean, 
       aes(x = logFC, y = disorder_percentage)) +
  geom_point(size=5, alpha = 0.6) +
  theme_classic() +
  labs(
    title = "Relationship between logFC and Disorder Percentage",
    x = "Average log2 fold change",
    y = "IDR percentage (%)"
  )
#Normality test for the degree--------------------------------------------------------------------------
shapiro.test(pediatric_IDR_clean$Degree) 
#correlation between logFC and Degree
cor.test(
  pediatric_IDR_clean$logFC,
  pediatric_IDR_clean$Degree,
  method = "spearman")
#visualization
ggplot(pediatric_IDR_clean, aes(x = logFC, y = Degree)) +
  geom_point(size = 4) +
  theme_classic() +
  labs(
    x = "log2 Fold Change",
    y = "STRING Degree",
    title = "logFC vs STRING Degree"
  )
#Normality test for the length---------------------------------------------------------------------------
shapiro.test(pediatric_IDR_clean$length) 
#correlation between length and IDR
cor.test(
  pediatric_IDR_clean$length,
  pediatric_IDR_clean$disorder_percentage,
  method = "spearman")
ggplot(pediatric_IDR_clean, aes(x = length, y = disorder_percentage)) +
  geom_point(size = 4) +
  theme_classic() +
  labs(
    x = "Protein length",
    y = "IDR percentage",
    title = "Protein length vs IDR percentage"
  )
#regression model between logFC and degree, length and IDR percentage
model1<-lm(logFC~disorder_percentage+length+Degree,data = pediatric_IDR_clean)
summary(model1)
coeff_table1 <- as.data.frame(summary(model1)$coefficients)
coeff_table1
write.csv(coeff_table1,"statistical test/pediatric/coeff_table.csv",row.names = F)
#statistical analysis between adult and pediatric----------------------------------------------------------
#adding corresponding pattern to each group
adult_data$Group <- "Adult"
pediatric_IDR_clean$Group <- "Pediatric"
#merging the two datasets
IDR_comparison <- rbind(
  adult_data,
  pediatric_IDR_clean
)
#-------------------------------------------------------------------------------------------------------------
#doing non-parametric test : wilcoxon: to see if the IDR-percentage differ between adult and pediatric
wilcox.test(disorder_percentage ~ Group, data = IDR_comparison) 

#visualization
ggplot(IDR_comparison,aes(x = Group,y = disorder_percentage,fill = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  theme_classic()

#calculating effect size
IDR_comparison %>%
  wilcox_effsize(disorder_percentage ~ Group) # effsize of 0.0510 <0.1 meaning small effect and slight difference between adult and pediatric group
#--------------------------------------------------------------------------------------------------------------
#doing non-parametric test : wilcoxon: to see if the protein length differ between adult and pediatric
wilcox.test(length ~ Group, data = IDR_comparison) 

#visualization
ggplot(IDR_comparison,aes(x = Group,y = length,fill = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  theme_classic()

#calculating effect size
IDR_comparison %>%
  wilcox_effsize( length~ Group) # effsize of 0.157 <0.1 meaning small effect and slight difference between adult and pediatric group
#-------------------------------------------------------------------------------------------------------------
#doing non-parametric test : wilcoxon: to see if the protein length differ between adult and pediatric
wilcox.test(Degree ~ Group, data = IDR_comparison) 

#visualization
ggplot(IDR_comparison,aes(x = Group,y = Degree,fill = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  theme_classic()

#calculating effect size
IDR_comparison %>%
  wilcox_effsize( Degree~ Group) # effsize of 0.238 <0.1 meaning small effect and slight difference between adult and pediatric group
save.image("statistical test/Stats.RData")
