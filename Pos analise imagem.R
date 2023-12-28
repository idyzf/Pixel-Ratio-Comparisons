library(dplyr)
library(ggplot2)
library(openxlsx)
library(survminer)
library(survival)
rm(list=ls())

pixels<- read.csv("ratio_med.csv", sep=";")

pixels<- na.omit(pixels)
rownames(pixels)<- pixels$REDCAP_ID
counts<- read.xlsx("mRNA_Counts.xlsx")
rownames(counts)<- counts$REDCAP_ID

data<- merge(pixels, counts, by=0)


###################
clinical <- read.xlsx("clinical_v4-22-05-2023.xlsx")
rownames(clinical) <- clinical$ID
#data<- read.csv("ratio_med.csv", sep=";")
data$ID <- data$Row.names
names(clinical)

data1 <- left_join(clinical[, c(1, 21, 26, 75, 81)], data)

modelo <- lm(med_pixels ~ Número.de.ciclos, data=data1)
summary(modelo)

plot(data1$med_pixels, data1$Número.de.ciclos, main= "Correlation between Cyclos of TMZ and med pixels",
     xlab= "pixels", ylab = "Overall Survival")

abline(modelo)

p <- ggplot() + geom_violin(data1, mapping = aes(x=quantile, y=med_pixels, fill=quantile))
p + ggtitle("Pixels within the survival groups") +
  xlab("Survival Groups") + ylab("Pixels")


ratios<- read.csv("Ratio_teste.csv", sep = ";")
ratios$patients
test_patients <- c("E09-11-A5", "E13-15568-B3", "E10-6562-B1", "HCB-151")

train <- subset(ratios, ratios$patients %in% test_patients)
train$Groups <- train$patients
train<- mutate(train, train$Groups)

train <- train %>% 
  mutate(Groups = recode(Groups, 
                         "E09-11-A5" = "Low",
                         "E10-6562-B1" = "Mid",
                         "E13-15568-B3" = "No",
                         "HCB-151" = "High"))

p <- ggplot() + geom_violin(train, mapping = aes(x=Groups, y=Ratio_N, fill=Groups))
p + ggtitle("Pixels within the survival groups") +
  xlab("Survival Groups") + ylab("Pixels")

data1<- na.omit(data1)
med<- median(data1$CDK7)
data1$median <- data1$CDK7
data1$median[which(data1$median > med)] <- "high"
data1$median[which(data1$median != "high")] <- "low"

p <- ggplot() + geom_violin(data1, mapping = aes(x=quantile, y=CDK7, fill=quantile))
p + ggtitle("Pixels within the CDK7 expression Groups") +
  xlab("CDK7 expression Groups") + ylab("Pixels")
#####################################################
################Kapplan Meier########################
med_pix <- median(data1$med_pixels)
data1<- na.omit(data1)
data1$median_pix <- data1$med_pixels
data1$median_pix[which(data1$median_pix > med_pix)] <- "high"
data1$median_pix[which(data1$median_pix != "high")] <- "low"

data1$status <- factor(data1$"Status", 
                           levels = c("Vivo com doenca", "obito por câncer", "obito por outras causas com doenca ativa"), 
                           labels = c("1", "2", "2"))



data1$status <- as.numeric(as.character(data1$status))
data1$median_pix <- factor(data1$median_pix)
surv_object <- Surv(time = data1$Sobrevida.Global.sem.arredondamento, 
                    event = data1$status)
surv_object 
data_esp <- data1$`clinical$OS`
fit1 <- survfit(surv_object ~ median_pix, data = data1)
summary(fit1)

ggsurvplot (fit1 , pval = TRUE, title= "", legend.labs=c("High","Low"), 
            legend.title = "Median of pixels")
