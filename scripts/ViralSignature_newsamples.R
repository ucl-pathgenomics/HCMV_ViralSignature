#Script to generate viral signature score from a new sample

library(tidyverse)
library(viridis)
library("ggsci")
library(here)
library(pROC)
library(ggbeeswarm)
library(ggpubr)
#Annotate data with cmvdrg (online or package)


#read data from paper 
#the data contains n of non-synonymous variants (minority) for each gene
DataPaper <-  read_csv(here("data","DataPaper_forModel_Clinical.csv"))
#select variables needed - 10 genes viral signature
listgenestotest <- c("sample_id","Patient_ID","outcome","UL97","UL74","UL7","UL11","UL20","UL75","UL54","UL121","UL8","UL37")
DataPaper_model <- DataPaper %>% select(listgenestotest)

#testing presence/absence of MV in the 10 genes - recode it variables
#(probably not needed)

#little function to transform n of variants into presence/absence 
update_columns <- function(data, columns) {
  for (col in columns) {
    data[[col]] <- ifelse(data[[col]] > 0, 1, 0)
  }
  return(data)
}

columns_to_update <- c("UL74", "UL75", "UL7", "UL121", "UL11", "UL20", "UL97", "UL54", "UL8", "UL37")
DataPaper_model <- update_columns(DataPaper_model, columns_to_update)
#recode outcome variable
DataPaper_model$status <- ifelse(DataPaper_model$outcome=="poor",1,0)
DataPaper_model$status <- as.factor(DataPaper_model$status)

#run model from the paper 
modsimple <- glm(status~UL97+UL74+UL7+UL8+UL37+UL11+UL20+UL75+UL121+UL54,DataPaper_model,family=binomial())
summary(modsimple)
#run model from the paper with only UL97 and UL54 
moddef <- glm(status~UL97+UL54,DataPaper_model,family=binomial())
summary(moddef)

#Calculate prob
prob=predict.glm(modsimple,type=c("response"))
DataPaper_model$prob=prob
#Calculate percentages
DataPaper_model$perc <- round(DataPaper_model$prob*100)

#plot % 
#ggplot(DataPaper_model,aes(status, perc,colour=status)) + geom_quasirandom() +
  #theme_pubr()

#New sample 
#Annotate new sample with cmvdrg online - see data_example for expected format
newsample <- read_csv(here("data_example","example_newsample.csv"))
newsample$Freq <- parse_number(newsample$freq) #parse freq as numbers

#Filter minority variants - the example was mapped to Merlin reference sequence so we are taking everything that is between 2% and 98%. 
#filtering also for non-synonymous and at least 5 reads supporting the variant
newsample_clean <- newsample %>% 
  filter(Freq > 2) %>% 
  filter(Freq < 98) %>% 
  filter(var_count > 4) %>%  
  filter(consequence=="nonsynonymous") %>% 
  select(change,var_count,Freq,consequence) %>% distinct()

newsample_clean$GENEID <- str_split_fixed(newsample_clean$change,"_",n=2)[,1]
genes <- c("UL97","UL74","UL7","UL11","UL20","UL75","UL54","UL121","UL8","UL37")
#select the 10 genes
newsample_clean_g <- newsample_clean[newsample_clean$GENEID %in% genes,]
newsample_clean_g$GENEID <- factor(newsample_clean_g$GENEID, levels= genes)
#count n of variants
newsample_count <- newsample_clean_g %>% group_by(GENEID,.drop=FALSE) %>%
  tally() %>%
  complete(GENEID)

#format wider and recode presence/absence
newsample_count_wider <- newsample_count %>% pivot_wider(names_from = GENEID, values_from = n)
newsample_count_wider <- update_columns(newsample_count_wider, columns_to_update)

#Calculate score for new sample using model created in the paper
pred.fit <-predict.glm(modsimple,newdata = newsample_count_wider ,type=c("response"),se.fit = TRUE)
newdata <- data_frame(sample_id="New_sample",Patient_ID="New_sample",outcome="New_sample",perc=pred.fit$fit*100,newsample_count_wider)
print(paste("The patient's score (HCMV viral signature) is ",newdata$perc))

#Make a plot
allData <- full_join(newdata,DataPaper_model)
allData$outcome <- factor(allData$outcome, levels = c("good","poor","New_sample"))
ggplot(allData,aes(outcome, perc,fill=outcome)) + 
  #geom_quasirandom() +
  geom_boxplot() +
  theme_pubr() +
  xlab("") +
  ylab("Viral signature (10 genes) score %") +
  theme(legend.position = "None") +
  scale_fill_jco() +
  scale_x_discrete(labels = c('Good Prognosis','Poor Prognosis','New Sample'))

