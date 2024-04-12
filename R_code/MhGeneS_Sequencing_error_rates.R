########################################################################################
########################################################################################
##### This R script processes Error rates for each position and sample     #############
#####      in the Amplicon Sequence produced by Seq2Sat                   ##############
########################################################################################
########################################################################################

setwd("YOUR_WORKING_DIRECTORY/")


library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library("ggrepel") 


##############################################################################
#### 1. Get mean Error rates per position within the amplicon ################
##############################################################################

# List all the Seq2Sat output files containing the sequencing error rates per position and the average error rate (per sample)

Testdata_filenames <- list.files(path = "~/PATH/testdata_ErrorRate_files/", pattern="*_error_rate.txt",full.names=T)

Testdata_error_files <- lapply(Testdata_filenames, fread)

# Drop first column of all files so we can combine them better

for (i in seq_along(Testdata_error_files))
{Testdata_error_files[[i]] <- select(Testdata_error_files[[i]], -1)
}

# extract the column with the error rates (per position)

Testdata_error_rates <- lapply(Testdata_error_files,function(i) i[,1]) %>% do.call(cbind,.) %>% data.frame

## rename the columns with the actual sample names

Testdata__samples <- read.table("Testdata_sample_list.txt")
Testdata__samples <- split(Testdata__samples, seq_len(nrow(Testdata__samples)))
Testdata__samples <- lapply(Testdata__samples, function(x) x[x != ""])


colnames(Testdata_error_rates) <- c(Testdata__samples)


# Switch around columns and rows and rename columns
#first add another column specifying Sample ID and Error rate 

Testdata_error_rates <- cbind(sample = "Error_rate", Testdata_error_rates)
Testdata_error_rates_long <- Testdata_error_rates %>%  pivot_longer(!sample)


Testdata_error_rates_long <- data.frame(t(Testdata_error_rates[-1]))

Testdata_error_rates_long <- cbind(rownames(Testdata_error_rates_long), data.frame(Testdata_error_rates_long, row.names=NULL))

colnames(Testdata_error_rates_long)[1] ="Sample_ID"

colnames(Testdata_error_rates_long)[2] ="Error_rate"


# remove rows with NA (inconclusive samples have NO error rate)

Testdata_error_rates_long <- Testdata_error_rates_long[!is.na(Testdata_error_rates_long$Error_rate),]

# Split 'string' in 2nd column into 367 columns (one column per position in the Amplicon sequence)

colnames <- paste("Position", 1:367) 

Testdata_error_rates_FINAL <- separate_wider_delim(Testdata_error_rates_long, cols = Error_rate, delim = ";", names = colnames)

# save the sequencing error rate overview across your dataset

write.table(Testdata_error_rates_FINAL, file = "~/PATH/Testdata_Error_rates_perSample.csv", sep = ",")

#############################################################################################################################################
########### Plot Error rate distribution within the amplicon sequence by calculating the mean sequencing error rate for each position ###
#############################################################################################################################################

# first drop the Sample_ID column (not important for the overal pattern of sequencing errors in the dataset)

Testdata_error_rates_mean <- Testdata_error_rates_FINAL[,-1] 

# make all columns numeric

Testdata_error_rates_mean <- Testdata_error_rates_mean %>% mutate_if(is.character, as.numeric)

# calculate mean and add as new row

Testdata_error_rates_mean <-  rbind(Testdata_error_rates_mean,colMeans(Testdata_error_rates_mean))

# Isolate the row with the mean seq. error rate per position for plotting

Testdata_error_rates_mean_plot <- tail(Testdata_error_rates_mean, n=1)

# rounding mean  values for a better display 
round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

Testdata_error_rates_mean_plot <- round_df(Testdata_error_rates_mean_plot, 2)

# add a column saying 'mean' for plotting

Mean <- sample('mean')

Testdata_error_rates_mean_plot <- cbind(Mean, Testdata_error_rates_mean_plot)

#switch columns and rows around for plotting

Testdata_error_rates_mean_plot <-  Testdata_error_rates_mean_plot %>%  pivot_longer(!Mean)

#label positions with a higher sequencing error rate than 1%

Testdata_error_rates_mean_plot <- Testdata_error_rates_mean_plot %>% mutate(label = ifelse(value >= 1, value, NA))

Testdata_error_rates_mean_plot$name <- factor(Testdata_error_rates_mean_plot$name, levels = Testdata_error_rates_mean_plot$name)

keeps<-c(1,50,100,150,200,250,300,350,367) # this is the length of the Prnp amplicon (367bp), adjust depending on amplcion length

Testdata_Error_rates_plot <- ggplot(Testdata_error_rates_mean_plot, aes(x = name, y = value, label= label)) + geom_bar(stat = "identity") + 
  ggtitle("Testdata Error Rates") +  xlab("Position in amplicon") + ylab("mean error rate (%)") + 
  geom_text_repel(aes(x = name, y = value, label= label),col="red",box.padding = 0.25, nudge_x = 5,nudge_y = 0) + 
  scale_x_discrete(breaks=levels(Testdata_error_rates_mean_plot$name)[keeps],labels=keeps)

Testdata_Error_rates_plot

ggsave("Testdata_meanErrorrates_perPosition.png")


########################################################################################################################################################


#######################################################################################
#### 2. Get the average error rates per sample within your data set   ################
######################################################################################

# extract the column with the error rates (per position)

Testdata_avg_error_rates <- lapply(Testdata_error_files,function(i) i[,2]) %>% do.call(cbind,.) %>% data.frame

## rename the columns with the actual sample names

colnames(Testdata_avg_error_rates) <- c(Testdata__samples)

# Switch around columns and rows and rename columns
#first add another column specifying sample and no. reads

Testdata_avg_error_rates <- cbind(sample = "Error_rate", Testdata_avg_error_rates)
Testdata_avg_error_rates_long <- Testdata_avg_error_rates %>%  pivot_longer(!sample)

# remove rows with NA (inconclusive samples have NO error rate)

Testdata_avg_error_rates_long <- Testdata_avg_error_rates_long[!is.na(Testdata_avg_error_rates_long$value),]


write.table(Testdata_avg_error_rates_long, file = "~/PATH/Testdata_AVERAGEErrorPERSAMPLE.csv", sep = ",")


########################################################################################################################################################
# Sort by average sequencing error rate to demonstrate the distribution of error rates


Testdata_avg_error_rates_sorted <- Testdata_avg_error_rates_long[order(Testdata_avg_error_rates_long$value),]

# Add new column with values 1 - 2883 and remove SAMPLE ID column (for plotting)


Error_rates_NOfilter_trim10bp_sorted  <- Error_rates_NOfilter_trim10bp_sorted  %>%
  mutate(Samples=c(1:2883))

Testdata_avg_error_rates_sorted  <- Testdata_avg_error_rates_sorted  %>%
  mutate(Samples=c(1:39))


# Plot for AVERGAE Error rate in each sample per batch

Testdata_avg_error_rates_plot <- ggplot(Testdata_avg_error_rates_sorted, aes(x = Samples, y = value)) + geom_bar(stat = "identity") + 
  ggtitle("Testdata Average Error Rate per sample") +  xlab("Samples") + ylab("Error rate (%)") + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

Testdata_avg_error_rates_plot


ggsave("Testdata_AVERAGEErrorPERSAMPLE.png")
