########################################################################################
########################################################################################
##### This script to plotb the read depth distribution in the dataset   ###############
########################################################################################
########################################################################################

setwd("~/YOUR_WORKING_DIRECTORY/")


library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)
library("ggrepel") 


# List all the output files containing the Microhaplotype information of the test dataset

Testdata_filenames<- list.files(path = "~/PATH/testdata_MH_files/", pattern="*_haplotype.txt",full.names=T)

Testdata_MH_files <- lapply(Testdata_filenames, fread)

# Drop first column of all files so we can combine them better (first column just describes that its PRNP!)

for (i in seq_along(Testdata_MH_files))
{Testdata_MH_files[[i]] <- select(Testdata_MH_files[[i]], -1)
}

# extract the column with the number of total reads of the Microhaplotypes

Testdata_no_reads <- lapply(Testdata_MH_files,function(i) i[,5]) %>% do.call(cbind,.) %>% data.frame

## rename the column names with the actual sample names

Testdata_samples <- read.table("Testdata_sample_list.txt")
Testdata_samples <- split(Testdata_samples, seq_len(nrow(Testdata_samples)))
Testdata_samples <- lapply(Testdata_samples, function(x) x[x != ""])


colnames(Testdata_no_reads) <- c(Testdata_samples)


#extract first row with all column values

Testdata_no_reads_plot <- head(Testdata_no_reads,1)

############################
## Plot the distribution of number of total reads among data set 

#switch columns and rows around for plotting the distribution of no. of reads

Testdata_reads_long <-  Testdata_no_reads_plot %>%  pivot_longer(!sample)


# remove Samples with missing read depth information (=NA)

Testdata_reads_long <- Testdata_reads_long[!is.na(Testdata_reads_long$value),]

#this is the recommended read depth filtering in the MhGeneS manuscript
aline <- 20


Testdata_reads_Plot <- ggplot(Testdata_reads_long, aes(x = name, y = value)) + geom_point(size = 0.5)  + ggtitle("Testdata- Read depth distribution") +
  xlab("Samples") + ylab("number of reads") + geom_hline(yintercept = aline, color="red") +theme(axis.text.x=element_blank(),
                                                                                                 axis.ticks.x=element_blank(),panel.background = element_blank())

Testdata_reads_Plot <- Testdata_reads_Plot + geom_text(aes( 0, aline, label = aline, hjust = -0.3, vjust= -0.6), size = 3, color="red")

Testdata_reads_Plot

ggsave("Testdata_read_depth_distribution.png")


