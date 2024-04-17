l########################################################################################
########################################################################################
##### This R script processes Microhaplotypes (produced by Seq2Sat) for each  ##########
#####  sample and creates a 'look-up' and 'dataset summary' table about       ##########
######      Microhaplotypes within the entire dataset.                    ##############
########################################################################################
########################################################################################

setwd("~/YOUR_WORKING_DIRECTORY/")

library(data.table)
library(dplyr)


##############################################################################
#### 1. Create empty dataframes to be filled, thresholds to be set  ##########
##############################################################################

# 1.1 This file shows why samples did not meet the criteria to be processed in run and why they are excluded
excludefileConn<-file("Testdataset_excluded_Samples_MH_screening.txt", "a")

# 1.2 Depth-of-read filter: less than # of reads will be skipped
####   5 is the default for Seq2Sat and was therefore considered to be  the 'unfiltered' setting.
####   in the MhGeneS pipelin we tested read filters of 10, 20 and 50.

num_reads_allowed = 5

# 1.3 create the empty dataframes for the 'look-up' and dataset summary table

tbl_lookup <- data.frame(
                 name_h=character(), 
                 prnp_string=character(),
                 sample_freq=numeric(),
                 dataset_freq=numeric(),
                 stringsAsFactors=FALSE)

tbl_history <- data.frame(sampleID=numeric(),
                          mh_a=character(), 
                          mh_b=character(),
                          total_reads=numeric(),
                          stringsAsFactors=FALSE)

# 1.4 this is to grab the element after last slash using gsub and grabbing sample ID (parse underscore)
#### that is specific to the sample names in the testdataset (e.g. 49934_S9_L001_snps_haplotype.txt). The first 5 digit number is the sample ID
#### this has to be adapted to the sample name structure

extract_sample_id <- function(fileN) {
  getLastElement<- gsub("^.*/", "", fileN)
  parse_first_element <- strsplit(getLastElement, '_')[[1]]
  
  return(parse_first_element[1])
}


####################################################################################################################
#### 2. Functions for parsing through all samples (sample by sample), comparing Microhaplotypes,   #################
####    getting the frequencies of the Microhaplotypes (MHs) within the dataset and populating    ##################
#####                the 'look-up' table and the 'dataset summary' table                         ###################
####################################################################################################################

# 2.1 Function for calculating the sample and dataset frequency.
#### Sample freq= how many samples do show a certain MH 
#### Dataset freq= how many times a MH is present in the dataset

update_frequency_count <- function(mha, mhb)
  
{
  #########################################################
  ###Frequency for sample freq####
  #homozygote counts as one increment for that MH in sample_count column
  if(mha == mhb){
    #print(paste ('match: so homozygote:' , mha , 'and' , mhb, sep = " "))
    hap_value <- tbl_lookup$sample_freq[tbl_lookup$name_h == mha] 
    tbl_lookup$sample_freq[tbl_lookup$name_h == mha] <<- 1 + hap_value
    #second variable does not matter since same MH present and will increment by one
  }else{
    #print(paste ('non homozygote:' , mha , 'and' , mhb, sep = " "))
    hap_value_a <- tbl_lookup$sample_freq[tbl_lookup$name_h == mha] 
    tbl_lookup$sample_freq[tbl_lookup$name_h == mha] <<- 1 + hap_value_a
    
    hap_value_b <- tbl_lookup$sample_freq[tbl_lookup$name_h == mhb] 
    tbl_lookup$sample_freq[tbl_lookup$name_h == mhb] <<- 1 + hap_value_b
    
  }
  #########################################################
  ###Frequency for dataset freq####
  dataset_hap_value_a <- tbl_lookup$dataset_freq[tbl_lookup$name_h == mha] 
  tbl_lookup$dataset_freq[tbl_lookup$name_h == mha] <<- 1 + dataset_hap_value_a
  
  dataset_hap_value_b <- tbl_lookup$dataset_freq[tbl_lookup$name_h == mhb] 
  tbl_lookup$dataset_freq[tbl_lookup$name_h == mhb] <<- 1 + dataset_hap_value_b
  
  #########################################################
}

# 2.2 Processing the microhaplotypes for each sample and assessing the zygoisty for each sample
#### this function starts with the first Microhaplotype string in a sample file, documents it and compares it to the other Microhaplotypes within the
#### 'look-up- table. It will then either give it a new name or  assign the an exisiting name of the same MH string.

process_zygosity_sample <- function(sequence1, sequence2, zygosity, sampleId, totalReadsVal)
{
  #processing the 1st sequence string
  if(!is.na(sequence1)){
    #print(paste0("Sample ID: ", sampleId," has not been seen before and will continue processing:" ))
    
    #will return the name_h value and NA if the string has not been seen
    seq_name1 <-tbl_lookup %>% filter(prnp_string == sequence1) %>% select(name_h)
    
    #check if sequence is found in look-up table
    #below if NA, then sequence has not been seen, so it needs to be added and create new Microhaplotype 'name_h' name (eg. MH1)
    if(nrow(seq_name1) == 0){
      
      print('Sequence 1 has not been seen and will be entered in lookup table')
      
      name_h_tempa <- paste('MH', nrow(tbl_lookup) + 1, sep="")
      tmp <- data.frame(name_h = name_h_tempa, prnp_string = sequence1, sample_freq = 0, dataset_freq= 0)
      tbl_lookup <<- rbind(tbl_lookup, tmp)
      
      mha_int_value <- as.numeric(gsub(".*?([0-9]+).*", "\\1", name_h_tempa))
      
    }else{
      print('Sequence 1 has been seen and will be put in history')
      mha_int_value <- as.numeric(gsub(".*?([0-9]+).*", "\\1", seq_name1))
      name_h_tempa <- seq_name1$name_h
    }
    
  }else{
    print('NA found in seq1, and will not be inserted')
    NA_Sequence_a <- TRUE
    
    writeLines(c("NA found in seq1, and will not be insertedd: "), excludefileConn)
    writeLines(c(sampleId), excludefileConn)
  }

#### this function will ONLY be active IF the second MH string for the sample is different to the first (indicating heterozygosity)
  
  #only process second string if heterozygous (= expressing 2 different microhaplotypes)
  if(zygosity == 'heterozygote'){
    #processing 2nd sequence string
    if(!is.na(sequence2)){
      seq_name2 <-tbl_lookup %>% filter(prnp_string == sequence2) %>% select(name_h)
      
      if(nrow(seq_name2) == 0){
        
        print('Sequence 2 has not been seen and will be entered in lookup table')
        
        name_h_tempb <- paste('MH', nrow(tbl_lookup) + 1, sep="")
        tmp <- data.frame(name_h = name_h_tempb, prnp_string = sequence2, sample_freq = 0, dataset_freq= 0)
        tbl_lookup <<- rbind(tbl_lookup, tmp)
        
        mhb_int_value <- as.numeric(gsub(".*?([0-9]+).*", "\\1", name_h_tempb))
        
      }else{
        print('Sequence 2 has been seen and will be put in history')
        mhb_int_value <- as.numeric(gsub(".*?([0-9]+).*", "\\1", seq_name2))
        name_h_tempb <- seq_name2$name_h
      }
    }else{
      print('NA found in seq2, and will not be inserted')
      NA_Sequence_b <- TRUE
      
      writeLines(c("NA found in seq2, and will not be insertedd: "), excludefileConn)
      writeLines(c(sampleId), excludefileConn)
    }
  }
  
  ##enter into history table
  
  if(zygosity == 'heterozygote'){
    if(mha_int_value < mhb_int_value){
      tmpa <- data.frame(sampleID = sampleId, mh_a = name_h_tempa, mh_b = name_h_tempb, total_reads= totalReadsVal)
      tbl_history <<- rbind(tbl_history, tmpa)
    }else{
      tmpb <- data.frame(sampleID = sampleId, mh_a = name_h_tempb, mh_b = name_h_tempa, total_reads= totalReadsVal)
      tbl_history <<- rbind(tbl_history, tmpb)
    }
  }else if(zygosity == 'homozygote'){
    tmpa <- data.frame(sampleID = sampleId, mh_a = name_h_tempa, mh_b = name_h_tempa, total_reads= totalReadsVal)
    tbl_history <<- rbind(tbl_history, tmpa)
    name_h_tempb = name_h_tempa; #makes it so it copies the exact name before running it in update_frequency_count since it is a homozygote
  }
  
  update_frequency_count(name_h_tempa, name_h_tempb)
  
}#end of function


# 2.3 Checking function for 3things:
#### 1. If sample files have information (empty files indicate a lack of data for Seq2Sat to work, most likely due to sequencing issues)
#### 2. If a sample has any microhaplotype conclusive (Conclusive column shows 2 Y if heterozygous, 1 Y if homozygous and 2 N if 'inconclusive)
#### 3. the number of reads for the Microhaplotypes and apply the read depth filter (from above)
####### this functions first checks if a sample is heterozygous (= both rows in the sample files have a Y in the conclusive column) and adds their
####### numreads vlaues (=number of microhaplotype reads) to populate the 'look-up' and 'dataset summary' table.
####### IF the sample is homozygous, then it only takes the first row into account, since the numReads value is already the number of total reads for this
####### microhaploytpe.


process_sample <- function(MH_files, sampleId)
{
    print(paste0("Processing Sample: ", sampleId))
    
    ##checks for empty rows in file
    if(nrow(MH_files[[1]]) != 0){
      
      microHaplo <- MH_files[[1]][["Sequence"]]
      
      conclusive <- MH_files[[1]][["Conclusive"]]
      numReads <- MH_files[[1]][["NumReads"]]
     
      totalReadsVal <- 0
      
      zygosity_derived <- ''
      
      #if both rows are Y for conclusive then add the number of reads from both (heterozygous)
      #else -> check if 2nd row is an N or NA then is considered homozygous. Take numReads from 1st row only
      if(conclusive[1] == 'Y' && conclusive[2] == 'Y'){
        totalReadsVal <-  numReads[1] + numReads[2] 
        zygosity_derived <- 'heterozygote'
      }
      else if(conclusive[1] == 'Y'){
        if(conclusive[2] == 'N' || is.na(conclusive) ){
          totalReadsVal <-  numReads[1]
          zygosity_derived <- 'homozygote'
          }
      }
      
      #print(totalReadsVal)
      
      #CHECK - make sure same sampleID is not proccess
      sampleID_found <-tbl_history %>% filter(sampleId == sampleID) %>% select(sampleID)
     
    
    if((totalReadsVal >= num_reads_allowed) & (!is.na(totalReadsVal))) {
      
      if(is.na(sampleID_found[1,])){
          if(zygosity_derived == "heterozygote" || zygosity_derived == 'homozygote') 
            process_zygosity_sample(microHaplo[1], microHaplo[2], zygosity_derived, sampleId, totalReadsVal)
        }
        else{
          print('Sample ID has been seen before, so will be skipped')
          writeLines(c("Sample ID has been seen before, so will be skipped: "), excludefileConn)
          writeLines(c(sampleId), excludefileConn)
        }
          
         
    }else{
        print('Does not need meet totalReadsVal min requirments and will be skipped')
        writeLines(c("Does not need meet totalReadsVal min requirments and will be skipped: "), excludefileConn)
        writeLines(c(sampleId), excludefileConn)
        }
  
      
    }else{
      print("Empty File, will be skipped")
      writeLines(c("Empty File, will be skipped: "), excludefileConn)
      writeLines(c(sampleId), excludefileConn)
    }
    
    print('---------------------------------------------------------')

  

}#end of function


#check if current sample being looked at is in samples_list and process if TRUE
start_run <- function(batch, sample_lists_select)
{
  for (k in 1:length(batch)){
    
    sampleId <- extract_sample_id(batch[k])
    MH_files <- lapply(batch[k], fread)
    # #print(sample[1,])
    
    if(sample_lists_select){
      if(any(samples_list == sampleId)){
        print(paste0("Processing Sample: ", sampleId))
        process_sample(MH_files, sampleId)
      }
    }else{
      process_sample(MH_files, sampleId)
    }
    
    
  }
}#end of function


##############################################################################
#### 3. Run the Microhaplotype screening on your dataset       ###############
##############################################################################

# 3.1 set the folder with the Microhaplotype output files from Seq2at

batchList <- c('testdata_MH_files')

### OPTIONAL
# IF you want to subset your dataset and only run certain samples (then 'samples_lists = TRUE)
#--------------------------------------------
sample_lists = FALSE

#Only samples in the text file will be ran.
samples_list <- read.table('samples_list_examples.txt', header = FALSE, sep = "")
#---------------------------------------------

# 3.2 process each sample within the dataset

for ( i in 1:length(batchList))
{
  print(batchList[i])
  batches <- list.files(path = batchList[i], pattern="*_haplotype.txt",full.names=T)
  start_run(batches, sample_lists)
  
  
}

# 3.3 closes the connection to the file which lists the excluded samples
close(excludefileConn)


##############################################################################
#### 4. Save the 'look-up' and 'dataset summary' tables        ###############
##############################################################################

write.table(tbl_history, file = "~/PATH/Testdataset_history.csv", sep = ",")

write.table(tbl_lookup, file = "~/PATH/Testdataset_lookup.csv", sep = ",")

