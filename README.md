# MhGeneS

## Description

**MhGeneS** (**M**icro**h**aplotype **Gene S**creening) works in tandem with Seq2Sat to help validate microhaplotype genotyping of the coding region of genes. Important quality metrics affecting microhaplotype calling, can be the sequencing error rate profile related to the overlap or non-overlap of paired-end reads as well as the read depth. Within the MhGeneS pipeline we give the user the opportunity to assess data quality in the form of sequencing error rates and read depth filter to increase the power and confidence in downstream genotype calling based on genic regions. MhGeneS can also be applied for genotype calling within non-coding regions. 

## Getting started

### Step 1. Run Seq2Sat

For detailed instructions please visit: https://github.com/ecogenomicscanada/Seq2Sat
*Please be aware that the **parameter text file** needs to be adjusted to your dataset (e.g. primers, trimming option, known positions within the amplicon)*

A test dataset with raw files which can be processed with Seq2Sat can be downloaded here:  **https://www.ecogenomicscanada.ca/wp-content/uploads/files/MhGeneS_testdata_fastq.tar.gz**

### Step 2. Microhaplotype screening and sequencing error rate exploration

Seq2Sat will output several files (for more details:https://github.com/ecogenomicscanada/Seq2Sat). Within the **MhGeneS** pipeline two output files are processed further.

#### Step 2.1. Sequencing error rate exploration

The sequencing error rate profile can be explored and plotted with the provided R code (see *R_code* and *testdata_ErrorRate_files*). Based on the seq. error rate profile the user has the option to implement trimming of the 5' and 3' end of the amplicon sequence. If trimming is a sensible step to improve data quality, it can be set as a parameter in the **parameter text file** in Seq2Sat (see *Seq2Sat* folder). In the columns 4 and 5 of the **parameter text file** the user has to set the amount of basepairs (bp) to be trimmed off (4th column = 5' end; 5th column = 3' end).

#### Step 2.2. Microhaplotype screening

In order to sreen for Microhaplotypes (MH) within a dataset, specific R code is provided (see *R_code* folder). In this step all samples in a dataset are processed (see *testdata_MH_files* folder) and the information on the top (two) Microhaplotypes per sample is processed. By comparing all present MHs in a dataset a *Look-up* table is created, documenting the frequency of every detected microhaplotype within the dataset. Additioanlly, a *Dataset summary table* is created to summarize the present MHs for each samples and the according read depth.

An ***optional* read depth filter** is implemented in the R code, so the User can decide what read depth is appropriate.

## Pipeline overview
![Figure1_pipeline_V3](https://github.com/ecogenomicscanada/MhGeneS/assets/70644096/ecf0713f-85a3-4224-b093-72e0a990ee60|width=100)
