# MhGeneS

## Description

**MhGeneS** (**M**icro**h**aplotype **Gene S**creening) works in tandem with Seq2Sat to help validate microhaplotype genotyping of the coding region of genes. Within the MhGeneS pipeline the user can assess data quality in the form of sequencing error rates and read depth filter to increase the power and confidence in downstream genotype calling based on genic regions. MhGeneS can also be applied for genotype calling within non-coding regions. 

## Getting started

### Step 1. Run Seq2Sat

For detailed instructions please visit: https://github.com/ecogenomicscanada/Seq2Sat
*Please be aware that the **parameter text file** needs to be adjusted to your dataset (e.g. primers, trimming option, known positions within the amplicon)*

A test dataset with raw files which can be processed with Seq2Sat can be downloaded here:  **https://www.ecogenomicscanada.ca/wp-content/uploads/files/MhGeneS_testdata_fastq.tar.gz**

### Step 2. Microhaplotype screening and sequencing error rate exploration

Seq2Sat will output several files (for more details:https://github.com/ecogenomicscanada/Seq2Sat). Within the **MhGeneS** pipeline two output files are processed further.

#### Step 2.1. Sequencing error rate exploration

#### Step 2.2. Microhaplotype screening


## Pipeline overview
![Figure1_pipeline_V3](https://github.com/ecogenomicscanada/MhGeneS/assets/70644096/ecf0713f-85a3-4224-b093-72e0a990ee60)
