# Overview
Chameleolyser is a bioinformatics tool to identify genetic variants in homologous regions using whole-exome sequencing (WES) data. These variants remain hidden in a regular WES analysis. The current implementation of our software is hg19-based and is tested on CentOS Linux 7. However, it should run on any Linux OS. The only required input is a CRAM or BAM file.

# Installation
It is highly recommended to install all dependencies by cloning the Chameleolyser repository onto your machine. This procedure takes around 5 minutes.
```
git clone https://github.com/Genome-Bioinformatics-RadboudUMC/Chameleolyser.git
cd Chameleolyser/
conda env create -f ChameleolyserEnvironment.yml
conda activate Chameleolyser
```

# Usage
## Prepare BED
The prepareBED function will download all necessary BED files. The working directory is the directory in which all intermediate and result files will be written. Choose an existing directory for this. The PREFIX option can be used to indicate whether or not the names of the chromosomes start with 'chr' (i.e. NCBI reference genome) in the reference sequence that was used to generate your input CRAM/BAM. The prepareBED function only need to be run once (also in case multiple samples are analysed in the same working directory). The OMIM option can be used if only known disease genes need to be analysed.
```
perl Chameleolyser.pl --PrepareBED --WORKING_DIR=<WORKING_DIRECTORY> --PREFIX=chr --OMIM=yes
```
## Mask reference genome
The MaskReferenceGenome function will download a copy of the hg19 reference genome. After completion, it will create a masked version of it. This option only need to be run once (also in case multiple samples are analysed in the same working directory).
```
perl Chameleolyser.pl --MaskReferenceGenome --WORKING_DIR=<WORKING_DIRECTORY> --PREFIX=chr --OMIM=yes
```

## Generate masked alignments and raw VCF
This function will extract reads in the homologous regions and re-align them to the masked reference sequence. Subsequently it will call variants with a sensitive method. The sample name is an identifier of choice. The alignment filepath is the full path of the CRAM/BAM file of your sample of interest which is stored on your machine.

```
perl Chameleolyser.pl --GenerateMaskedAlignmentAndVcf --WORKING_DIR=<WORKING_DIRECTORY> --PREFIX=chr --SAMPLE_NAME=<SAMPLE_NAME> --ALIGNMENT_FP=<ALIGNMENT_FP> --NR_OF_THREADS=<NR_OF_THREADS> --OMIM=yes
```
## Filter raw variants

```
perl Chameleolyser.pl --FilterRawVariants --WORKING_DIR=<WORKING_DIRECTORY> --PREFIX=chr --SAMPLE_NAME=<SAMPLE_NAME> --OMIM=yes
```

# Demo
## Download public exome sample
Download the BAM file (and the associated index) of HG002 into a directory of choice on your system.

```
wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam
ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bai
```

# License
[GNU GPLv3](https://choosealicense.com/licenses/gpl-3.0/)

Copyright (c) 2023 Wouter Steyaert















