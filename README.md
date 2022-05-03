# Chameleolyser
Chameleolyser is a bioinformatics tool to identify genetic variants in homologous regions using whole-exome sequencing (WES) data. These variants remain hidden in a regular WES analysis. The current implementation of our software is hg19-based and is tested on CentOS Linux 7. However, it should run on any Linux OS.

## Installation
It is highly recommended to install all dependencies by cloning the Chameleolyser repository onto your machine. 
```
git clone https://github.com/Genome-Bioinformatics-RadboudUMC/Chameleolyser.git
cd Chameleolyser/
conda env create -f ChameleolyserEnvironment.yml
conda activate Chameleolyser
```

## Usage
### Prepare BED
The prepareBED function will download all necessary BED files. The working directory is the directory in which all intermediate and result files will be written to. Choose an existing directory for this. This option only need to be run once (also in case multiple samples are analysed in the same working directory).
```
perl Chameleolyser.pl --PrepareBED --WORKING_DIR=<WORKING_DIRECTORY>
```
### Mask reference genome
The MaskReferenceGenome function will download a copy of the hg19 reference genome (source=1000G). After completion, it will create a masked version of it. This option only need to be run once (also in case multiple samples are analysed in the same working directory).
```
perl Chameleolyser.pl --MaskReferenceGenome --WORKING_DIR=<WORKING_DIRECTORY>
```

### Generate masked alignments and raw VCF

```
perl Chameleolyser.pl --GenerateMaskedAlignmentAndVcf --WORKING_DIR=<WORKING_DIRECTORY> --SAMPLE_NAME=<SAMPLE_NAME> --ALIGNMENT_FP=<ALIGNMENT_FP> --NR_OF_THREADS=<NR_OF_THREADS>
```

