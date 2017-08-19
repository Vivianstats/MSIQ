MSIQ: Joint modeling of multiple RNA-seq samples for accurate isoform quantification
================
Wei Vivian Li, Anqi Zhao, Jingyi Jessica Li
2017-08-17

<!-- README.md is generated from README.Rmd. Please edit that file -->
Introduction
------------

`MSIQ` is developed for more accurate and robust isoform quantification, by integrating multiple RNA-seq samples under a Bayesian framework. Our method aims to

-   identify the consistent group of samples with homogeneous quality
-   improve isoform quantification accuracy by jointly modeling multiple RNA-seq samples with more weights on the consistent group.

Any suggestions on the package are welcome! For technical problems, please report to [Issues](https://github.com/Vivianstats/MSIQ/issues). For suggestions and comments on the method, please contact Wei (<liw@ucla.edu>) or Dr. Jessica Li (<jli@stat.ucla.edu>).

Installation
------------

The package is not on CRAN yet. For installation please use the following codes in `R`

``` r
install.packages("devtools")
library(devtools)

install_github("Vivianstats/MSIQ")
```

Usage
-----

The package needs two types of input: GTF file and BAM files.

#### GTF file

MSIQ is an annotation-based method, so the users need to supply a GTF file from which MSIQ extracts the gene structures.

#### BAM files

MSIQ starts from alignment file in bam format. Each bam file represents one RNA-seq sample. The bam files need to be sorted and indexed, and the index files should be in the same folder as their corresponding bam files. Currently, the package only supports paired-end reads.

#### Example

Suppose the genes of interest are specified in the GTF file named as "example.gtf". We would like to estimate the transcript expression of these genes based on 6 RNA-seq samples. The bam files are named as "sample1.bam", ... , "sample4.bam". Suppose all the files are in the working directory, then we can use the following example codes to run the program

``` r
library(MSIQ)
D = 4 # number of samples
gtf_path = "example.gtf"
bam_path = c("sample1.bam", "sample2.bam", "sample3.bam", "sample4.bam")
ncores = 5
result = msiq(D, gtf_path, bam_path, ncores)
```

`MSIQ` summarizes the estimation results in MSIQ\_result/transcript\_summary.txt. The text file includes the following the columns:

-   transcriptID. The transcript ID extracted from the GTF file.
-   geneID. The gene ID extracted from the GTF file.
-   FPKM\_1, ..., FPKM\_D. The FPKM values of the corresponding transcript in samples 1, ..., *D*.
-   sample1, ..., sampleD. The estimated indicator variables of samples 1, ..., *D*. 1 means that the sample belongs to the consistent group; 0 means that the sample is not in the consistent group.

Below is an example output for one gene:

``` r
transcriptID geneID frac FPKM_1 FPKM_2 FPKM_3 FPKM_4 sample1 sample2 sample3 sample4
ENST00000366899 ENSG00000143507 0.9866 0.5164 0.6148 0.1041 0.1691 1 1 0 0 
ENST00000468085 ENSG00000143507 0.0028 0.0025 0.003 5e-04 9e-04 1 1 0 0
ENST00000477026 ENSG00000143507 0.0096 0.015 0.0178 0.0055 0.015 1 1 1 0
ENST00000494642 ENSG00000143507 0.001 0.0021 0.0025 4e-04 8e-04 1 1 1 0
```
