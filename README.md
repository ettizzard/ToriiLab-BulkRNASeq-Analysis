# *A. thaliana* Bulk RNA-Seq Analysis

## Deliverables

1) Process the [bulk RNA-Seq data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107018) with 2 separate software pipelines:
    1) FastQC -> STAR (alignment and feature counts) -> DESeq2 (DGE)
    2) FastQC -> STAR (alignment) -> Cufflinks (feature counts) -> Cuffdiff (DGE)

2) Report findings regarding specific classes of genes induced or repressed by induced MUTE overexpression. Compare the results of each pipeline.


<!--
could you please process the following bulk-RNA-seq data, and then provide us your pipeline? Also report the findings regarding the specific classes of genes induced or repressed by induced MUTE overexpression. You are welcome to do further elaborated analysis

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107018

These are classic bulk-RNA seq data, so does not require expert level of skills for analysis. If you could compare different suite of pipelines (as I indicated in the previous e-mail) and find out whether the results change based on the different pipelines, that would also be interested in.

FAST-QC -> HISAT/TOPHAT2/STAR -> Feature Counts -> DEseq2 (or EdgeR; we used DEseq2 for normalization)
Alternatively

FAST-QC -> HISAT/TOPHAT2/STAR ->CUFFLINKS -> CuffDiff2


-->

## Prerequisites

The corresponding ```fastq``` files were accessed and downloaded from the experiment [SRA Repository](https://www.ncbi.nlm.nih.gov/sra?term=SRP125136). There were 9 runs in total - 3 sets of estradiol treatments with 2 runs per set and 3 sets of DMSO treatments (control) with 1 run per set. The ```fastq``` files were paired-end, so they were also split into their Read_1 and Read_2 files with a repo of 18 individual ```fastq``` files. The 2 estradiol treatment run files were then concatenated per accession number, resulting in 12 final ```fastq``` files.

FastQC is a requirement of both pipelines, so it was applied to the ```fastq``` files prior to any alignment steps of either workflow. The FastQC results were quite positive, with all 18 files exhibiting high-quality Phred scores >28 and no flagged low-quality sequences. This was expected, and further data analysis can proceed without any initial filtering steps.

Prior to the alignment steps of either pipeline, an *A. thaliana* reference genome ```fasta``` and gene annotation ```gtf``` file were accessed and downloaded from the [NCBI RefSeq assembly submitted by TAIR](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001735.4/).


#### Alignment to Reference Genome with STAR

##### Genome Index Generation:
```
$ STAR \
--runThreadN 8 \
--runMode genomeGenerate \
--genomeDir /Users/evan/bioinfo/ToriiLab-BulkRNASeq-Analysis/star/genome_indices/ \
--genomeFastaFiles /Users/evan/bioinfo/ToriiLab-BulkRNASeq-Analysis/star/genome_files/RefSeq/GCF_000001735.4_TAIR10.1_genomic.fasta \
--sjdbGTFfile /Users/evan/bioinfo/ToriiLab-BulkRNASeq-Analysis/star/genome_files/RefSeq/genomic.gtf \
--genomeSAindexNbases 12
```

##### Mapping and Filtering:

```
$ STAR \
--runThreadN 8 \
--runMode alignReads \
--genomeDir /Users/evan/bioinfo/ToriiLab-BulkRNASeq-Analysis/star/genome_indices/ \
--readFilesIn /Users/evan/bioinfo/ToriiLab-BulkRNASeq-Analysis/fastqs/SRX3401147_DMSO_S1_R1.fastq /Users/evan/bioinfo/ToriiLab-BulkRNASeq-Analysis/fastqs/SRX3401147_DMSO_S1_R2.fastq \
--outFileNamePrefix /Users/evan/bioinfo/ToriiLab-BulkRNASeq-Analysis/star/alignment_output/
```



## Pipeline 1


