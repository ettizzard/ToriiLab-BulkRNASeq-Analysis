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

### Data Accession

The corresponding ```fastq``` files were accessed and downloaded from the experiment [SRA Repository](https://www.ncbi.nlm.nih.gov/sra?term=SRP125136).

There were 9 runs in total - 3 sets of estradiol treatments with 2 runs per set and 3 sets of DMSO treatments (control) with 1 run per set. The ```fastq``` files were paired-end, so they were also split into their Read_1 and Read_2 files with a repo of 18 individual ```fastq``` files. 
```
$ prefetch SRR6300219
$ prefetch SRR6300220
$ prefetch SRR6300221
$ prefetch SRR6300222
$ prefetch SRR6300223
$ prefetch SRR6300224
$ prefetch SRR6300225
$ prefetch SRR6300226
$ prefetch SRR6300227
```

```
$ fastq-dump --outdir fastqs --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip /projects/bgmp/tizzard/bioinfo/ToriiLab-BulkRNASeq-Analysis/SRA_cache/SRR6300219/SRR6300219.sra
...
$ fastq-dump --outdir fastqs --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip /projects/bgmp/tizzard/bioinfo/ToriiLab-BulkRNASeq-Analysis/SRA_cache/SRR6300227/SRR6300227.sra
```

The 2 estradiol treatment run files were then concatenated per replicate, resulting in 12 final ```fastq``` files.

```
$ cat SRR6300222_pass_1.fastq.gz SRR6300223_pass_1.fastq.gz > Estradiol_S1_R1.fastq.gz
$ cat SRR6300222_pass_2.fastq.gz SRR6300223_pass_2.fastq.gz > Estradiol_S1_R2.fastq.gz

$ cat SRR6300224_pass_1.fastq.gz SRR6300225_pass_1.fastq.gz > Estradiol_S2_R1.fastq.gz
$ cat SRR6300224_pass_2.fastq.gz SRR6300225_pass_2.fastq.gz > Estradiol_S2_R2.fastq.gz

$ cat SRR6300226_pass_1.fastq.gz SRR6300227_pass_1.fastq.gz > Estradiol_S3_R1.fastq.gz
$ cat SRR6300226_pass_2.fastq.gz SRR6300227_pass_2.fastq.gz > Estradiol_S3_R2.fastq.gz
```

## Quality Control

FastQC is a requirement of both pipelines, so it was applied to the ```fastq``` files prior to any alignment steps of either workflow. The [FastQC results](https://github.com/ettizzard/ToriiLab-BulkRNASeq-Analysis/tree/main/fastqc) were quite positive, with all 12 files exhibiting high-quality average Phred scores >28 and no flagged low-quality sequences. Although many samples exhibit high levels of duplication, the distribution of duplication appears to be within an acceptable range and does not indicate immeidate issues with methodology. Further analysis can thus proceed.
```
fastqc *.fastq.gz
```

## STAR Genome File Accession

Prior to the alignment steps of either pipeline, an *A. thaliana* reference genome ```fasta``` and gene annotation ```gtf``` file were accessed and downloaded from the [NCBI RefSeq assembly submitted by TAIR](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001735.4/).



## Pipeline 1: STAR and DESeq2


### Reference Genome Index Generation:

After installing STAR and acquiring both the reference genome assembly ```fasta``` file and the gene annotation ```gtf```, reference genome indexing can begin.

```
$ STAR \
--runThreadN 14 \
--runMode genomeGenerate \
--genomeDir /projects/bgmp/tizzard/bioinfo/ToriiLab-BulkRNASeq-Analysis/star/genome_indices/ \
--genomeFastaFiles /projects/bgmp/tizzard/bioinfo/ToriiLab-BulkRNASeq-Analysis/star/genome_files/RefSeq/GCF_000001735.4_TAIR10.1_genomic.fasta \
--sjdbGTFfile /projects/bgmp/tizzard/bioinfo/ToriiLab-BulkRNASeq-Analysis/star/genome_files/RefSeq/genomic.gtf \
--genomeSAindexNbases 12
```

After the reference genome was successfully indexed, the ```fastq``` files can be mapped to the reference.

### Mapping and Filtering:

```
$ STAR \
--runThreadN 14 \
--runMode alignReads \
--genomeDir /projects/bgmp/tizzard/bioinfo/ToriiLab-BulkRNASeq-Analysis/star/genome_indices/ \
--readFilesCommand zcat \
--quantMode GeneCounts \
--outSAMtype BAM SortedByCoordinate \
--readFilesIn /projects/bgmp/tizzard/bioinfo/ToriiLab-BulkRNASeq-Analysis/fastqs/DMSO_S1_R1.fastq.gz /projects/bgmp/tizzard/bioinfo/ToriiLab-BulkRNASeq-Analysis/fastqs/DMSO_S1_R2.fastq.gz \
--outFileNamePrefix /projects/bgmp/tizzard/bioinfo/ToriiLab-BulkRNASeq-Analysis/star/alignment_output/
...
$ STAR \
--runThreadN 14 \
--runMode alignReads \
--genomeDir /projects/bgmp/tizzard/bioinfo/ToriiLab-BulkRNASeq-Analysis/star/genome_indices/ \
--readFilesCommand zcat \
--quantMode GeneCounts \
--outSAMtype BAM SortedByCoordinate \
--readFilesIn /projects/bgmp/tizzard/bioinfo/ToriiLab-BulkRNASeq-Analysis/fastqs/Estradiol_S3_R1.fastq.gz /projects/bgmp/tizzard/bioinfo/ToriiLab-BulkRNASeq-Analysis/fastqs/Estradiol_S3_R2.fastq.gz \
--outFileNamePrefix /projects/bgmp/tizzard/bioinfo/ToriiLab-BulkRNASeq-Analysis/star/alignment_output/
```
After examining the output of STAR, namely ```ReadsPerGene.out.tab```, we see successful mapping of reads to genes.

```
$ head ReadsPerGene.out.tab 
N_unmapped      1812843 1812843 1812843
N_multimapping  1159457 1159457 1159457
N_noFeature     134421  44907870        294181
N_ambiguous     1669996 34627   106157
AT1G01010       457     0       458
AT1G01020       762     1       761
AT1G03987       9       9       0
AT1G01030       372     1       372
AT1G01040       2260    71      2575
AT1G03993       0       307     1
```
Here, Column 1 responds to the Ensembl Gene ID, whereas Columns 2, 3, and 4 represent the number of reads mapped per gene. After consulting the [original paper](https://www.sciencedirect.com/science/article/pii/S1534580718302855?via%3Dihub) this data was generated for, it was determined that the mRNA libraries were prepared by using Illumina TruSeq Stranded mRNA. Since our reads are stranded, we are concerned with read counts within columns 3 and 4, with Col3 representing the first strand aligned and Col4 representing the second read strand aligned. However, since both columns 3 and 4 contain read counts, we must determine which column contains the majority of mapped reads.

```
$ tail -n +5 /projects/bgmp/tizzard/bioinfo/ToriiLab-BulkRNASeq-Analysis/star/alignment_output/DMSO_S1_output/ReadsPerGene.out.tab | cut -f 3 -d $'\t' | awk '{ sum += $1 } END { print sum }'
1692801

$ tail -n +5 /projects/bgmp/tizzard/bioinfo/ToriiLab-BulkRNASeq-Analysis/star/alignment_output/DMSO_S1_output/ReadsPerGene.out.tab | cut -f 4 -d $'\t' | awk '{ sum += $1 } END { print sum }'
46234960
```

On average, column 4 has 27-28x the number of reads that column 3 does. Thus, column 4 contains the vast majority of our properly aligned reads and will be used for further downstream analysis.
