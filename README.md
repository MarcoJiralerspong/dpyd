# DPYD Loss of Function Analysis
Analysis of DPYD variants with loss of function characteristics and their frequency across different populations.

## Purpose
Use the **gnomAD** (genome and exome) and **ClinVar** VCF files to identify variants with potential loss of function. The VCFs are parsed in python using [cyvcf2](https://github.com/brentp/cyvcf2) (very fast parsing of VCF with region-queries) to find all variants which belong to the gene DPYD. These variants are then filtered to only include those with no filter and that either exhibit one of the three following characteristics:
 - One of the INESSS variants
 - Low or high confidence LOF
 - Likely pathogenic or pathogenic according to ClinVar

Once filtered, we compare the allele frequency of these variants as well as the relative frequency of the INESSS variants (compared to all the variants) for each population provided in gnomAD.

## Required Packages

 - numpy
 - pandas
 - cyvcf2
 - matplotlib
 - seaborn (for prettier graphs)

## Structure
```bash
├───data
│   ├───Put data files here
│───scripts
│   ├───Python scripts are here
│───analysis
│   ├───.tsv files are generated here
```


## Required Data

 - [gnomAD exome VCF (chromosome 1)](https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.1.vcf.bgz)
	 - [Tabix](https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.1.vcf.bgz.tbi)
 - [gnomAD genome VCF (chromosome 1)](https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.1.vcf.bgz)
	 - [Tabix](https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.1.vcf.bgz.tbi)
 - [ClinVar VCF (all chromosomes)](ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20190513.vcf.gz)
	 - [Tabix](ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20190513.vcf.gz.tbi)

## Use

 1. Create directory with `data` and `scripts` folders.
 2. Download all data files into

