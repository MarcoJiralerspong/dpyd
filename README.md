# DPYD Loss of Function Analysis
Analysis of DPYD variants with loss of function characteristics and their frequency across different populations.

## Purpose
Use the **gnomAD** (genome and exome) and **ClinVar** VCF files to identify variants with potential loss of function. The VCFs are parsed in python using [cyvcf2](https://github.com/brentp/cyvcf2) (very fast parsing of VCF with region-queries) to find all variants which belong to the gene DPYD. These variants are then filtered to only include those with no filter and that either exhibit one of the three following characteristics:
 - One of the INESSS variants
 - No function or decreased function according to CPIC
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
│───analysis
│   ├───.tsv files are generated here
├───data
│   ├───Put data files here
│───scripts
│   ├───Python scripts are here
```


## Required Data

 - [gnomAD exome VCF (chromosome 1)](https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.1.vcf.bgz)
	 - [Tabix](https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.1.vcf.bgz.tbi)
 - [gnomAD genome VCF (chromosome 1)](https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.1.vcf.bgz)
	 - [Tabix](https://storage.googleapis.com/gnomad-public/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.1.vcf.bgz.tbi)
 - [ClinVar VCF (all chromosomes)](ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20190513.vcf.gz)
	 - [Tabix](ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_20190513.vcf.gz.tbi)
 - CPIC recommendations, reformatted version included in `data` directory (original from [here](https://api.pharmgkb.org/v1/download/file/attachment/DPYD_Allele_Functionality_Table.xlsx), downloaded on May 24 2019)

## Use

 1. Navigate in your command line to directory where you want to put project.
 2. Run `git clone https://github.com/MarcoJiralerspong/dpyd.git`
 3. `cd dpyd` and create folder called `analysis`
 4. Download all data files into `data` folder, for example with `wget url` (see structure above).
 5. Navigate to the `scripts` folder and run `python create_tsvs.py` to create .tsv files and then run `python create_graphs.py` to create the graph.
 6. All results will be placed in the `analysis` folder.

## Generated files
 - `clean_gnomad.tsv`: All DPYD variants with their associated ID, clinical signature, allele count, etc.
 - `inesss_variants.tsv`: Same fields as above but only INESSS variants.
 - `clin_variants.tsv`: Only variants that have a "Likely pathogenic", "Likely pathogenic/pathogenic" or "Pathogenic" CLIN SIG and that are not in INESSS.
- `cpic_variants.tsv`: Only variants that have a "No function" or "Decreased" rating by CPIC
 - `lof_variants.tsv`: Only variants that are not in the previous 2 files and have either a "HC" or "LC" LOF.
 - `all_variants.tsv`: The union of the variants from the previous 4 files.
 - `all_cpic_variants.tsv`: Same as above but also includes CPIC variants rated as "Normal".

## Process

We begin by downloading our data which comes from 3 sources:
- **gnomAD**: It provides us with the allele counts/homozygote counts of variations, separated in to 7 different populations. For each variant, it also specifies the confidence that there is some Loss of Function (LOF) as well as the Filter. This information is available in a VCF file.
- **ClinVar**: Has the clinical significance of certain variants available once again in a VCF file.
- **CPIC**: Has a list of variants for which fluoropyrimidines might be harmful. Gives us the Allele Functional Status for these variants which is either normal, decreased or no function.

We then combine all 3 data sources together in to a central .tsv file (similar to a .csv except uses tabs instead of commas) where for each variant, we include the relevant fields from each data source.

After, we separate the variants into 4 independent categories (i.e. we exclude from category 2 the variants in category 1, we exclude from category 3 the variants in category 1 and 2, etc.):
- **INESSS**: The 4 variants that INESSS recommends that we screen.
- **CPIC**: Variants with either "Decreased" or "No function" according to the CPIC.
- **ClinVar**: Variants with either "Likely_pathogenic", "Pathogenic/Likely_pathogenic" or "Pathogenic" clinical significance.
- **LOF**: Variants with either "HC" or "LC" as LOF.

With the separated data, we can then compare and analyse the relative allele frequencies of different categories for different populations. Specifically, we can examine with plots how effective the 4 variants chosen by INESSS are at dealing with all cases of variants causing defects.


## Implementation
 - Download the .vcf.gz files and their associated tabix specified above.
 - Using a region query around the areas specified [here](https://gnomad.broadinstitute.org/gene/ENSG00000188641) with padding of 100 000 on both sides, get all variants associated to the gene DPYD.
   - The gene of the variant is obtained in the INFO.GENEINFO field (1st position) of the VCF file for ClinVar
   - For gnomAD, it is in the INFO.vep field (4th position)
 - Create dataframe with the following fields:
    - `VAR_ID`: CHR-POS-REF-ALT.
    - `CHR`: CHR from gnomAD.
    - `REF`: REF from gnomAD.
    - `ALT`: ALT from gnomAD.
    - `QUAL`: QUAL from gnomAD.
    - `FILTER`: FILTER from gnomAD.
    - `AC`: AC from gnomAD.
    - `NHOMALT`: nhomalt from gnomAD.
    - `LOF`: LOF from gnomAD.
    - `INESSS`: `'True'` if VAR_ID in `['1-97915614-C-T', '1-97547947-T-A', '1-97981343-A-C', '1-98039419-C-T']` else `'False'`
    - `ALL_FUNCT-STATUS`: Status provided by CPIC if found, else `'NA'`
    - `CLIN_SIG`: CLNSIG from ClinVar if found, else `'NA'`.
    - For each population in `["eas", "afr", "amr", "asj", "sas", "nfe", "fin"]`, `AC_pop`, `AN_pop` and `NHOMALT_pop` from gnomAD.
    - `TRANSCRIPT`: Using the VAR_ID, query the variantvalidator.org API and get the (:c.) transcript.
 - Filter the dataframe by only including variants whose `FILTER` is `NONE`.
 - Split the filtered dataframe into four distinct dataframes:
    - INESSS variants which only includes those specified by INESSS.
    - CPIC variants which have "No function" or "Decreased" for `ALL_FUNCT_STATUS`
    - CLINVAR variants which have a "Likely pathogenic", "Likely pathogenic/pathogenic" or "Pathogenic" `CLIN SIG`.
    - LOF variants which have "HC" or "LC" for `LOF`.
 - Export the dataframes as .tsv files.


## Clean CPIC

Unfortunately, the data from CPIC is given in an Excel sheet with poor formatting for parsing. To allow for parsing, some manual formatting is necessary. We create a clean .csv with the important fields from CPIC as follows:
- Create new spreadsheet
- Copy over column names
- Copy over rows which are either "Strong Evidence supporting function", "Moderate Evidence supporting function (in vitro and clinical/ex vivo data)", "In vitro data only and/or limited clinical/ex vivo data"
- Delete letters which are used for footnotes in rsID column
- Delete row rs115632870 and rs6668296 since they are linked to previous variant which probably causes decrease in function (as per footnote)
- Transform row with 2 rsIDs (rs1801267 and rs1801265) into 2 rows, using the value before the comma for the first row and the value after the comma for the second row (for fields without comma, put same in both)
- For the row with HapB3, only keep the middle variant, deleting the information to the left and right of the commas
- Delete rows rs143154602, rs72549303, rs72549309, rs1801268, rs111858276, rs137999090 and rs72547601 since they are not in gnomAD.
- Change rs72549309 to rs539032572 (version compatible with gnomAD)
- Save as .tsv (to do so in excel, save as .txt and then change the file extension to .tsv)

