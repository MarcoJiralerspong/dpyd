import pandas as pd
import requests
from cyvcf2 import VCF


def get_gene_clin(variant):
    """
    :param variant: cyvcf2 variant object from ClinVar
    :return: Gene corresponding to variant, False if not found
    """
    if variant.INFO.get('GENEINFO') is None:
        return False

    return variant.INFO.get('GENEINFO').split(':')[0]


def get_gene_gnomad(variant):
    """
    :param variant: cyvcf2 variant object from gnomAD
    :return: Gene corresponding to variant
    """
    return variant.INFO.get('vep').split('|')[3]


def create_var_dict(var_list):
    """
    :param var_list: List of variants
    :return: Dictionary with following key/value pair: <CHROM>-<POS>-<REF>-<ALT>: variant
    """
    dict = {}

    for variant in var_list:
        var_id = str(variant.CHROM) + "-" + str(variant.POS) + "-" + variant.REF + "-" + ''.join(variant.ALT)

        dict[var_id] = variant

    return dict


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def get_transcript(ids):
    """
    :param id: List of variant IDs
    :return: Dict of var_id -> transcript from variantvalidator.org (API)
    """

    # Make request
    url = "https://rest.variantvalidator.org/variantformatter/GRCh37/%s/all/None/False" % ('|').join(ids)
    response = requests.get(url)

    dict = response.json()
    transcript_dict = {}

    # Create dictionary
    for id in ids:
        id_dict = dict[id][id]['hgvs_t_and_p']
        transcript_dict[id] = id_dict[list(id_dict.keys())[0]]['t_hgvs']

    return transcript_dict


def add_transcript_column(df):
    """
    :param df: Dataframe to add column to with ('VAR_ID') column
    :return: Dataframe with added column
    """
    transcripts = {}

    # Splits the list of variants into lists of 10 for VariantValidator (recommends 10 at a time)
    for id_list in list(chunks(df['VAR_ID'].tolist(), 10)):

        # Get corresponding dictionary to list
        transcript_dict = get_transcript(id_list)

        # Update transcripts dict
        for id in id_list:
            transcripts[id] = transcript_dict[id]

    # Add column to df
    df = df.assign(TRANSCRIPT = pd.Series(transcripts))

    return df


def create_row_dict(gnomad_dict, clinvar_dict):
    """
    Function to construct a row of the final .tsv
    :param gnomad_dict: A dictionary of var_id -> variant from gnomad
    :param clinvar_dict: A dictionary of var_id -> variant from clinvar
    :return: Dictionary with aggregated counts between what was already in output dict
    """

    output_dict = {}

    for id, variant in gnomad_dict.items():

        # Basic fields
        output_dict[id] = {
            'VAR_ID': id,
            'CHR': variant.POS,
            'REF': variant.REF,
            'ALT': ''.join(variant.ALT),  # Formatted as list for some reason
            'QUAL': variant.QUAL,
            'FILTER': 'None' if variant.FILTER == None else variant.FILTER,
            'AC': variant.INFO.get('AC'),
            'AN': variant.INFO.get('AN'),
            'NHOMALT': variant.INFO.get('nhomalt'),
            'LOF': variant.INFO.get('vep').split('|')[64],  # LOF in vep field
            'INESSS': 'True' if id in inesss_var_id else 'False'
        }

        # If there is a Clin_SIG from ClinVar, use that one, else use provided one (always none so far)
        if id in clinvar_dict.keys():
            output_dict[id]['CLIN_SIG'] = clinvar_dict[id].INFO.get('CLNSIG')
        else:
            output_dict[id]['CLIN_SIG'] = 'NA'

        # Add AC, AN for each population
        for population in populations:
            output_dict[id]['AC' + '_' + population] = variant.INFO.get('AC' + '_' + population)
            output_dict[id]['AN' + '_' + population] = variant.INFO.get('AN' + '_' + population)
            output_dict[id]['NHOMALT' + '_' + population] = variant.INFO.get('nhomalt' + '_' + population)


    return output_dict


""" GET DPYD VARIANTS """

# Create lists of dpyd variants from ClinVar and gnomAD vcfs
clin_vcf = VCF('../data/clinvar_20190513.vcf.gz')
gnomad_exome_vcf = VCF('../data/gnomad.exomes.r2.1.1.sites.1.vcf.bgz')

clin_dpyd_variants = []
for variant in clin_vcf('1:97443300-98486606'):
    if get_gene_clin(variant) == 'DPYD':
        clin_dpyd_variants.append(variant)

gnomad_exome_dpyd_variants = []
for variant in gnomad_exome_vcf('1:97443300-98486606'):
    if get_gene_gnomad(variant) == 'DPYD':
        gnomad_exome_dpyd_variants.append(variant)

# Create dictionaries from list using var_id as key
clin_dict = create_var_dict(clin_dpyd_variants)
gnomad_exome_dict = create_var_dict(gnomad_exome_dpyd_variants)


""" FORMAT DICTIONARIES """
inesss_var_id = ['1-97915614-C-T',
                 '1-97547947-T-A',
                 '1-97981343-A-C',
                 '1-98039419-C-T'] # Variants recommended for testing by INESSS

populations = ["eas", "afr", "amr", "asj", "sas", "nfe", "fin"] # Populations available in gnomAD


# Create dictionary with relevant fields for each variant found in gnomad
clean_dict = create_row_dict(gnomad_exome_dict, clin_dict) # Update entries from exome

pre_filter_variants = pd.DataFrame.from_dict(clean_dict, orient='index')

# Export pre_filter version to .tsv
export_tsv = pre_filter_variants.to_csv (r'../analysis/clean_gnomad.tsv', index = None, header=True, sep='\t')


""" SPLIT DATA """

LOFs = ["HC", "LC"]
CLIN_SIGs = ["Pathogenic", "Likely_pathogenic" , "Pathogenic/Likely_pathogenic"]

# Only look through those with no filter
filtered_variants = pre_filter_variants.loc[(pre_filter_variants['FILTER'] == 'None')]


# Only INESSS variants
inesss_df = filtered_variants.loc[filtered_variants['INESSS'] == 'True']
inesss_df = add_transcript_column(inesss_df)

# Only likely or pathogenic variants according to ClinVar (not INESSS)
clin_df = filtered_variants.loc[filtered_variants['CLIN_SIG'].isin(CLIN_SIGs) & (filtered_variants['INESSS'] != 'True')]
clin_df = add_transcript_column(clin_df)

# Only HC or LC LOF not in first 2
lof_df = filtered_variants.loc[filtered_variants['LOF'].isin(LOFs) &
                                 (filtered_variants['INESSS'] != 'True') &
                                 ~filtered_variants['CLIN_SIG'].isin(CLIN_SIGs)
]
lof_df = add_transcript_column(lof_df)

# Those in either of the first 3
all_df = filtered_variants.loc[filtered_variants['LOF'].isin(LOFs) |
                                 (filtered_variants['INESSS'] == 'True') |
                                 filtered_variants['CLIN_SIG'].isin(CLIN_SIGs)
]
all_df = add_transcript_column(all_df)


""" EXPORT TO .TSV """

inesss_df.to_csv (r'../analysis/inesss_variants.tsv', index = None, header=True, sep='\t')
clin_df.to_csv (r'../analysis/clin_variants.tsv', index = None, header=True, sep='\t')
lof_df.to_csv (r'../analysis/lof_variants.tsv', index = None, header=True, sep='\t')
all_df.to_csv (r'../analysis/all_variants.tsv', index = None, header=True, sep='\t')



