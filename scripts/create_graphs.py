import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def create_df(tsv, extra):
    """
    :param tsv: Path to .tsv
    :param extra: Tuple for extra column to include in dictionary. (name, tsv_column_name)
    :return: Dataframe where each row has columns "pop", "var_id", extra and "AF"
    """

    variants = pd.read_csv(tsv, sep='\t')

    dicts = []
    for index, row in variants.iterrows():
        for population in populations:
            dict = {}
            dict['Population'] = population
            dict['Allele Frequency'] = row['AC_' + population] / row['AN_' + population]
            dict[extra[0]] = row[extra[1]]
            dicts.append(dict)

    graph_inesss = pd.DataFrame(dicts)

    return graph_inesss


populations = ["eas", "afr", "amr", "asj", "sas", "nfe", "fin"]


""" SEABORN SETTINGS """

custom = ["#ec6b2d", "#00aabb","#ff5959", "#545454", "#fbae17", "#34495e", "#2ecc71"]
sns.set_palette(sns.color_palette(custom))
sns.set_style("darkgrid")


""" INESSS GRAPH """

graph_inesss = create_df('../analysis/inesss_variants.tsv', ('VAR_ID', 'VAR_ID'))

ax = sns.catplot(data=graph_inesss, x='Population', y='Allele Frequency', hue="VAR_ID", kind='bar', legend=False, aspect=2, height=9, hue_order=['1-97981343-A-C', '1-97547947-T-A' ,'1-97915614-C-T', '1-98039419-C-T'])

plt.legend(loc='upper left')
plt.show()


"""" CLIN GRAPH """

graph_clin = create_df('../analysis/clin_variants.tsv', ('CLIN_SIG', 'CLIN_SIG'))

ax = sns.catplot(data=graph_clin, x='Population', y='Allele Frequency', hue="CLIN_SIG", kind='swarm', legend=False, aspect=1.5, height=5)
plt.ylim(0,0.002)
plt.legend(loc='upper right')
plt.show()


"""" LOF GRAPH """

graph_LOF = create_df('../analysis/lof_variants.tsv', ('LOF', 'LOF'))

ax = sns.catplot(data=graph_inesss, x='POP', y='AF', hue="LOF", kind='swarm', legend=False, aspect=1.5, height=5)
plt.ylim(0,0.0005)
plt.legend(loc='upper right')
plt.show()


""" COMBINED GRAPH """

all_variants = pd.read_csv('../analysis/all_variants.tsv', sep='\t')

dicts = []
for index, row in all_variants.iterrows():
    for population in populations:
        dict = {}
        dict['Population'] = population
        dict['VAR_ID'] = row['VAR_ID']

        if row['INESSS'] == True:
            second_cat = 'INESSS'
        elif row['CLIN_SIG'] in ["Pathogenic", "Likely_pathogenic" , "Pathogenic/Likely_pathogenic"]:
            second_cat = 'CLIN_VAR'
        else:
            second_cat = 'LOF'

        dict['SECOND'] = second_cat
        dict['Allele Frequency'] = row['AC_' + population] / row['AN_' + population]
        dicts.append(dict)

graph_inesss = pd.DataFrame(dicts)

ax = sns.catplot(data=graph_inesss, x='Population', y='Allele Frequency', hue="SECOND", kind='swarm', legend=False, aspect = 1.5)
plt.legend(loc='upper right')
plt.ylim(0,0.03)
plt.show()


""" INESSS VS CLINVAR VS LOF BY POPULATION """

# Get variants
inesss_variants = pd.read_csv('../analysis/inesss_variants.tsv', sep='\t')
clin_variants = pd.read_csv('../analysis/clin_variants.tsv', sep='\t')
lof_variants = pd.read_csv('../analysis/lof_variants.tsv', sep='\t')

dicts = []

# Add each inesss variant
for index, row in inesss_variants.iterrows():
    for population in populations:
        dict = {}
        dict['Population'] = population
        dict['Source'] = row['VAR_ID']
        dict['Allele Frequency'] = row['AC_' + population] / row['AN_' + population]
        dicts.append(dict)

# Count total AF for the ClinVar variants
totals = {population: 0 for population in populations}
for index, row in clin_variants.iterrows():
    for population in populations:
        totals[population] += row['AC_' + population] / row['AN_' + population]

# Add entries
for population in populations:
    dicts.append({
        'Source': 'ClinVar',
        'Population': population,
        'Allele Frequency': totals[population]
    })

# Count total AF for LOF variants
totals = {population: 0 for population in populations}
for index, row in lof_variants.iterrows():
    for population in populations:
        totals[population] += row['AC_' + population] / row['AN_' + population]

# Add entries
for population in populations:
    dicts.append({
        'Source': 'LOF',
        'Population': population,
        'Allele Frequency': totals[population]
    })

ax = sns.catplot(data=pd.DataFrame(dicts), x='Population', y='Allele Frequency', hue='Source', kind='bar', legend=False, aspect=2)
plt.legend(loc='upper left')
plt.ylim(0,0.03)
plt.show()


""" ENTIRE POPULATION GRAPH """

dicts = []

# Add each inesss variant
for index, row in inesss_variants.iterrows():
    dict = {}
    dict['ID'] = row['VAR_ID']
    dict['Allele Frequency'] = row['AC'] / row['AN']
    dicts.append(dict)

# Get total for ClinVar variants
total = 0
for index, row in clin_variants.iterrows():
    total += (row['AC'] / row['AN'])

# Add entry
dicts.append({
    'ID': 'ClinVar',
    'AF': total
})

# Get total for LOF variants
total = 0
for index, row in lof_variants.iterrows():
    total += row['AC'] / row['AN']

# Add entry
dicts.append({
    'ID': 'LOF',
    'AF': total
})

graph_inesss = pd.DataFrame(dicts)

ax = sns.catplot(data=graph_inesss, x='ID', y='Allele Frequency', kind='bar', legend=False, aspect=2)
plt.legend(loc='upper right')
plt.ylim(0,0.015)
plt.show()