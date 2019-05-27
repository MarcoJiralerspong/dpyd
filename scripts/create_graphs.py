import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def create_basic_df(variants, extra):
    """
    :param tsv: Df from .tsv
    :param extra: Tuple for extra column to include in dictionary. (name, tsv_column_name)
    :return: Dataframe where each row has columns "pop", "var_id", extra and "AF"
    """

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


def create_cat_df(categories):
    """
    :param categories: List of categories as ([tsv1, tsv2], category_name)
    :return: Concatenated dataframe of all entries of all categories, with category field corresponding to category name
    """

    df = pd.DataFrame()

    for category in categories:

        for source in category[0]:
            variants = pd.read_csv(source, sep='\t')
            variants['CATEGORY'] = category[1]
            df = df.append(variants)

    return df


def create_sum_df(df, population):
    """
    :param df:
    :return:
    """

    pop_totals = {}
    dicts = []

    if population:
        for population in populations:

            totals = {category: 0 for category in df.CATEGORY.unique()}

            for index, row in df.iterrows():
                totals[row['CATEGORY']] += row['AC_' + population] / row['AN_' + population]

            pop_totals[population] = totals

        for population in populations:
            for category in df.CATEGORY.unique():
                dicts.append({
                    'CATEGORY': category,
                    'Population': population,
                    'Allele Frequency': pop_totals[population][category]
                })

    else:
        totals = {category: 0 for category in df.CATEGORY.unique()}

        for index, row in df.iterrows():
            totals[row['CATEGORY']] += row['AC'] / row['AN']

        for category in df.CATEGORY.unique():
            dicts.append({
                'CATEGORY': category,
                'Allele Frequency': totals[category]
            })

    return pd.DataFrame(dicts)


populations = ["eas", "afr", "amr", "asj", "sas", "nfe", "fin"]

""" SEABORN SETTINGS """

custom = ["#ec6b2d", "#00aabb", "#ff5959", "#545454", "#fbae17", "#34495e", "#2ecc71"]
sns.set_palette(sns.color_palette(custom))
sns.set_style("darkgrid")

""" INESSS GRAPH """

graph_inesss = create_basic_df(pd.read_csv('../analysis/inesss_variants.tsv', sep='\t'), ('VAR_ID', 'VAR_ID'))

ax = sns.catplot(data=graph_inesss, x='Population', y='Allele Frequency', hue="VAR_ID", kind='swarm', legend=False,
                 aspect=1.5, height=5,
                 hue_order=['1-97981343-A-C', '1-97547947-T-A', '1-97915614-C-T', '1-98039419-C-T'])

plt.legend(loc='upper left')
plt.ylim(0, 0.025)
plt.show()

""" CPIC GRAPH """
graph_cpic = create_basic_df(pd.read_csv('../analysis/cpic_variants.tsv', sep='\t'),
                             ('ALL_FUNCT_STATUS', 'ALL_FUNCT_STATUS'))

ax = sns.catplot(data=graph_cpic, x='Population', y='Allele Frequency', hue="ALL_FUNCT_STATUS", kind='swarm',
                 legend=False, aspect=1.5, height=5)

plt.ylim(0, 0.025)
plt.legend(loc='upper left')
plt.show()

"""" CLIN GRAPH """

graph_clin = create_basic_df(pd.read_csv('../analysis/clin_variants.tsv', sep='\t'), ('CLIN_SIG', 'CLIN_SIG'))

ax = sns.catplot(data=graph_clin, x='Population', y='Allele Frequency', hue="CLIN_SIG", kind='swarm', legend=False,
                 aspect=1.5, height=5)
plt.ylim(0, 0.002)
plt.legend(loc='upper right')
plt.show()

"""" LOF GRAPH """

graph_LOF = create_basic_df(pd.read_csv('../analysis/lof_variants.tsv', sep='\t'), ('LOF', 'LOF'))

ax = sns.catplot(data=graph_LOF, x='Population', y='Allele Frequency', hue="LOF", kind='swarm', legend=False,
                 aspect=1.5, height=5)
plt.ylim(0, 0.0005)
plt.legend(loc='upper right')
plt.show()

""" COMBINED GRAPH """

all_variants = create_cat_df([
    (['../analysis/inesss_variants.tsv'], 'INESSS'),
    (['../analysis/clin_variants.tsv'], 'CLINVAR'),
    (['../analysis/lof_variants.tsv'], 'LOF'),
    (['../analysis/cpic_variants.tsv'], 'CPIC')
])

all_variants = create_basic_df(all_variants, ('CATEGORY', 'CATEGORY'))

ax = sns.catplot(data=all_variants, x='Population', y='Allele Frequency', hue="CATEGORY", kind='swarm', legend=False,
                 aspect=1.5, height=5)
plt.legend(loc='upper right')
plt.ylim(0, 0.03)
plt.show()

""" INESSS VS CLINVAR VS LOF BY POPULATION """

all_variants = create_cat_df([
    (['../analysis/inesss_variants.tsv'], 'INESSS'),
    (['../analysis/clin_variants.tsv'], 'CLINVAR'),
    (['../analysis/lof_variants.tsv'], 'LOF'),
    (['../analysis/cpic_variants.tsv'], 'CPIC')
])

all_variants_pop = create_sum_df(all_variants, True)

ax = sns.catplot(data=all_variants_pop, x='Population', y='Allele Frequency', hue='CATEGORY', kind='bar', legend=False,
                 aspect=2)
plt.legend(loc='upper left')
plt.ylim(0, 0.03)
plt.show()

""" ENTIRE POPULATION GRAPH """

all_variants_total = create_sum_df(all_variants, False)

ax = sns.catplot(data=all_variants_total, x='CATEGORY', y='Allele Frequency', kind='bar', legend=False, aspect=1.5)
plt.legend(loc='upper right')
plt.ylim(0, 0.03)
plt.show()

""" INESSS VS NON-INESSS """

inesss_vs_non = create_cat_df([
    (['../analysis/inesss_variants.tsv'], 'INESSS'),
    (['../analysis/clin_variants.tsv', '../analysis/lof_variants.tsv', '../analysis/cpic_variants.tsv'], 'NON-INESSS'),
])

inesss_pop = create_sum_df(inesss_vs_non, True)

ax = sns.catplot(data=inesss_pop, x='Population', y='Allele Frequency', hue='CATEGORY', kind='swarm', legend=False,
                 aspect=2)
plt.legend(loc='upper left')
plt.ylim(0, 0.04)
plt.show()

""" Y-axis is proportion of frequency """

# Get total population frequencies
population_frequencies = {population: 0 for population in populations}

for key, dict in inesss_pop.iterrows():
    population_frequencies[dict['Population']] += dict['Allele Frequency']

# Create dict for frequencies of both
freq_inesss = {}
freq_noninesss = {}

for key, dict in inesss_pop.iterrows():
    population = dict['Population']

    if dict['CATEGORY'] == 'INESSS':
        freq_inesss[population] = (dict['Allele Frequency'] / population_frequencies[population]) * 100
    else:
        freq_noninesss[population] = (dict['Allele Frequency'] / population_frequencies[population]) * 100

freq_dict = {
    'Population': populations,
    'Proportion INESSS': list(freq_inesss.values()),
    'Proportion NONINESSS': list(freq_noninesss.values()),
}
freq_df = pd.DataFrame.from_dict(freq_dict)
width = 0.5

p1 = plt.bar(populations, list(freq_inesss.values()), width, color="#ec6b2d")
p2 = plt.bar(populations, list(freq_noninesss.values()), width,
             bottom=list(freq_inesss.values()), color="#00aabb")

plt.ylabel('Proportion')
plt.xlabel('Population')
plt.legend((p1[0], p2[0]), ('INESSS', 'NON INESSS'))

plt.show()

# sns.barplot(x='Population',y = 'Total Frequency', data=freq_df, color = "darkblue")
# sns.barplot(x='Population',y = 'Proportion INESSS', data=freq_df, color = "skyblue")
# sns.barplot(x='Population',y = 'Proportion NONINESSS', data=freq_df, color = "pink")
#
# plt.xlabel('Population')
# plt.ylabel('Proportion')
#
# for index, row in freq_df.iterrows():
#     txt = str(round(row['Proportion INESSS'], 2)) + " %"
#     plt.text(index, row['Proportion INESSS'], txt, ha = "center")
# plt.show()
