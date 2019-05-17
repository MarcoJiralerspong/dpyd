import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns



""" SETUP DATA """

populations = ["eas", "afr", "amr", "asj", "sas", "nfe", "fin"]
filtered_variants = pd.read_csv('../data/filtered_variants.tsv', sep='\t')
filtered_variants_inesss = filtered_variants.loc[filtered_variants['INESSS'] == True] # Only inesss variants

freq_lof = {}
freq_inesss = {}

# Get individual population allele frequencies
for population in populations:
    ac_pop = 'AC_' + population
    an_pop = 'AN_' + population

    # To deal with NaN, get indices where they are NaN and only sum when not NaN (issue with sas pop)
    null_index = filtered_variants[an_pop].map(np.isnan)
    notnull_index = [not i for i in null_index]

    # Compute frequency over all variants for only inesss as well as all
    freq_inesss[population] = np.sum(filtered_variants_inesss[ac_pop].map(int) / filtered_variants_inesss[an_pop].map(int))
    freq_lof[population] = np.sum(filtered_variants[ac_pop][notnull_index].map(int) / filtered_variants[an_pop][notnull_index].map(int))

# Transform to dataframe
freq_dict = {
    'populations': populations,
    'freq_lof': list(freq_lof.values()),
    'freq_inesss': list(freq_inesss.values())
}
freq_df = pd.DataFrame.from_dict(freq_dict)

# Get relative proportion of inesss over all variants
freq_df['inesss_lof_freq'] = freq_df['freq_inesss']/freq_df['freq_lof']


""" PLOT DATA """

sns.set_style("darkgrid")
sns.barplot(x='populations',y = 'freq_lof', data=freq_df, color = "skyblue")
sns.barplot(x='populations',y = 'freq_inesss', data=freq_df, color = "pink")

plt.xlabel('Population')
plt.ylabel('Frequency')
# plt.title('Allele frequence of LOF variants of DPYD in different populations and proportion covered by INESSS variants')

for index, row in freq_df.iterrows():
    txt = str(round(row['inesss_lof_freq'] * 100, 2)) + " %"
    plt.text(index, row['freq_lof'], txt, ha = "center")
plt.show()
