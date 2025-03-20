import pandas as pd
import os

dataframe = pd.read_csv('/home/kolos100/PHD/cancer_clustering_app/source_data/approved_symbols_sep.txt',
                        dtype=str)

# get all the unique symbols
symbols = list(dataframe['approved_symbol'].unique().tolist() +
              dataframe['previous_symbol_1'].dropna().unique().tolist() +
              dataframe['previous_symbol_2'].dropna().unique().tolist() +
              dataframe['previous_symbol_3'].dropna().unique().tolist() +
              dataframe['previous_symbol_ 4'].dropna().unique().tolist() +
              dataframe['previous_symbol_ 5'].dropna().unique().tolist() +
              dataframe['previous_symbol_6'].dropna().unique().tolist() +
              dataframe['previous_symbol_7'].dropna().unique().tolist() +
              dataframe['previous_symbol_ 8'].dropna().unique().tolist() +
              dataframe['previous_symbol_9'].dropna().unique().tolist() +
              dataframe['previous_symbol_10'].dropna().unique().tolist() +
              dataframe['previous_symbol_11'].dropna().unique().tolist() +
              dataframe['previous_symbol_12'].dropna().unique().tolist() +
              dataframe['alias_symbol_1'].dropna().unique().tolist() +
              dataframe['alias_symbol_2'].dropna().unique().tolist() +
              dataframe['alias_symbol_3'].dropna().unique().tolist() +
              dataframe['alias_symbol_4'].dropna().unique().tolist() +
              dataframe['alias_symbol_5'].dropna().unique().tolist() +
              dataframe['alias_symbol_6'].dropna().unique().tolist())

print("Gene dict load , Done!")

all_gene_names = pd.DataFrame()
all_gene_names['Gene_Names'] = symbols
all_gene_names.set_index('Gene_Names', inplace=True)

# create a dictionary with gene names as keys and associated names as values
gene_dict = {}
for gene_name in symbols:
    associations = dataframe[(dataframe['approved_symbol'] == gene_name) |
                             (dataframe['previous_symbol_1'] == gene_name) |
                             (dataframe['previous_symbol_2'] == gene_name) |
                             (dataframe['previous_symbol_3'] == gene_name) |
                             (dataframe['previous_symbol_ 4'] == gene_name) |
                             (dataframe['previous_symbol_ 5'] == gene_name) |
                             (dataframe['previous_symbol_6'] == gene_name) |
                             (dataframe['previous_symbol_7'] == gene_name) |
                             (dataframe['previous_symbol_ 8'] == gene_name) |
                             (dataframe['previous_symbol_9'] == gene_name) |
                             (dataframe['previous_symbol_10'] == gene_name) |
                             (dataframe['previous_symbol_11'] == gene_name) |
                             (dataframe['previous_symbol_12'] == gene_name) |
                             (dataframe['alias_symbol_1'] == gene_name) |
                             (dataframe['alias_symbol_2'] == gene_name) |
                             (dataframe['alias_symbol_3'] == gene_name) |
                             (dataframe['alias_symbol_4'] == gene_name) |
                             (dataframe['alias_symbol_5'] == gene_name) |
                             (dataframe['alias_symbol_6'] == gene_name)]['approved_symbol'].tolist()
    gene_dict[gene_name] = associations

# use map() to populate the new column with associations
all_gene_names['associations'] = all_gene_names.index.map(gene_dict)
all_gene_names['count_associations'] = all_gene_names['associations'].apply(lambda x: len(x))

print("Stage 2 Done")

george_df = pd.read_csv('/home/kolos100/PHD/cancer_clustering_app/source_data/data_SCLC NCI-DTP_exp.txt', delimiter="\t", index_col=3)
george_df = george_df.drop(columns=['TranscriptClusterID',"UNIT_ID", "GeneName", "GeneAccession", "EntrezID",
                                    "Chromosome", "Cytoband", "Start", "Stop", "Strand",
                                    "CrossHybridization", "ProbesetType"], axis=1)
george_df = george_df.assign(row_sum=george_df.sum(axis=1))
george_df = george_df.sort_values(by=['row_sum'], ascending=False)
george_df = george_df.drop_duplicates(subset=['GeneSymbol'], keep='first')
#george_df.set_index("Hugo_Symbol", inplace=True)
george_df = george_df.drop(columns=['row_sum'])
george_df['Avg_Expr'] = george_df.mean(axis=1)

select_list = pd.DataFrame()
select_list["Names"] = george_df['Hugo_Symbol']
select_list['Expr_Mean'] = george_df['Avg_Expr']
select_list_name = 'select_list_name.csv'
select_list.to_csv(select_list_name, index=False, header=False, encoding="utf")

with open(select_list_name, "r") as f:
    gene_ids = [row.split(',')[0] for row in f.read().splitlines()]
with open(select_list_name, "r") as f:
    gene_value = [row.split(',')[1] for row in f.read().splitlines()]

# get the list of missing gene names
missing_gene_names = list(set(gene_ids) - set(all_gene_names.index.tolist()))

all_gene_names = all_gene_names[~all_gene_names.index.duplicated(keep='first')]
all_gene_names = all_gene_names.reindex(all_gene_names.index.union(missing_gene_names), fill_value=0)

# add missing gene names to the index and fill NaN values with 0 in 'associations' column
all_gene_names = all_gene_names.reindex(all_gene_names.index.union(missing_gene_names),
                                        fill_value={'associations': 0, 'count_associations': 0})

# select only rows from the dataframe where the 'Gene' column matches an ID in the list
selected_genes = all_gene_names.loc[gene_ids]

# add the 'gene_value' column to the selected_genes dataframe
selected_genes['gene_value'] = gene_value

# define function to remove brackets from list elements
def remove_brackets(x):
    if isinstance(x, list):
        if len(x) == 1:
            return x[0].strip('[]')
        else:
            return ','.join(x)
    else:
        return x

# apply function to col2
selected_genes['associations'] = selected_genes['associations'].apply(remove_brackets)

selected_genes_fals = selected_genes[selected_genes.index!= selected_genes['associations']]
selected_genes_corr = selected_genes[selected_genes.index == selected_genes['associations']]
selected_genes_corr = selected_genes_corr.drop(columns=['associations', 'count_associations'])

# save the selected genes dataframe to a CSV file
selected_genes_fals.to_csv("/home/kolos100/PHD/cancer_clustering_app/results_data/mapped_names_fals.csv", header=True, index=True)
selected_genes_corr.to_csv("/home/kolos100/PHD/cancer_clustering_app/results_data/mapped_names_corr.csv", header=True, index=True)

os.remove(select_list_name)

print("Done")