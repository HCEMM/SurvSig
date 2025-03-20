#from unittest.mock import inplace

import pandas as pd
import streamlit as st
#import openpyxl
#import pyxlsb
#from memory_profiler import profile


@st.cache_resource(show_spinner='Gene Dictionary is Loading')
def gene_autocorrect():
    # Load the gene name dictionary into a DataFrame
    dataframe = pd.read_csv('/home/hcemm-user/SurvSig/source_data/approved_symbols_sep.txt', dtype=str)

    # create an empty dictionary to store the names
    name_dictionary = {}

    # loop over the dataframe and add the names to the dictionary
    for i in range(len(dataframe)):
        target_name = dataframe.loc[i, 'approved_symbol']
        for col_name in dataframe.columns:
            if col_name.startswith('approved_symbol') or col_name.startswith('previous_symbol_') or col_name.startswith(
                    'alias_symbol_'):
                old_name = dataframe.loc[i, col_name]
                if pd.notnull(old_name):
                    if old_name not in name_dictionary:
                        name_dictionary[old_name] = []
                    name_dictionary[old_name].append(target_name)

    return name_dictionary, dataframe


@st.cache_data(show_spinner="NCI-CELL 60 Line loading ...")
def nci_60_xai_data_load():
    nci_sclc_cell = pd.read_csv('/home/kolos/Cancer_Cluster_Application/source_data/data_NCI-60_xai.txt',
                                delimiter="\t", index_col=1)
    nci_sclc_cell = nci_sclc_cell.drop(["Probe id", "Entrez gene id", "Chromosome", "Start", "End", "Cytoband",
                                        "RefSeq(protein)", "RefSeq(mRNA)"], axis=1)
    nci_sclc_cell = nci_sclc_cell.assign(row_sum=nci_sclc_cell.sum(axis=1))
    nci_sclc_cell = nci_sclc_cell.sort_values(by=['row_sum'], ascending=False)
    nci_sclc_cell.index.rename("GeneSymbol", inplace=True)
    nci_sclc_cell = nci_sclc_cell.reset_index()
    nci_sclc_cell = nci_sclc_cell.drop_duplicates(subset="GeneSymbol", keep='first')
    nci_sclc_cell.set_index("GeneSymbol", inplace=True)
    nci_sclc_cell = nci_sclc_cell.drop(columns=['row_sum'])

    return nci_sclc_cell


@st.cache_data(show_spinner="GDSC XSQ Line loading ...")
def gdsc_xsq_data_load():
    gdsc_xsq = pd.read_csv('/home/kolos/Cancer_Cluster_Application/source_data/data_GDSC-MGH-Sanger_xsq.txt',
                                delimiter="\t", index_col=0)
    st.write(gdsc_xsq.head())
    gdsc_xsq = gdsc_xsq.drop(["ensembl_gene_id", "hgnc_id", "entrez_id", "refseq_id", "uniprot_id"], axis=1)
    gdsc_xsq = gdsc_xsq.assign(row_sum=gdsc_xsq.sum(axis=1))
    gdsc_xsq = gdsc_xsq.sort_values(by=['row_sum'], ascending=False)
    gdsc_xsq.index.rename("GeneSymbol", inplace=True)
    gdsc_xsq = gdsc_xsq.reset_index()
    gdsc_xsq = gdsc_xsq.drop_duplicates(subset="GeneSymbol", keep='first')
    gdsc_xsq.set_index("GeneSymbol", inplace=True)
    gdsc_xsq = gdsc_xsq.drop(columns=['row_sum'])

    return gdsc_xsq


@st.cache_data(show_spinner="GDSC EXP Line loading ...")
def gdsc_exp_data_load():
    gdsc_xsq = pd.read_csv('/home/kolos/Cancer_Cluster_Application/source_data/data_GDSC-MGH-Sanger_exp.txt',
                                delimiter="\t", index_col=0)
    gdsc_xsq = gdsc_xsq.drop(["id", "id_type"], axis=1)
    gdsc_xsq = gdsc_xsq.assign(row_sum=gdsc_xsq.sum(axis=1))
    gdsc_xsq = gdsc_xsq.sort_values(by=['row_sum'], ascending=False)
    gdsc_xsq.index.rename("GeneSymbol", inplace=True)
    gdsc_xsq = gdsc_xsq.reset_index()
    gdsc_xsq = gdsc_xsq.drop_duplicates(subset="GeneSymbol", keep='first')
    gdsc_xsq.set_index("GeneSymbol", inplace=True)
    gdsc_xsq = gdsc_xsq.drop(columns=['row_sum'])

    return gdsc_xsq


@st.cache_data(show_spinner="CCLE-CELL EXP Line loading ...")
def ccle_cell_exp_data_load():
    ccle_exp = pd.read_csv('/home/kolos/Cancer_Cluster_Application/source_data/data_CCLE-Broad-MIT_exp.txt',
                           delimiter="\t", index_col=0)
    ccle_exp = ccle_exp.assign(row_sum=ccle_exp.sum(axis=1))
    ccle_exp = ccle_exp.sort_values(by=['row_sum'], ascending=False)
    ccle_exp.index.rename("GeneSymbol", inplace=True)
    ccle_exp = ccle_exp.reset_index()
    ccle_exp = ccle_exp.drop_duplicates(subset="GeneSymbol", keep='first')
    ccle_exp.set_index("GeneSymbol", inplace=True)
    ccle_exp = ccle_exp.drop(columns=['row_sum'])

    st.write(ccle_exp)

    return ccle_exp


@st.cache_data(show_spinner="CCLE-CELL XSQ Line loading ...")
def ccle_cell_xsq_data_load():
    ccle_xsq = pd.read_csv('/home/kolos/Cancer_Cluster_Application/source_data/data_CCLE-Broad-MIT_xsq.txt',
                           delimiter="\t", index_col=2)
    ccle_xsq = ccle_xsq.drop(["Name", "ID"], axis=1)
    ccle_xsq = ccle_xsq.assign(row_sum=ccle_xsq.sum(axis=1))
    ccle_xsq = ccle_xsq.sort_values(by=['row_sum'], ascending=False)
    ccle_xsq.index.rename("GeneSymbol", inplace=True)
    ccle_xsq = ccle_xsq.reset_index()
    ccle_xsq = ccle_xsq.drop_duplicates(subset="GeneSymbol", keep='first')
    ccle_xsq.set_index("GeneSymbol", inplace=True)
    ccle_xsq = ccle_xsq.drop(columns=['row_sum'])

    return ccle_xsq


@st.cache_data(show_spinner="NCI-CELL 60 Line loading ...")
def nci_60_xsq_data_load():
    nci_sclc_cell = pd.read_csv('/home/kolos/Cancer_Cluster_Application/source_data/data_NCI-60_xsq.txt',
                                delimiter="\t", index_col=0)
    nci_sclc_cell = nci_sclc_cell.drop(["chr", "txStart", "txEnd", "entrez.gene.id"], axis=1)
    nci_sclc_cell = nci_sclc_cell.assign(row_sum=nci_sclc_cell.sum(axis=1))
    nci_sclc_cell = nci_sclc_cell.sort_values(by=['row_sum'], ascending=False)
    nci_sclc_cell.index.rename("GeneSymbol", inplace=True)
    nci_sclc_cell = nci_sclc_cell.reset_index()
    nci_sclc_cell = nci_sclc_cell.drop_duplicates(subset="GeneSymbol", keep='first')
    nci_sclc_cell.set_index("GeneSymbol", inplace=True)
    nci_sclc_cell = nci_sclc_cell.drop(columns=['row_sum'])

    st.write(nci_sclc_cell)

    return nci_sclc_cell


@st.cache_data(show_spinner="NCI-CELL Line loading ...")
def nci_sclc_data_load():
    nci_sclc_cell = pd.read_csv('/home/kolos100/PHD/Cancer_Cluster_Application/source_data/data_SCLC NCI-DTP_exp.txt',
                                delimiter="\t", index_col=0)
    nci_sclc_cell = nci_sclc_cell.drop(["TranscriptClusterID", "UNIT_ID", "GeneName", "GeneAccession", "EntrezID",
                                        "Chromosome", "Cytoband", "Start", "Stop", "Strand",
                                        "CrossHybridization", "ProbesetType"], axis=1)
    nci_sclc_cell = nci_sclc_cell.assign(row_sum=nci_sclc_cell.sum(axis=1))
    nci_sclc_cell = nci_sclc_cell.sort_values(by=['row_sum'], ascending=False)
    nci_sclc_cell = nci_sclc_cell.reset_index()
    nci_sclc_cell = nci_sclc_cell.drop_duplicates(subset="GeneSymbol", keep='first')
    # nci_sclc_cell.set_index("Hugo_Symbol", inplace=True)
    nci_sclc_cell = nci_sclc_cell.drop(columns=['row_sum'])

    return nci_sclc_cell


@st.cache_data(show_spinner="SCLC-GDSC Cell Line loading ...")
def gdsc_sclc_data_load():
    nci_sclc_cell = pd.read_csv('/home/kolos100/PHD/Cancer_Cluster_Application/source_data/data_SCLC '
                                'GDSC-MGH-Sanger_exp.txt',
                                delimiter="\t", index_col=0, header=0)
    nci_sclc_cell = nci_sclc_cell.drop(['id', 'id_type'], axis=1)
    nci_sclc_cell.index.rename('GeneSymbol', inplace=True)
    nci_sclc_cell = nci_sclc_cell.assign(row_sum=nci_sclc_cell.sum(axis=1))
    nci_sclc_cell = nci_sclc_cell.sort_values(by=['row_sum'], ascending=False)
    nci_sclc_cell = nci_sclc_cell.reset_index()
    nci_sclc_cell = nci_sclc_cell.drop_duplicates(subset="GeneSymbol", keep='first')
    nci_sclc_cell.set_index("GeneSymbol", inplace=True)
    nci_sclc_cell = nci_sclc_cell.drop(columns=['row_sum'])
    nci_sclc_cell.dropna(how='all', axis=1, inplace=True)
    st.write(nci_sclc_cell.head(50))

    return nci_sclc_cell


@st.cache_data(show_spinner="SCLC-UTSW Cell Line loading ...")
def utsw_sclc_data_load():
    utsw_sclc_cell = pd.read_csv('/home/kolos100/PHD/Cancer_Cluster_Application/source_data/data_SCLC UTSW_xsq.txt',
                                 delimiter="\t", index_col=3, header=0)
    utsw_sclc_cell = utsw_sclc_cell.drop(['ID', 'Gene Symbol', 'Gene ID', 'Entrez ID', 'Gene Type',
                                          'Locus'], axis=1)
    utsw_sclc_cell = utsw_sclc_cell.assign(row_sum=utsw_sclc_cell.sum(axis=1))
    utsw_sclc_cell.sort_values(by=['row_sum'], ascending=False)
    utsw_sclc_cell.reset_index(inplace=True)
    utsw_sclc_cell.dropna(how='all', axis=1, inplace=True)
    utsw_sclc_cell.dropna(how='any', axis=0, inplace=True)
    utsw_sclc_cell.rename(columns={"HUGO Symbol": "GeneSymbol"}, inplace=True)
    utsw_sclc_cell.drop_duplicates(subset="GeneSymbol", keep='first', inplace=True)
    utsw_sclc_cell.set_index("GeneSymbol", inplace=True)
    utsw_sclc_cell.drop(columns=['row_sum'])
    utsw_sclc_cell.dropna(how='all', axis=1, inplace=True)
    st.write(utsw_sclc_cell.head(50))

    return utsw_sclc_cell


@st.cache_data(show_spinner="SCLC-CCLE Cell Line loading ...")
def ccle_sclc_data_load():
    ccle_sclc_cell = pd.read_csv(
        '/home/kolos100/PHD/Cancer_Cluster_Application/source_data/data_SCLC CCLE-Broad-MIT_xsq.txt',
        delimiter="\t", index_col=0, header=0)
    ccle_sclc_cell = ccle_sclc_cell.drop(['Description', 'Name'], axis=1)
    ccle_sclc_cell = ccle_sclc_cell.assign(row_sum=ccle_sclc_cell.sum(axis=1))
    ccle_sclc_cell.sort_values(by=['row_sum'], ascending=False)
    ccle_sclc_cell.reset_index(inplace=True)
    ccle_sclc_cell.dropna(how='all', axis=1, inplace=True)
    ccle_sclc_cell.dropna(how='any', axis=0, inplace=True)
    ccle_sclc_cell.rename(columns={"ID": "GeneSymbol"}, inplace=True)
    ccle_sclc_cell.drop_duplicates(subset="GeneSymbol", keep='first', inplace=True)
    ccle_sclc_cell.set_index("GeneSymbol", inplace=True)
    ccle_sclc_cell.drop(columns=['row_sum'])
    ccle_sclc_cell.dropna(how='all', axis=1, inplace=True)
    st.write(ccle_sclc_cell.head(50))

    return ccle_sclc_cell


@st.cache_data(show_spinner="Anish loading ...")
def anish_sclc_data_load():
    anish_sclc_data = pd.read_excel('/home/kolos100/PHD/Cancer_Cluster_Application/source_data'
                                    '/Thomas_RNASeq_log2_TMM_FPKM_only_for_NE_project_patient_tumor_20201004.xlsx',
                                    index_col=0)

    anish_sclc_data = anish_sclc_data.assign(row_sum=anish_sclc_data.sum(axis=1))
    anish_sclc_data = anish_sclc_data.sort_values(by=['row_sum'], ascending=False)
    anish_sclc_data = anish_sclc_data.reset_index()
    anish_sclc_data = anish_sclc_data.drop_duplicates(subset="Sample_id", keep='first')
    anish_sclc_data.set_index("Sample_id", inplace=True)
    anish_sclc_data = anish_sclc_data.drop(columns=['row_sum'])

    return anish_sclc_data

@st.cache_data(show_spinner="LNEC loading ...")
def lcnec_sclc_data_load():
    lcnec_sclc_data = pd.read_excel('/home/kolos/website_test_ccapp_v2_lite/source_data/LCNEC_George_expr.xlsx',
                                    index_col=0)

    lcnec_sclc_data = lcnec_sclc_data.assign(row_sum=lcnec_sclc_data.sum(axis=1))
    lcnec_sclc_data = lcnec_sclc_data.sort_values(by=['row_sum'], ascending=False)
    lcnec_sclc_data = lcnec_sclc_data.reset_index()
    lcnec_sclc_data = lcnec_sclc_data.drop_duplicates(subset="Gene", keep='first')
    lcnec_sclc_data.rename(columns={"Gene": "GeneSymbol"}, inplace=True)
    lcnec_sclc_data.set_index("GeneSymbol", inplace=True)
    lcnec_sclc_data = lcnec_sclc_data.drop(columns=['row_sum'])

    st.write(lcnec_sclc_data.head(50))

    return lcnec_sclc_data

@st.cache_data(show_spinner="Carcinoid loading ...")
def carcinoid_data_load():
    carcinoid_data = pd.read_excel('/home/kolos/website_test_ccapp_v2_lite/source_data/41467_2014_BFncomms4518_MOESM1020_ESM.xlsx',
                                    index_col=0)
    carcinoid_data.drop(columns=['ID'], inplace=True)
    carcinoid_data = carcinoid_data.assign(row_sum=carcinoid_data.sum(axis=1))
    carcinoid_data = carcinoid_data.sort_values(by=['row_sum'], ascending=False)
    carcinoid_data = carcinoid_data.reset_index()
    carcinoid_data = carcinoid_data.drop_duplicates(subset="Gene", keep='first')
    carcinoid_data.rename(columns={"Gene": "GeneSymbol"}, inplace=True)
    carcinoid_data.set_index("GeneSymbol", inplace=True)
    carcinoid_data = carcinoid_data.drop(columns=['row_sum'])

    st.write(carcinoid_data.head(50))

    return carcinoid_data


@st.cache_data(show_spinner="LCNEC cell loading ...")
def cell_line_szandi():
    lcenc_cell_data = pd.read_csv('/home/hcemm-user/SurvSig/source_data/lcnec_cell.tsv',  sep='\t')
    lcenc_cell_data = lcenc_cell_data.transpose()

    lcenc_cell_data = lcenc_cell_data.assign(row_sum=lcenc_cell_data.sum(axis=1))
    #lcenc_cell_data = lcenc_cell_data.sort_values(by=['row_sum'], ascending=False)
    #lcenc_cell_data = lcenc_cell_data.reset_index()
    #lcenc_cell_data = lcenc_cell_data.drop_duplicates(subset="index", keep='first')
    #lcenc_cell_data.rename(columns={"index": "GeneSymbol"}, inplace=True)
    #lcenc_cell_data.set_index("GeneSymbol", inplace=True)
    #lcenc_cell_data = lcenc_cell_data.drop(columns=['row_sum'])

    st.write(lcenc_cell_data.head(50))

    return lcenc_cell_data



@st.cache_data(show_spinner="Alcala loading ...")
def alcala_data_load():
    alcala_data = pd.read_csv('/home/kolos/website_test_ccapp_v2_lite/Carcinoid_Alcala_pan_sample/Carcinoid_pan_sample.tsv',
                                    index_col=0, sep="\t")
    alcala_data.index.rename("Gene", inplace=True)
    alcala_data = alcala_data.assign(row_sum=alcala_data.sum(axis=1))
    alcala_data = alcala_data.sort_values(by=['row_sum'], ascending=False)
    alcala_data = alcala_data.reset_index()
    alcala_data = alcala_data.drop_duplicates(subset="Gene", keep='first')
    alcala_data.rename(columns={"Gene": "GeneSymbol"}, inplace=True)
    alcala_data.set_index("GeneSymbol", inplace=True)
    alcala_data = alcala_data.drop(columns=['row_sum'])

    st.write(alcala_data.head(50))

    return alcala_data


def convert_string_to_float(value):
    if isinstance(value, str) and ',' in value:
        # Replace comma with dot and convert to float
        return float(value.replace(',', '.'))
    else:
        # Return the original value if not a string or no comma present
        return value

@st.cache_data(show_spinner="Liu loading ...")
def liu_data_load():

    liu_data = pd.read_csv('/home/kolos/website_test_ccapp_v2_lite/source_data/liu_expression_unmapped.csv',
                                    index_col=0)
    liu_data.index.rename("Gene", inplace=True)
    liu_data.fillna(0, inplace=True)
    columns_to_drop = [col for col in liu_data.columns if col.startswith('N')]
    liu_data = liu_data.drop(columns=columns_to_drop, axis=1)
    #liu_data.rename(columns=lambda x: x.replace('T', ''), inplace=True)
    liu_data = liu_data.applymap(convert_string_to_float)
    liu_data.astype('float')

    liu_data = liu_data.assign(row_sum=liu_data.sum(axis=1))
    liu_data = liu_data.sort_values(by=['row_sum'], ascending=False)
    liu_data = liu_data.reset_index()
    liu_data = liu_data.drop_duplicates(subset="Gene", keep='first')
    liu_data.rename(columns={"Gene": "GeneSymbol"}, inplace=True)
    liu_data.set_index("GeneSymbol", inplace=True)
    liu_data = liu_data.drop(columns=['row_sum'])

    st.write(liu_data.head(50))

    return liu_data


@st.cache_data(show_spinner="Liu loading ...")
def liu_data_load_proteo():

    liu_data = pd.read_excel('/home/kolos/website_test_ccapp_v2_lite/source_data/Table S1_LIU.xlsb',
                             sheet_name=5, index_col=0, header=0)
    liu_data.index.rename("Gene", inplace=True)
    liu_data = liu_data.fillna(liu_data.min())
    columns_to_drop = [col for col in liu_data.columns if col.startswith('N')]
    liu_data = liu_data.drop(columns=columns_to_drop, axis=1)
    #liu_data.rename(columns=lambda x: x.replace('T', ''), inplace=True)
    liu_data = liu_data.applymap(convert_string_to_float)
    liu_data.astype('float')

    liu_data = liu_data.assign(row_sum=liu_data.sum(axis=1))
    liu_data = liu_data.sort_values(by=['row_sum'], ascending=False)
    liu_data = liu_data.reset_index()
    liu_data = liu_data.drop_duplicates(subset="Gene", keep='first')
    liu_data.rename(columns={"Gene": "GeneSymbol"}, inplace=True)
    liu_data.set_index("GeneSymbol", inplace=True)
    liu_data = liu_data.drop(columns=['row_sum'])

    st.write(liu_data.head(50))

    return liu_data


@st.cache_data(show_spinner="Rousseaux loading ...")
def rousseaux_data_load():
    rousseaux_data = pd.read_csv('/home/kolos/website_test_ccapp_v2_lite/source_data/GSE30219.rma_norm.hgu133plus2_jetset.tsv',
                                    index_col=0, sep="\t")

    st.write(rousseaux_data.shape)
    rousseaux_data.index.rename("GeneSymbol", inplace=True)
    rousseaux_data = rousseaux_data.assign(row_sum=rousseaux_data.sum(axis=1))
    rousseaux_data = rousseaux_data.sort_values(by=['row_sum'], ascending=False)
    rousseaux_data = rousseaux_data.reset_index()
    rousseaux_data = rousseaux_data.drop_duplicates(subset="GeneSymbol", keep='first')
    rousseaux_data.set_index("GeneSymbol", inplace=True)
    rousseaux_data = rousseaux_data.drop(columns=['row_sum'])

    st.write(rousseaux_data.head(50))
    st.write(rousseaux_data.shape)

    return rousseaux_data



@st.cache_data(show_spinner="TCGA loading ...")
def tcga_data_load():
    tcga_df = pd.read_csv("/home/kolos100/PHD/Cancer_Cluster_Application/TCGA_prep"
                          "/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.renamed.tsv", delimiter="\t")
    unnamed_col = tcga_df.iloc[:, 0]
    unnamed_col.name = 'Genes'
    tcga_df.rename(columns={tcga_df.columns[0]: unnamed_col.name}, inplace=True)
    tcga_df[['Gene', 'Gene_ID']] = tcga_df['Genes'].str.split('|', expand=True)
    tcga_df = tcga_df.drop('Genes', axis=1)
    tcga_df = tcga_df.drop('Gene_ID', axis=1)
    tcga_df.set_index("Gene", inplace=True)

    tcga_df = tcga_df.assign(row_sum=tcga_df.sum(axis=1))
    tcga_df = tcga_df.sort_values(by=['row_sum'], ascending=False)
    tcga_df = tcga_df.reset_index()
    tcga_df = tcga_df.drop_duplicates(subset="Gene", keep='first')
    tcga_df.set_index("Gene", inplace=True)
    tcga_df = tcga_df.drop(columns=['row_sum'])

    return tcga_df


def app(name_dict, df_app, gene_names, file_name):
    # Set the app title
    st.title('Gene Name Suggestions')

    corrected_names = []
    original_names = []
    not_found_names = set()

    # Get the gene names, select first column, drop the duplicates, handle whitespaces, and convert into list
    gene_names.replace('^\s+', '', regex=True, inplace=True)  # front
    gene_names.replace('\s+$', '', regex=True, inplace=True)  # end whitespace remove

    # Create a list to store the corrected gene names and a set to store the not found gene names
    corrected_names = []
    not_found_names = set()

    # Loop over the gene names and check for suggestions
    for gene_name in gene_names.index:
        if gene_name not in name_dict:
            # If gene_name is not in the name_dict, add it to the not_found_names set
            not_found_names.add(gene_name)
            continue
        elif gene_name in name_dict and gene_name in df_app['approved_symbol'].values:
            # If gene_name is in approved_symbol but also has synonyms or alias names,
            # use the approved_symbol name automatically and don't suggest any correction
            gene_name = gene_name
        elif gene_name in name_dict:
            target_names = name_dict[gene_name]
            if len(target_names) == 1:
                suggestion = f'{gene_name} (consider using {target_names[0]})'
                options = [target_names[0], 'Remove gene']
                option = st.selectbox(suggestion, options, key=gene_name)
                if option == 'Remove gene':
                    continue
                original_names.append(gene_name)  # Add this line
                gene_name = option
                corrected_names.append(gene_name)
            elif len(target_names) == 2:
                suggestion = f'{gene_name} (consider using {", ".join(target_names)})'
                options = target_names + ['Remove gene']
                option = st.selectbox(suggestion, options)
                if option == 'Remove gene':
                    continue
                original_names.append(gene_name)  # Add this line
                gene_name = option
                corrected_names.append(gene_name)
            elif len(target_names) > 2:
                suggestion = f'{gene_name} (consider using {", ".join(target_names)})'
                options = target_names + ['Remove gene']
                option = st.selectbox(suggestion, options)
                if option == 'Remove gene':
                    continue
                original_names.append(gene_name)  # Add this line
                gene_name = option
                corrected_names.append(gene_name)

        else:
            st.warning(f"No suggestion available for gene name: {gene_name}")
            not_found_names.add(gene_name)
    # Remove the not found gene names from the input DataFrame
    gene_names = gene_names[~gene_names.index.isin(not_found_names)]

    # Display the corrected input DataFrame
    st.write('## Corrected Input Data')
    st.write("Number of Corrected Names of Genes:    " + str(len(corrected_names)))
    # gene_names.set_index('Gene', inplace=True)

    # Display the not found gene names
    if len(not_found_names) != 0:
        st.warning("Number of not founded genes:    " + str(len(not_found_names)))
        st.warning(f"The following names/terms are not in the gene dictionary: {', '.join(sorted(not_found_names))}")
    else:
        pass
    # Create a new DataFrame with the original and corrected gene names
    # Create a DataFrame with corrected names
    corrected_df = pd.DataFrame({'original_name': original_names, 'corrected_name': corrected_names})
    not_found_names_df = pd.DataFrame({'original_name': list(not_found_names), 'corrected_name': 'not found'})
    summary_df = pd.concat([corrected_df, not_found_names_df], ignore_index=True, sort=False)
    # summary_df.to_csv("/home/kolos100/PHD/Cancer_Cluster_Application/source_data/summary_df.csv", header=True,
    #                  index=False)

    # create a dictionary from the second dataframe
    mapping_dict = {k: v for k, v in zip(corrected_df['original_name'], corrected_df['corrected_name'])}
    gene_names.reset_index(inplace=True)
    # map the index of the first dataframe using the dictionary
    # replace the values in the 'Hugo_Symbol' column using the mapping dictionary
    gene_names['GeneSymbol'].replace(mapping_dict, inplace=True)
    gene_names.set_index('GeneSymbol', inplace=True)
    # count the number of duplicated index values
    num_duplicates = gene_names.copy().index.duplicated().sum()
    st.write(num_duplicates)
    # find the duplicate index values
    duplicates = gene_names.index[gene_names.index.duplicated(keep=False)]

    # select the rows corresponding to the duplicate index values
    duplicate_rows = gene_names.loc[duplicates]
    st.write(duplicates)

    # Display the original vs. changed names DataFrame
    st.write(summary_df.transpose())

    gene_names = gene_names.assign(row_sum=gene_names.sum(axis=1))
    gene_names = gene_names.sort_values(by=['row_sum'], ascending=False)
    gene_names = gene_names.reset_index()
    gene_names = gene_names.drop_duplicates(subset="GeneSymbol", keep='first')
    gene_names.set_index("GeneSymbol", inplace=True)
    gene_names = gene_names.drop(columns=['row_sum'])

    num_duplicates = gene_names.copy().index.duplicated().sum()
    st.write(num_duplicates)
    # find the duplicate index values
    duplicates = gene_names.index[gene_names.index.duplicated(keep=False)]

    # select the rows corresponding to the duplicate index values
    duplicate_rows = gene_names.loc[duplicates]
    st.write(duplicates)

    st.write("Shape of Mapped Dataframe:")
    st.write(gene_names.shape[0])
    st.write(gene_names.shape[1])

    # Download the corrected input DataFrame as a CSV file
    if st.button('Accept Correction'):
        gene_names.to_csv("/home/hcemm-user/SurvSig/source_data" + file_name, header=True,
                          index=True)


gene_names = cell_line_szandi()
st.write(gene_names)
name_dictionary, df = gene_autocorrect()
file_name = "lcnec_cell_mapped.csv"
app(name_dictionary, df, gene_names, file_name)
