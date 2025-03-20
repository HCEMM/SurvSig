import streamlit as st
import pandas as pd
import numpy as np
from Scripts.descriptions import rousseaux_surv_type, tcga_description


# Load George's SCLC Dataset and its Annotation
@st.cache_data
def george_sclc_data_load():
    """
    Load and process George's SCLC dataset and its annotation.

    Returns:
    - george_df (DataFrame): Processed gene expression data.
    - george_clinic (DataFrame): Processed clinical data.
    - george_st_clinic (DataFrame): Processed subset of clinical data.
    """
    try:
        george_df = pd.read_csv('source_data/george_mapped.csv', index_col=0)
        george_df = george_df + 1
        george_df = george_df.transform(lambda x: np.log2(x))
        pat_ids = list(george_df.transpose().index)

        # Read survival data
        george_clinic_original = pd.read_table(
            "source_data/sclc_ucologne_2015_clinical_data.tsv", index_col=1, delimiter="\t")

        # Rename columns and values
        george_clinic_original = george_clinic_original.rename(columns={"UICC Tumor Stage": "UICC"})
        george_clinic_original['Smoker'] = george_clinic_original['Smoker'].apply(
            lambda x: 'Former' if pd.notnull(x) and 'Former' in x else x)
        george_clinic_original = george_clinic_original.replace(
            ['Ib', 'Ia', 'IIa', 'IIb', 'IIIa', 'IIIb'],
            ['IB', 'IA', 'IIA', 'IIB', 'IIIA', 'IIIB'])
        george_clinic_original.sort_values(by=['UICC'], inplace=True)

        george_clinic = george_clinic_original[["Overall Survival (Months)", "Overall Survival Status"]]
        george_clinic = george_clinic.copy()  # Explicitly copy before modifying
        george_clinic.rename(
            columns={"Overall Survival (Months)": "Time", "Overall Survival Status": "Status"},
            inplace=True
        )

        george_st_clinic = george_clinic_original[['Sex', 'Smoker', 'UICC', 'Ethnicity Category']]
        george_st_clinic = george_st_clinic.copy()  # Explicitly copy before modifying
        george_st_clinic.fillna('NA', inplace=True)
        rename_mapping = {"NA": "NA_UICC"}
        george_st_clinic.loc[:, "UICC"] = george_st_clinic["UICC"].replace(rename_mapping)
        rename_mapping = {"NA": "NA_Smoker"}
        george_st_clinic.loc[:, "Smoker"] = george_st_clinic["Smoker"].replace(rename_mapping)
        rename_mapping = {"NA": "NA_EC"}
        george_st_clinic.loc[:, "Ethnicity Category"] = george_st_clinic["Ethnicity Category"].replace(rename_mapping)

        # Filter names based on patient IDs
        george_clinic = george_clinic[george_clinic.index.isin(pat_ids)]

        # Drop NaN values
        george_clinic = george_clinic.dropna(axis=0)

        # Encode "Status" -> 1: alive, 0: dead
        george_clinic = george_clinic.replace(['0:LIVING', '1:DECEASED'], ['0', '1'])
        george_clinic = george_clinic.astype("int")

        george_df.index = george_df.index.str.upper()

        return george_df, george_clinic, george_st_clinic
    except Exception as e:
        st.error(f"Error loading George's SCLC data: {e}")
        return None, None, None


@st.cache_data(show_spinner="Dataset is Loading...", persist=True)
def anish_sclc_data_load():
    """
    Load and process Anish's SCLC dataset and its annotation.

    Returns:
    - anish_df (DataFrame): Processed gene expression data.
    - anish_clinic (DataFrame): Processed clinical data.
    - anish_st_clinic (DataFrame): Processed subset of clinical data.
    """
    try:
        anish_df = pd.read_csv('source_data/anish_mapped.csv', index_col=0)
        pat_ids = list(anish_df.transpose().index)

        # Read survival data
        anish_clinic_original = pd.read_excel(
            "source_data/Supplementary_table_patient_clinical_characteristics_2021015.xlsx", index_col=2)
        anish_clinic = anish_clinic_original[["Overall survival since diagnosis", "Mortality status"]]
        anish_clinic.rename(
            columns={"Overall survival since diagnosis": "Time", "Mortality status": "Status"},
            inplace=True
        )
        anish_st_clinic = anish_clinic_original[['Sex', 'Race', 'Smoking history', 'SCLC staging/initial diagnosis ', 'Disease']]

        anish_st_clinic.fillna('NA', inplace=True)
        rename_mapping = {"NA": "NA_Sex"}
        anish_st_clinic["Sex"] = anish_st_clinic["Sex"].replace(rename_mapping)
        rename_mapping = {"NA": "NA_Race"}
        anish_st_clinic["Race"] = anish_st_clinic["Race"].replace(rename_mapping)
        rename_mapping = {"NA": "NA_Smoking"}
        anish_st_clinic["Smoking history"] = anish_st_clinic["Smoking history"].replace(rename_mapping)

        # Filter names based on patient IDs
        anish_clinic = anish_clinic[anish_clinic.index.isin(pat_ids)]

        # Drop NaN values
        anish_clinic = anish_clinic.dropna(axis=0)

        # Encode "Status" -> 1: dead, 0: alive
        anish_clinic = anish_clinic.replace(['Dead', 'Alive'], ['1', '0'])
        anish_clinic = anish_clinic.astype("int")

        anish_df.index = anish_df.index.str.upper()

        return anish_df, anish_clinic, anish_st_clinic
    except Exception as e:
        st.error(f"Error loading Anish's SCLC data: {e}")
        return None, None, None


@st.cache_data(show_spinner="Dataset is Loading...", persist=True)
def data_load_anno_george(selected_genes_anno):
    """
    Load and filter annotation data for George's SCLC dataset.

    Parameters:
    - selected_genes_anno (DataFrame): Selected genes annotation data.

    Returns:
    - df_anno_filtered (DataFrame): Filtered annotation data.
    """
    try:
        # Load the annotation dataframe from a CSV file
        df_anno_filtered = pd.read_csv("source_data/george_sclc_anno.csv", header=0)
        df_anno_filtered = df_anno_filtered.rename(columns={"Unnamed: 0": "Names"})
        df_anno_filtered.set_index("Names", inplace=True)

        # Save the sample names from the selected genes annotation dataframe
        sample_names = list(selected_genes_anno.index)

        # Reindex the annotation dataframe using only the common sample names
        df_anno_filtered = df_anno_filtered.reindex(index=sample_names)

        # Remove any sample names from the annotation dataframe that are not in the gene expression dataframe
        df_anno_filtered = df_anno_filtered.loc[df_anno_filtered.index.isin(sample_names)]

        return df_anno_filtered
    except Exception as e:
        st.error(f"Error loading annotation data for George's SCLC: {e}")
        return None


@st.cache_data(show_spinner="Dataset is Loading...", persist=True)
def data_load_anno_anish(selected_genes_anno):
    """
    Load and filter annotation data for Anish's SCLC dataset.

    Parameters:
    - selected_genes_anno (DataFrame): Selected genes annotation data.

    Returns:
    - df_anno_filtered (DataFrame): Filtered annotation data.
    """
    try:
        # Load the annotation dataframe from a CSV file
        df_anno_filtered = pd.read_csv("source_data/anish_anno.csv", index_col=0)
        df_anno_filtered = df_anno_filtered[["NAPY", "NE", 'EMT']]

        # Save the sample names from the selected genes annotation dataframe
        sample_names = list(selected_genes_anno.index)

        # Reindex the annotation dataframe using only the common sample names
        df_anno_filtered = df_anno_filtered.reindex(index=sample_names)

        # Remove any sample names from the annotation dataframe that are not in the gene expression dataframe
        df_anno_filtered = df_anno_filtered.loc[df_anno_filtered.index.isin(sample_names)]

        return df_anno_filtered
    except Exception as e:
        st.error(f"Error loading annotation data for Anish's SCLC: {e}")
        return None


@st.cache_data(show_spinner="TCGA Dataset is Loading", persist=True)
def data_read():
    """
    Load and process TCGA dataset.

    Returns:
    - tcga_df (DataFrame): Processed TCGA gene expression data.
    - sample_type (DataFrame): Sample types.
    - cancer_types_list (list): List of cancer types.
    - tcga_cancer_types (DataFrame): TCGA cancer types data.
    - tcga_subtypes (DataFrame): TCGA subtypes data.
    - tcga_surv (DataFrame): TCGA survival data.
    """
    try:
        tcga_df = pd.read_csv("source_data/tcga_mapped_fillna.csv", index_col=0, header=0)
        tcga_df.index = tcga_df.index.str.upper()
        tcga_df = tcga_df + 1
        tcga_df = tcga_df.transform(lambda x: np.log2(x).clip(lower=0))

        tcga_subtypes = pd.read_csv("TCGA_prep/TCGAbiolinks_subtypes.tsv", delimiter="\t")
        tcga_subtypes["pan.samplesID"] = tcga_subtypes["pan.samplesID"].str[:12]
        tcga_subtypes["pan.samplesID"] = tcga_subtypes["pan.samplesID"].str.ljust(12)[12:]

        tcga_cancer_types = pd.read_csv("TCGA_prep/TCGA-CDR_clinical_data.tsv", delimiter="\t")

        if 'OS.time' in tcga_cancer_types.columns and pd.api.types.is_numeric_dtype(tcga_cancer_types['OS.time']):
            tcga_cancer_types['OS.time'] = tcga_cancer_types['OS.time'] / 30.4375  # Convert days to months
        if 'PFI.time' in tcga_cancer_types.columns and pd.api.types.is_numeric_dtype(tcga_cancer_types['PFI.time']):
            tcga_cancer_types['PFI.time'] = tcga_cancer_types['PFI.time'] / 30.4375
        if 'DSS.time' in tcga_cancer_types.columns and pd.api.types.is_numeric_dtype(tcga_cancer_types['DSS.time']):
            tcga_cancer_types['DSS.time'] = tcga_cancer_types['DSS.time'] / 30.4375
        if 'DFI.time' in tcga_cancer_types.columns and pd.api.types.is_numeric_dtype(tcga_cancer_types['DFI.time']):
            tcga_cancer_types['DFI.time'] = tcga_cancer_types['DFI.time'] / 30.4375

        tcga_cancer_types['clinical_stage'] = tcga_cancer_types['clinical_stage'].str.replace('Stage ', '')
        tcga_cancer_types = tcga_cancer_types.replace(['IIa', 'IIb', 'IVa', 'IVb'], ['IIA', 'IIB', 'IVA', 'IVB'])
        tcga_surv = pd.DataFrame()
        tcga_surv[["Patient ID", "OS.time", "OS", 'DSS', 'DSS.time', 'PFI', 'PFI.time', 'DFI', 'DFI.time']] = \
            tcga_cancer_types[["bcr_patient_barcode", "OS.time", "OS", 'DSS.time', 'DSS', 'PFI.time', 'PFI', 'DFI.time', 'DFI']]

        cancer_types_list = tcga_cancer_types["type"].copy().unique()
        sample_type = pd.DataFrame()
        sample_type["Sample_ID"] = tcga_cancer_types["bcr_patient_barcode"].copy()
        sample_type["Cancer_type"] = tcga_cancer_types["type"].copy()
        tcga_df = tcga_df.transpose()

        # Renaming supplementary dataframe names using mapping
        name_mapping_tcga = pd.read_excel("source_data/name_mapping_tcga.xlsx", dtype=str, header=None, index_col=0)
        name_mapping_tcga.drop_duplicates(keep='first', inplace=True)
        name_mapping_tcga.dropna(axis=1, how="all", inplace=True)
        name_mapping_tcga.fillna('nan', axis=1, inplace=True)
        tcga_dict_supp = dict(zip(name_mapping_tcga.index, name_mapping_tcga.iloc[:, 0]))

        tcga_subtypes.replace(tcga_dict_supp, inplace=True)
        tcga_cancer_types.replace(tcga_dict_supp, inplace=True)

        # Specified columns to work with
        columns_to_modify = [
            'Subtype_protein',
            'Subtype_mRNA',
            'Subtype_miRNA',
            'Subtype_Integrative',
            'Subtype_DNAmeth',
            'Subtype_CNA'
        ]

        # String numbers to check for
        string_numbers = ['1', '2', '3', '4', '5', '6', '7']

        # Loop through each column and apply the transformation
        for col in columns_to_modify:
            suffix = col.replace('Subtype_', '')  # Remove "Subtype_" from the column name for the suffix
            tcga_subtypes[col] = tcga_subtypes[col].apply(lambda x: transform_value(x, suffix, string_numbers))

        return tcga_df, sample_type, cancer_types_list, tcga_cancer_types, tcga_subtypes, tcga_surv
    except Exception as e:
        st.error(f"Error loading TCGA data: {e}")
        return None, None, None, None, None, None


def cancer_type_separator(tcga_df, sample_type, cancer_types_list, gene_sig_exp, key):
    """
    Separate data by cancer type.

    Parameters:
    - tcga_df (DataFrame): TCGA gene expression data.
    - sample_type (DataFrame): Sample types.
    - cancer_types_list (list): List of cancer types.
    - gene_sig_exp (Streamlit selector): Streamlit selector component.
    - key (str): Key for the selector.

    Returns:
    - selected_cancer_df (DataFrame): Selected cancer type data.
    - cancer_opt (str): Selected cancer type.
    """
    try:
        if key == "feature":
            cancer_opt = gene_sig_exp.selectbox(
                ":blue[Select a ] :red[Type of Cancer]", options=cancer_types_list, key=key,
                help=tcga_description()
            )
        else:
            cancer_opt = st.sidebar.selectbox(
                ":blue[Select a ] :red[Type of Cancer]", options=cancer_types_list, key=key,
                help=tcga_description()
            )
        st.sidebar.header(" ", divider="blue")
        globals()[f"TCGA_{cancer_opt}"] = sample_type[sample_type['Cancer_type'] == cancer_opt]
        globals()[f"TCGA_{cancer_opt}"].set_index("Sample_ID", inplace=True)

        selected_cancer_df = tcga_df[tcga_df.index.isin(globals()[f"TCGA_{cancer_opt}"].index)]
        selected_cancer_df.dropna(axis=1, inplace=True, how='all')

        selected_cancer_df = selected_cancer_df.transpose()

        del sample_type, cancer_types_list

        selected_cancer_df.index = selected_cancer_df.index.str.upper()

        return selected_cancer_df, cancer_opt
    except Exception as e:
        st.error(f"Error in cancer type separator: {e}")
        return None, None


@st.cache_data(show_spinner="Dataset is Loading...", persist=True)
def jiang_sclc_data_load():
    """
    Load and process Jiang's SCLC dataset and its annotation.

    Returns:
    - jiang_df (DataFrame): Processed gene expression data.
    - jiang_clinic (DataFrame): Processed clinical data.
    - jiang_st_clinic (DataFrame): Processed subset of clinical data.
    """
    try:
        jiang_df = pd.read_csv('source_data/jiang_mapped.csv', index_col=0)
        jiang_df.index.rename("Sample_id", inplace=True)
        pat_ids_j = list(jiang_df.transpose().index)

        # Read survival data
        jiang_clinic_original = pd.read_csv("source_data/Jiang_SCLC_annotation.csv", index_col=0)
        jiang_clinic_original.rename(columns={'chemotherapy-exposed biopsy': 'Chemotherapy'}, inplace=True)
        jiang_clinic = jiang_clinic_original[["Months", "Survival Status"]]
        jiang_clinic.rename(columns={"Months": "Time", "Survival Status": "Status"}, inplace=True)
        jiang_st_clinic = jiang_clinic_original[['Gender', 'Smoke', 'UICC stage', 'Chemotherapy', 'EMT', 'NE', 'NAPY']]

        # Filter names based on patient IDs
        jiang_clinic = jiang_clinic[jiang_clinic.index.isin(pat_ids_j)]

        # Drop NaN values
        jiang_clinic = jiang_clinic.dropna(axis=0)
        jiang_df.index = jiang_df.index.str.upper()

        return jiang_df, jiang_clinic, jiang_st_clinic
    except Exception as e:
        st.error(f"Error loading Jiang's SCLC data: {e}")
        return None, None, None


@st.cache_data(show_spinner="Dataset is Loading...", persist=True)
def george_lcnec_data_load():
    """
    Load and process George's LCNEC dataset and its annotation.

    Returns:
    - george_lcnec_df (DataFrame): Processed gene expression data.
    - george_lcnec_clinic (DataFrame): Processed clinical data.
    - george_lcnec_st_clinic (DataFrame): Processed subset of clinical data.
    """
    try:
        george_lcnec_df = pd.read_csv('source_data/lcnec_mapped.csv', index_col=0)
        george_lcnec_df.columns = george_lcnec_df.columns.str.replace("LCNEC_", "", regex=True)
        george_lcnec_df.index.rename("Sample_id", inplace=True)
        george_lcnec_df = george_lcnec_df + 1
        george_lcnec_df = george_lcnec_df.transform(lambda x: np.log2(x))
        pat_ids_lcnec = list(george_lcnec_df.columns)

        # Read survival data
        george_lcnec_clinic_original = pd.read_csv("source_data/LCNEC_George_clinical_anno.csv", index_col=0)
        george_lcnec_clinic = george_lcnec_clinic_original[["Survival_months", "Survival_censor"]]
        george_lcnec_clinic.rename(columns={'Survival_months': 'Time', 'Survival_censor': 'Status'}, inplace=True)

        george_lcnec_st_clinic = george_lcnec_clinic_original[['Sex', 'Smoking_Hx', 'Stage_UICC', 'EMT', 'NE', 'NAPY', 'Classification', 'Molecular_Subtypes']]
        george_lcnec_st_clinic.rename(columns={'Sex': 'Gender', 'Smoking_Hx': 'Smoke', 'Stage_UICC': 'UICC Stage'}, inplace=True)

        george_lcnec_st_clinic.fillna('NA', inplace=True)
        rename_mapping = {"NA": "NA_UICC"}
        george_lcnec_st_clinic["UICC Stage"] = george_lcnec_st_clinic["UICC Stage"].replace(rename_mapping)
        rename_mapping = {"NA": "NA_Smoke"}
        george_lcnec_st_clinic["Smoke"] = george_lcnec_st_clinic["Smoke"].replace(rename_mapping)
        rename_mapping = {"NA": "NA_Mol_ST"}
        george_lcnec_st_clinic["Molecular_Subtypes"] = george_lcnec_st_clinic["Molecular_Subtypes"].replace(rename_mapping)

        # Filter names based on patient IDs
        george_lcnec_clinic = george_lcnec_clinic[george_lcnec_clinic.index.isin(pat_ids_lcnec)]

        # Drop NaN values
        george_lcnec_clinic = george_lcnec_clinic.dropna(axis=0, how="all")
        george_lcnec_df.index = george_lcnec_df.index.str.upper()
        george_lcnec_clinic = george_lcnec_clinic.replace(['dead', 'alive'], ['1', '0'])
        george_lcnec_clinic["Status"] = george_lcnec_clinic["Status"].astype("float").astype('int')

        george_lcnec_df.index = george_lcnec_df.index.str.upper()

        return george_lcnec_df, george_lcnec_clinic, george_lcnec_st_clinic
    except Exception as e:
        st.error(f"Error loading George's LCNEC data: {e}")
        return None, None, None


@st.cache_data(show_spinner="Dataset is Loading...", persist=True)
def carcinoid_data_load():
    """
    Load and process carcinoid dataset and its annotation.

    Returns:
    - carcinoid_df (DataFrame): Processed gene expression data.
    - carcinoid_clinic (DataFrame): Processed clinical data.
    - carcinoid_st_clinic (DataFrame): Processed subset of clinical data.
    """
    try:
        carcinoid_df = pd.read_csv('source_data/carcinoid_mapped.csv', index_col=0)
        carcinoid_df.index.rename("Sample_id", inplace=True)
        carcinoid_df = carcinoid_df + 1
        carcinoid_df = carcinoid_df.transform(lambda x: np.log2(x))
        pat_ids_carci = list(carcinoid_df.columns)

        carcinoid_clinic_original = pd.read_csv("source_data/carcinoid_clinical_anno.csv", index_col=0)
        carcinoid_clinic = carcinoid_clinic_original[["Survival_months", "Survival_censor"]]
        carcinoid_clinic.rename(columns={'Survival_months': 'Time', 'Survival_censor': 'Status'}, inplace=True)
        carcinoid_st_clinic = carcinoid_clinic_original[['Sex', 'Smoking_status', 'Stage_UICC', 'EMT', 'NE', 'NAPY', 'cluster_LNET', 'cluster_LNEN', 'Histopathology']]
        carcinoid_st_clinic.rename(columns={'Sex': 'Gender', 'Smoking_status': 'Smoke', 'Stage_UICC': 'UICC Stage'}, inplace=True)

        carcinoid_st_clinic.fillna('NA', inplace=True)
        rename_mapping = {"NA": "NA_UICC"}
        carcinoid_st_clinic["UICC Stage"] = carcinoid_st_clinic["UICC Stage"].replace(rename_mapping)
        rename_mapping = {"NA": "NA_Smoke"}
        carcinoid_st_clinic["Smoke"] = carcinoid_st_clinic["Smoke"].replace(rename_mapping)

        # Filter names based on patient IDs
        carcinoid_clinic = carcinoid_clinic[carcinoid_clinic.index.isin(pat_ids_carci)]

        # Drop NaN values
        carcinoid_clinic = carcinoid_clinic.dropna(axis=0, how="all")
        carcinoid_clinic.index = carcinoid_clinic.index.str.upper()
        carcinoid_clinic = carcinoid_clinic.replace(['dead', 'alive'], ['1', '0'])
        carcinoid_clinic["Status"] = carcinoid_clinic["Status"].astype("float").astype('int')

        carcinoid_df.index = carcinoid_df.index.str.upper()

        return carcinoid_df, carcinoid_clinic, carcinoid_st_clinic
    except Exception as e:
        st.error(f"Error loading carcinoid data: {e}")
        return None, None, None


@st.cache_data(show_spinner="Dataset is Loading...", persist=True)
def alcala_data_load(exclusion):
    """
    Load and process Alcala dataset and its annotation.

    Parameters:
    - exclusion (bool): Flag to exclude specific samples.

    Returns:
    - carcinoid_df (DataFrame): Processed gene expression data.
    - carcinoid_clinic (DataFrame): Processed clinical data.
    - carcinoid_st_clinic (DataFrame): Processed subset of clinical data.
    """
    try:
        carcinoid_df = pd.read_csv('source_data/alcala_mapped.csv', index_col=0)
        carcinoid_df.index.rename("Sample_id", inplace=True)

        if exclusion:
            carcinoid_df = carcinoid_df.loc[:, ~carcinoid_df.columns.str.startswith('SRR')]
        pat_ids_alcala = list(carcinoid_df.columns)

        carcinoid_clinic_original = pd.read_csv("source_data/alcala_clinical_anno.csv", index_col=0)
        carcinoid_clinic = carcinoid_clinic_original[["Survival_months", "Survival_censor"]]
        carcinoid_clinic.rename(columns={'Survival_months': 'Time', 'Survival_censor': 'Status'}, inplace=True)
        carcinoid_st_clinic = carcinoid_clinic_original[['Sex', 'Smoking_status', 'EMT', 'NE', 'NAPY', 'Histopathology_simplified', 'Molecular_clusters', 'Histopathology']]
        carcinoid_st_clinic.rename(columns={'Sex': 'Gender', 'Smoking_status': 'Smoke'}, inplace=True)

        rename_mapping = {"Supra_carcinoid": "Supra_carcinoid_HS", "Atypical": "Atypical_HS", "Typical": "Typical_HS"}
        carcinoid_st_clinic["Histopathology_simplified"] = carcinoid_st_clinic["Histopathology_simplified"].replace(rename_mapping)

        rename_mapping = {"Supra_carcinoid": "Supra_carcinoid_MC"}
        carcinoid_st_clinic["Molecular_clusters"] = carcinoid_st_clinic["Molecular_clusters"].replace(rename_mapping)

        conditions = [
            carcinoid_st_clinic['Molecular_clusters'].isin(['Carcinoid-A2', 'LC3']),
            carcinoid_st_clinic['Molecular_clusters'].isin(['Carcinoid-A1', 'LC1']),
            carcinoid_st_clinic['Molecular_clusters'].isin(['Carcinoid-B', 'LC2'])
        ]

        choices = ['A2_mol_clust', 'A1_mol_sclust', 'B_mol_clust']

        # Create a new column using numpy.select
        carcinoid_st_clinic['Merged_Mol_clusters'] = np.select(conditions, choices, default='Supra_carcinoid_mol_clust')

        # Filter names based on patient IDs
        carcinoid_clinic = carcinoid_clinic[carcinoid_clinic.index.isin(pat_ids_alcala)]
        carcinoid_st_clinic = carcinoid_st_clinic[carcinoid_st_clinic.index.isin(pat_ids_alcala)]

        # Drop NaN values
        carcinoid_clinic = carcinoid_clinic.dropna(axis=0, how="all")
        carcinoid_clinic.index = carcinoid_clinic.index.str.upper()
        carcinoid_clinic = carcinoid_clinic.replace(['dead', 'alive'], ['1', '0'])
        carcinoid_clinic["Status"] = carcinoid_clinic["Status"].astype("float").astype('int')

        carcinoid_df.index = carcinoid_df.index.str.upper()

        return carcinoid_df, carcinoid_clinic, carcinoid_st_clinic
    except Exception as e:
        st.error(f"Error loading Alcala data: {e}")
        return None, None, None


def rousseaux_surv_type_selector():
    """
    Select survival type for Rousseaux dataset.

    Returns:
    - survival_type (str): Selected survival type.
    """
    try:
        # Display a selectbox in the sidebar for selecting the type of survival time
        survival_type = st.sidebar.selectbox(
            ':blue[Please Select the ] :red[Type of Survival Time]', index=0,
            options=["Overall", "Disease Free Interval"], key="surv_type_rousseaux",
            help=rousseaux_surv_type()
        )
        st.sidebar.subheader(" ", divider="blue")

        return survival_type
    except Exception as e:
        st.error(f"An error occurred while selecting the survival type: {e}")
        return None


@st.cache_data(show_spinner="Dataset is Loading...", persist=True)
def rousseaux_data_load(exclusion, survival_type):
    """
    Load and process Rousseaux dataset and its annotation.

    Parameters:
    - exclusion (bool): Flag to exclude specific samples.
    - survival_type (str): Type of survival time to use.

    Returns:
    - rousseaux_df (DataFrame): Processed gene expression data.
    - rousseaux_clinic (DataFrame): Processed clinical data.
    - rousseaux_st_clinic (DataFrame): Processed subset of clinical data.
    """
    try:
        rousseaux_df = pd.read_csv('source_data/rousseaux_mapped.csv', index_col=0)
        rousseaux_df.index.rename("Sample_id", inplace=True)

        rousseaux_clinic_original = pd.read_csv("source_data/rousseaux_clinical_anno.csv", index_col=0)
        rousseaux_clinic_original.rename(columns={"T": "T_Stage"}, inplace=True)

        if exclusion:
            samples_to_exclude = rousseaux_clinic_original[rousseaux_clinic_original['histology'] == 'NTL'].index
            rousseaux_df = rousseaux_df.loc[:, ~rousseaux_df.columns.isin(samples_to_exclude)]

        pat_ids_rousseaux = list(rousseaux_df.columns)

        rousseaux_clinic = rousseaux_clinic_original[["os_surv", "os_status", "relapse_surv", "relapse_status"]]

        if survival_type == "Disease Free Interval":
            rousseaux_clinic.drop(columns=["os_surv", "os_status"], inplace=True)
            rousseaux_clinic.rename(columns={'relapse_status': 'Status', 'relapse_surv': 'Time'}, inplace=True)
        else:
            rousseaux_clinic.drop(columns=["relapse_surv", "relapse_status"], inplace=True)
            rousseaux_clinic.rename(columns={'os_status': 'Status', 'os_surv': 'Time'}, inplace=True)

        rousseaux_st_clinic = rousseaux_clinic_original[['gender', 'EMT', 'NE', 'NAPY', 'histology', 'TNM', 'T_Stage', 'N', 'M', 'SurvSig_Mol.Subtype']]
        rousseaux_st_clinic.rename(columns={'gender': 'Gender', 'histology': 'Histology'}, inplace=True)

        rename_mapping = {"NTL": "NTL_Gender"}
        rousseaux_st_clinic["Gender"] = rousseaux_st_clinic["Gender"].replace(rename_mapping)
        rename_mapping = {"NTL": "NTL_TNM"}
        rousseaux_st_clinic["TNM"] = rousseaux_st_clinic["TNM"].replace(rename_mapping)
        rename_mapping = {"NTL": "NTL_T_Stage"}
        rousseaux_st_clinic["T_Stage"] = rousseaux_st_clinic["T_Stage"].replace(rename_mapping)
        rename_mapping = {"NTL": "NTL_N"}
        rousseaux_st_clinic["N"] = rousseaux_st_clinic["N"].replace(rename_mapping)
        rename_mapping = {"NTL": "NTL_M"}
        rousseaux_st_clinic["M"] = rousseaux_st_clinic["M"].replace(rename_mapping)

        rousseaux_st_clinic["TNM"].fillna('NA', inplace=True)
        rename_mapping = {'NA': "NTL_TNM"}
        rousseaux_st_clinic["TNM"] = rousseaux_st_clinic["TNM"].replace(rename_mapping)

        rousseaux_st_clinic["SurvSig_Mol.Subtype"].fillna('NA', inplace=True)
        rename_mapping = {'NA': "non_CARCI"}
        rousseaux_st_clinic["SurvSig_Mol.Subtype"] = rousseaux_st_clinic["SurvSig_Mol.Subtype"].replace(rename_mapping)

        rousseaux_clinic = rousseaux_clinic[rousseaux_clinic.index.isin(pat_ids_rousseaux)]
        rousseaux_st_clinic = rousseaux_st_clinic[rousseaux_st_clinic.index.isin(pat_ids_rousseaux)]

        rousseaux_clinic = rousseaux_clinic.dropna(axis=0, how="all")
        rousseaux_clinic.index = rousseaux_clinic.index.str.upper()
        rousseaux_clinic["Status"] = rousseaux_clinic["Status"].astype("float").astype('int')

        rousseaux_df.index = rousseaux_df.index.str.upper()

        return rousseaux_df, rousseaux_clinic, rousseaux_st_clinic
    except Exception as e:
        st.error(f"Error loading Rousseaux data: {e}")
        return None, None, None


@st.cache_data(show_spinner="Dataset is Loading...", persist=True)
def data_load_anno_rousseaux(selected_genes_anno):
    """
    Load and filter annotation data for Rousseaux dataset.

    Parameters:
    - selected_genes_anno (DataFrame): Selected genes annotation data.

    Returns:
    - df_anno_filtered (DataFrame): Filtered annotation data.
    """
    try:
        df_anno_filtered = pd.read_csv('source_data/rousseaux_clinical_anno.csv', index_col=0)
        df_anno_filtered.rename(columns={"T": "T_Stage"}, inplace=True)
        df_anno_filtered = df_anno_filtered[["NAPY", "NE", "EMT", 'TNM', 'T_Stage', 'N', 'M','SurvSig_Mol.Subtype']]

        df_anno_filtered["TNM"].fillna('NA', inplace=True)
        rename_mapping = {'NA': "NTL_TNM"}
        df_anno_filtered["TNM"] = df_anno_filtered["TNM"].replace(rename_mapping)

        df_anno_filtered["SurvSig_Mol.Subtype"].fillna('NA', inplace=True)
        rename_mapping = {'NA': "non_LCNEC"}
        df_anno_filtered["SurvSig_Mol.Subtype"] = df_anno_filtered["SurvSig_Mol.Subtype"].replace(rename_mapping)

        rename_mapping = {"NTL": "NTL_T_Stage"}
        df_anno_filtered["T_Stage"] = df_anno_filtered["T_Stage"].replace(rename_mapping)

        rename_mapping = {"NTL": "NTL_N"}
        df_anno_filtered["N"] = df_anno_filtered["N"].replace(rename_mapping)

        rename_mapping = {"NTL": "NTL_M"}
        df_anno_filtered["M"] = df_anno_filtered["M"].replace(rename_mapping)

        sample_names = list(selected_genes_anno.index)

        df_anno_filtered = df_anno_filtered.reindex(index=sample_names)
        df_anno_filtered = df_anno_filtered.loc[df_anno_filtered.index.isin(sample_names)]

        return df_anno_filtered
    except Exception as e:
        st.error(f"Error loading annotation data for Rousseaux: {e}")
        return None


def convert_string_to_float(value):
    """
    Convert string to float, replacing commas with dots if necessary.

    Parameters:
    - value (str): Input string value.

    Returns:
    - float: Converted float value.
    """
    try:
        if isinstance(value, str) and ',' in value:
            return float(value.replace(',', '.'))
        else:
            return value
    except ValueError as e:
        st.error(f"Error converting string to float: {e}")
        return value


@st.cache_data(show_spinner="Dataset is Loading...", persist=True)
def liu_data_load():
    """
    Load and process Liu's dataset and its annotation.

    Returns:
    - liu_df (DataFrame): Processed gene expression data.
    - liu_clinic (DataFrame): Processed clinical data.
    - liu_st_clinic (DataFrame): Processed subset of clinical data.
    """
    try:
        liu_df = pd.read_csv('source_data/liu_mapped.csv', index_col=0)
        liu_df.index.rename("Sample_id", inplace=True)
        liu_df.index = liu_df.index.astype('str')
        pat_ids_liu = list(liu_df.columns)

        liu_clinic_original = pd.read_csv("source_data/liu_clinical_anno.csv", index_col=0)
        liu_clinic_original.index = liu_clinic_original.index.astype('str')

        liu_clinic = liu_clinic_original[["Survial.(months)", "Status."]]
        liu_clinic.rename(columns={'Survial.(months)': 'Time', 'Status.': 'Status'}, inplace=True)
        liu_st_clinic = liu_clinic_original[['Gender', 'Smoking', 'EMT', 'NE', 'NAPY', 'Tumor.site', 'Histologic.type', 'TNM.Stage']]
        liu_st_clinic.rename(columns={'Smoking': 'Smoke', 'TNM.Stage': 'TNM_Stage', 'Tumor.site': 'Tumor_site', 'Histologic.type': 'Histologic_type'}, inplace=True)

        liu_clinic = liu_clinic[liu_clinic.index.isin(pat_ids_liu)]
        liu_st_clinic = liu_st_clinic[liu_st_clinic.index.isin(pat_ids_liu)]

        liu_clinic = liu_clinic.dropna(axis=0, how="all")
        liu_clinic.index = liu_clinic.index.str.upper()
        liu_clinic['Time'] = liu_clinic['Time'].str.split(r' \(').str[0]
        liu_clinic = liu_clinic.applymap(convert_string_to_float)
        liu_clinic['Time'] = liu_clinic['Time'].astype("float")
        liu_clinic = liu_clinic.replace(['dead', 'alive'], ['1', '0'])
        liu_clinic["Status"] = liu_clinic["Status"].astype("float").astype('int')
        liu_df.index = liu_df.index.str.upper()

        return liu_df, liu_clinic, liu_st_clinic
    except Exception as e:
        st.error(f"Error loading Liu's data: {e}")
        return None, None, None


@st.cache_data(show_spinner="Dataset is Loading...", persist=True)
def data_load_anno_liu(selected_genes_anno):
    """
    Load and filter annotation data for Liu's dataset.

    Parameters:
    - selected_genes_anno (DataFrame): Selected genes annotation data.

    Returns:
    - df_anno_filtered (DataFrame): Filtered annotation data.
    """
    try:
        df_anno_filtered = pd.read_csv('source_data/liu_clinical_anno.csv', index_col=0)
        df_anno_filtered.index = df_anno_filtered.index.astype('str')
        df_anno_filtered = df_anno_filtered[['Gender', 'Smoking', 'EMT', 'NE', 'NAPY', 'Tumor.site', 'Histologic.type', 'TNM.Stage']]
        df_anno_filtered.rename(columns={'Smoking': 'Smoke', 'TNM.Stage': 'TNM_Stage', 'Tumor.site': 'Tumor_site', 'Histologic.type': 'Histologic_type'}, inplace=True)

        sample_names = list(selected_genes_anno.index)

        df_anno_filtered = df_anno_filtered.reindex(index=sample_names)
        df_anno_filtered = df_anno_filtered.loc[df_anno_filtered.index.isin(sample_names)]

        return df_anno_filtered
    except Exception as e:
        st.error(f"Error loading annotation data for Liu: {e}")
        return None


@st.cache_data(show_spinner="Dataset is Loading...", persist=True)
def data_load_anno_alcala(selected_genes_anno):
    """
    Load and filter annotation data for Alcala dataset.

    Parameters:
    - selected_genes_anno (DataFrame): Selected genes annotation data.

    Returns:
    - df_anno_filtered (DataFrame): Filtered annotation data.
    """
    try:
        df_anno_filtered = pd.read_csv('source_data/alcala_clinical_anno.csv', index_col=0)
        df_anno_filtered = df_anno_filtered[["NAPY", "NE", "EMT", 'Histopathology_simplified', 'Molecular_clusters', 'Histopathology']]

        rename_mapping = {"Supra_carcinoid": "Supra_carcinoid_HS", "Atypical": "Atypical_HS", "Typical": "Typical_HS"}
        df_anno_filtered["Histopathology_simplified"] = df_anno_filtered["Histopathology_simplified"].replace(rename_mapping)

        rename_mapping = {"Supra_carcinoid": "Supra_carcinoid_MC"}
        df_anno_filtered["Molecular_clusters"] = df_anno_filtered["Molecular_clusters"].replace(rename_mapping)

        conditions = [
            df_anno_filtered['Molecular_clusters'].isin(['Carcinoid-A2', 'LC3']),
            df_anno_filtered['Molecular_clusters'].isin(['Carcinoid-A1', 'LC1']),
            df_anno_filtered['Molecular_clusters'].isin(['Carcinoid-B', 'LC2'])
        ]

        choices = ['A2_mol_clust', 'A1_mol_sclust', 'B_mol_clust']

        # Create a new column using numpy.select
        df_anno_filtered['Merged_Mol_clusters'] = np.select(conditions, choices, default='carcinoid_st_clinic')

        sample_names = list(selected_genes_anno.index)

        df_anno_filtered = df_anno_filtered.reindex(index=sample_names)
        df_anno_filtered = df_anno_filtered.loc[df_anno_filtered.index.isin(sample_names)]

        return df_anno_filtered
    except Exception as e:
        st.error(f"Error loading annotation data for Alcala: {e}")
        return None


@st.cache_data(show_spinner="Dataset is Loading...", persist=True)
def data_load_anno_carcinoid(selected_genes_anno):
    """
    Load and filter annotation data for carcinoid dataset.

    Parameters:
    - selected_genes_anno (DataFrame): Selected genes annotation data.

    Returns:
    - df_anno_filtered (DataFrame): Filtered annotation data.
    """
    try:
        df_anno_filtered = pd.read_csv('source_data/carcinoid_clinical_anno.csv', index_col=0)
        df_anno_filtered = df_anno_filtered[["NAPY", "NE", "EMT", 'cluster_LNET', 'cluster_LNEN', 'Histopathology']]

        sample_names = list(selected_genes_anno.index)

        df_anno_filtered = df_anno_filtered.reindex(index=sample_names)
        df_anno_filtered = df_anno_filtered.loc[df_anno_filtered.index.isin(sample_names)]

        return df_anno_filtered
    except Exception as e:
        st.error(f"Error loading annotation data for carcinoid: {e}")
        return None


@st.cache_data(show_spinner="Dataset is Loading...", persist=True)
def data_load_anno_george_lcnec(selected_genes_anno):
    """
    Load and filter annotation data for George's LCNEC dataset.

    Parameters:
    - selected_genes_anno (DataFrame): Selected genes annotation data.

    Returns:
    - df_anno_filtered (DataFrame): Filtered annotation data.
    """
    try:
        df_anno_filtered = pd.read_csv('source_data/LCNEC_George_clinical_anno.csv', index_col=0)
        df_anno_filtered = df_anno_filtered[["NAPY", "NE", "EMT", "Classification", "Molecular_Subtypes"]]

        df_anno_filtered.fillna('NA', inplace=True)
        rename_mapping = {"NA": "NA_Mol_ST"}
        df_anno_filtered["Molecular_Subtypes"] = df_anno_filtered["Molecular_Subtypes"].replace(rename_mapping)

        sample_names = list(selected_genes_anno.index)

        df_anno_filtered = df_anno_filtered.reindex(index=sample_names)
        df_anno_filtered = df_anno_filtered.loc[df_anno_filtered.index.isin(sample_names)]

        return df_anno_filtered
    except Exception as e:
        st.error(f"Error loading annotation data for George's LCNEC: {e}")
        return None


@st.cache_data(show_spinner="Dataset is Loading...", persist=True)
def data_load_anno_jiang(selected_genes_anno):
    """
    Load and filter annotation data for Jiang's SCLC dataset.

    Parameters:
    - selected_genes_anno (DataFrame): Selected genes annotation data.

    Returns:
    - df_anno_filtered (DataFrame): Filtered annotation data.
    """
    try:
        df_anno_filtered = pd.read_csv('source_data/Jiang_SCLC_annotation.csv', index_col=0)
        df_anno_filtered = df_anno_filtered[["NAPY", "NE", "EMT"]]

        sample_names = list(selected_genes_anno.index)

        df_anno_filtered = df_anno_filtered.reindex(index=sample_names)
        df_anno_filtered = df_anno_filtered.loc[df_anno_filtered.index.isin(sample_names)]

        return df_anno_filtered
    except Exception as e:
        st.error(f"Error loading annotation data for Jiang's SCLC: {e}")
        return None


def tcga_surv_time_corrector(tcga_surv, target_cancer_exp_df, tcga_clinics):
    """
    Correct survival time data for TCGA dataset.

    Parameters:
    - tcga_surv (DataFrame): TCGA survival data.
    - target_cancer_exp_df (DataFrame): Target cancer expression data.
    - tcga_clinics (DataFrame): TCGA clinics data.

    Returns:
    - tcga_surv (DataFrame): Corrected TCGA survival data.
    - tcga_clinics (DataFrame): Corrected TCGA clinics data.
    """
    try:
        tcga_surv.set_index("Patient ID", inplace=True)
        tcga_surv = tcga_surv[tcga_surv.index.isin(target_cancer_exp_df.columns)]
        tcga_surv.dropna(axis=1, how='all', inplace=True)

        return tcga_surv, tcga_clinics
    except Exception as e:
        st.error(f"Error in TCGA survival time corrector: {e}")
        return None, None


def tcga_survival_selector(tcga_surv, tcga_clinics):
    """
    Select survival type for TCGA dataset.

    Parameters:
    - tcga_surv (DataFrame): TCGA survival data.
    - tcga_clinics (DataFrame): TCGA clinics data.

    Returns:
    - tcga_surv (DataFrame): Selected survival data.
    - tcga_clinics (DataFrame): Updated clinics data.
    """
    try:
        # Define columns to check and corresponding items to add
        columns_and_items = {
            'OS': 'Overall',
            'DFI': 'Disease-free Interval',
            'PFI': 'Progression-free Interval',
            'DSS': 'Disease Specific'
        }
        # Initialize an empty list for survival options
        items_to_add = []

        # Iterate through the columns and add corresponding item if column exists
        for column, item in columns_and_items.items():
            if column in tcga_surv.columns:
                items_to_add.append(item)

        survival_type = st.sidebar.selectbox(
            ':blue[Please Select the ] :red[Type of Survival Time]', index=0,
            options=items_to_add, key="surv_type"
        )
        st.sidebar.subheader(" ", divider="blue")

        if survival_type == "Disease Specific":
            tcga_clinics.drop(columns=["OS", "OS.time", 'DFI', 'DFI.time', 'PFI', 'PFI.time'], inplace=True)
            tcga_clinics.rename(columns={'DSS': 'Status', 'DSS.time': 'Time'}, inplace=True)
            columns_to_drop = ["OS", "OS.time", 'DFI', 'DFI.time', 'PFI', 'PFI.time']
            tcga_surv = tcga_surv.drop(columns=[col for col in columns_to_drop if col in tcga_surv.columns])
            tcga_surv.rename(columns={'DSS': 'Time', 'DSS.time': 'Status'}, inplace=True)

        elif survival_type == "Progression-free Interval":
            tcga_clinics.drop(columns=["OS", "OS.time", 'DFI', 'DFI.time', 'DSS', 'DSS.time'], inplace=True)
            tcga_clinics.rename(columns={'PFI': 'Status', 'PFI.time': 'Time'}, inplace=True)
            columns_to_drop = ["OS", "OS.time", 'DFI', 'DFI.time', 'DSS', 'DSS.time']
            tcga_surv = tcga_surv.drop(columns=[col for col in columns_to_drop if col in tcga_surv.columns])
            tcga_surv.rename(columns={'PFI': 'Status', 'PFI.time': 'Time'}, inplace=True)

        elif survival_type == "Disease-free Interval":
            tcga_clinics.drop(columns=["OS", "OS.time", 'DSS', 'DSS.time', 'PFI', 'PFI.time'], inplace=True)
            tcga_clinics.rename(columns={'DFI': 'Status', 'DFI.time': 'Time'}, inplace=True)
            columns_to_drop = ["OS", "OS.time", 'DSS', 'DSS.time', 'PFI', 'PFI.time']
            tcga_surv = tcga_surv.drop(columns=[col for col in columns_to_drop if col in tcga_surv.columns])
            tcga_surv.rename(columns={'DFI': 'Time', 'DFI.time': 'Status'}, inplace=True)
        else:
            tcga_clinics.drop(columns=["DSS", "DSS.time", 'DFI', 'DFI.time', 'PFI', 'PFI.time'], inplace=True)
            columns_to_drop = ["DSS", "DSS.time", 'DFI', 'DFI.time', 'PFI', 'PFI.time']
            tcga_surv = tcga_surv.drop(columns=[col for col in columns_to_drop if col in tcga_surv.columns])
            tcga_surv.rename(columns={'OS': 'Status', 'OS.time': 'Time'}, inplace=True)
            tcga_clinics.rename(columns={'OS': 'Status', 'OS.time': 'Time'}, inplace=True)

        tcga_surv.dropna(axis=0, how='any', inplace=True)

        return tcga_surv, tcga_clinics
    except Exception as e:
        st.error(f"Error in TCGA survival selector: {e}")
        return None, None


def transform_value(value, suffix, string_numbers):
    """
    Transform the value by appending the suffix if the value is in the predefined list.

    Parameters:
    - value (str): The value to be transformed.
    - suffix (str): The suffix to be appended.
    - string_numbers (list): List of predefined string numbers to check against.

    Returns:
    - transformed_value (str): Transformed value with the suffix appended if applicable.
    """
    try:
        # Check if the value is in string_numbers list, if yes, append the suffix
        if str(value) in string_numbers:
            return f"{value}_{suffix}"
        else:
            return value
    except Exception as e:
        st.error(f"An error occurred while transforming the value: {e}")
        return value
