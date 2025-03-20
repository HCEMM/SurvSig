import streamlit as st
import pandas as pd
import Scripts.descriptions as dsc


def non_unique_values(selected_subtypes):
    # Identify non-unique values across columns
    value_columns = {}
    for column in selected_subtypes.columns:
        for value in selected_subtypes[column].unique():
            if value not in value_columns:
                value_columns[value] = []
            value_columns[value].append(column)

    non_unique_values = {value for value, columns in value_columns.items() if len(columns) > 1}

    # Rename non-unique values across columns
    for idx, column in enumerate(selected_subtypes.columns, start=1):
        # Get unique values in the column
        unique_values = selected_subtypes[column].unique()
        # Create a mapping for non-unique values
        value_mapping = {
            value: f"{value}_{idx}"
            for value in unique_values
            if value in non_unique_values and pd.notna(value)
        }

        # Apply renaming within the column
        selected_subtypes[column] = selected_subtypes[column].apply(
            lambda x: value_mapping[x] if x in value_mapping else x
        )

    return selected_subtypes


def create_multiselect(df, column, i):
    """
    Create a multiselect widget for a given column in the DataFrame.

    Args:
        df (pd.DataFrame): The input DataFrame.
        column (str): The column name for which to create the multiselect widget.
        i (int): Index for unique key generation.

    Returns:
        list: List of selected values from the multiselect widget.
    """
    value_counts = df[column].value_counts(dropna=False).reset_index()
    value_counts.columns = [column, 'Count']
    value_counts_df = value_counts.set_index(column)
    unique_values = df[column].unique()
    df_multitype = pd.DataFrame(unique_values, columns=[column])
    include_list = [False for _ in range(len(df_multitype))]
    df_multitype.set_index(column, inplace=True)
    df_multitype['Count'] = value_counts_df["Count"].astype('str')
    df_multitype["Include"] = include_list
    df_multitype.reset_index(inplace=True)

    with st.container():
        st.subheader(f":blue[{str(column)}]")
        subtype_select = st.data_editor(df_multitype,
                                        column_config={"Include": st.column_config.CheckboxColumn(
                                            "Select Subtypes",
                                            help="Select Subtypes",
                                            default=True,
                                            width='small')},
                                        hide_index=True,
                                        use_container_width=True,
                                        key="subt_de" + str(i),
                                        height=200)

    st.header(" ", divider="blue")
    subtype_select = subtype_select[(subtype_select == True).any(axis=1)]
    subtype_select = subtype_select[column].values.tolist()

    return subtype_select


def subtype_select(selected_genes, dataset_anno, selected_genes_anno, anno_file_name, dataset_option, clinic_df):
    """
    Select subtypes based on user input and filter the data accordingly.

    Args:
        selected_genes (pd.DataFrame): DataFrame of selected genes.
        dataset_anno (pd.DataFrame): DataFrame of dataset annotations.
        selected_genes_anno (pd.DataFrame): DataFrame of selected genes annotations.
        anno_file_name (str): Path to save the filtered annotation file.
        dataset_option (str): Dataset option selected by the user.
        clinic_df (pd.DataFrame): DataFrame of clinical data.

    Returns:
        tuple: Filtered selected genes DataFrame, list of names, and merged DataFrame.
    """
    columns_of_interest = None

    dataset_anno.index.rename("Samples", inplace=True)
    clinic_df.index.rename("Samples", inplace=True)

    if dataset_option == 'Anish-SCLC':
        columns_of_interest = ["NAPY", "NE", 'Sex', 'Race', 'Smoke', 'Tumor_Stage', 'Disease', 'EMT']
        clinic_df.rename(columns={"SCLC staging/initial diagnosis ": "Tumor_Stage",
                                  "Smoking history": "Smoke"}, inplace=True)
    elif dataset_option == 'George-SCLC':
        columns_of_interest = ["NAPY", "NE", 'Sex', 'Smoker', 'UICC', 'Ethnicity Category', 'EMT']
    elif dataset_option == 'Jiang-SCLC':
        columns_of_interest = ["NAPY", "NE", 'Gender', 'Smoke', 'UICC stage', 'EMT', 'Chemotherapy']
    elif dataset_option == 'George-LCNEC':
        columns_of_interest = ["NAPY", "NE", 'Gender', 'Smoke', 'UICC Stage', 'EMT', 'Classification',
                               'Molecular_Subtypes']
    elif dataset_option == 'Fernandez-Carcinoid':
        columns_of_interest = ["NAPY", "NE", 'Gender', 'Smoke', 'UICC Stage', 'EMT', 'cluster_LNET', 'cluster_LNEN',
                               'Histopathology']
    elif dataset_option == 'Liu-SCLC':
        columns_of_interest = ["NAPY", "NE", 'Gender', 'Smoke', 'TNM_Stage', 'EMT', 'Tumor_site', 'Histologic_type']
    elif dataset_option == 'Alcala-Carcinoid':
        columns_of_interest = ["NAPY", "NE", 'Gender', 'Smoke', 'EMT', 'Histopathology_simplified',
                               'Molecular_clusters', 'Histopathology', 'Merged_Mol_clusters']
    elif dataset_option == 'Rousseaux-Mixed(non-NE)':
        columns_of_interest = ["NAPY", "NE", 'Gender', 'EMT', 'Histology', 'T_Stage', 'N', 'M', 'TNM', 'SurvSig_Mol.Subtype']
    elif dataset_option == 'TCGA':
        pass

    if dataset_option in ['Jiang-SCLC', 'George-LCNEC', 'Fernandez-Carcinoid', 'Liu-SCLC', 'Alcala-Carcinoid',
                          'Rousseaux-Mixed(non-NE)']:
        merged_df = clinic_df
    else:
        merged_df = dataset_anno.join(clinic_df)

    df_multi_subtype = []

    and_or_container = st.container(key="subtype_and_or")
    anno_and_or = and_or_container.radio(':blue[You Can Choose] :red["AND"] :blue[or] :red["OR"] :blue[Logic for Filtering ]'
                             ' :blue[Subtypes]', options=("OR", "AND"), horizontal=True, index=0,
                             key='or_and', help=dsc.sub_text())
    st.header(" ", divider="blue")

    with st.container():
        cols = st.columns(4)  # Create 3 columns
        for i, column in enumerate(columns_of_interest):
            with cols[i % 4]:
                subtype_select = create_multiselect(merged_df, column, i)
                df_multi_subtype += subtype_select

    filtered_selected_genes = merged_df[merged_df.isin(df_multi_subtype)]
    filtered_selected_genes.dropna(how='all', inplace=True)

    if len(df_multi_subtype) == 0:
        filtered_selected_genes = merged_df
    else:
        if anno_and_or == 'OR':
            conditions = [merged_df[column].isin(df_multi_subtype) for column in columns_of_interest]
            filtered_selected_genes = merged_df[pd.concat(conditions, axis=1).any(axis=1)]
        elif anno_and_or == 'AND':
            filtered_selected_genes = merged_df[
                merged_df.isin(df_multi_subtype).sum(axis=1) == len(df_multi_subtype)]

    filtered_selected_genes.to_csv(anno_file_name)
    name_list = filtered_selected_genes.index
    selected_genes_filtered = selected_genes_anno[selected_genes_anno.index.isin(name_list)]
    selected_genes_filtered = selected_genes_filtered.transpose()

    st.info(f":blue[Number of samples after filtering: ] "
            f" :red[{str(filtered_selected_genes.shape[0])}]  :blue[of ] :red[{str(selected_genes.shape[1])}]")

    return selected_genes_filtered, name_list, merged_df


def subtypes_selector_tcga(selected_genes, tcga_clinic, tcga_subtypes, anno_path, key):
    """
    Select subtypes for TCGA data based on user input and filter the data accordingly.

    Args:
        selected_genes (pd.DataFrame): DataFrame of selected genes.
        tcga_clinic (pd.DataFrame): DataFrame of TCGA clinical data.
        tcga_subtypes (pd.DataFrame): DataFrame of TCGA subtypes data.
        anno_path (str): Path to save the filtered annotation file.
        key (str): Suffix for unique key generation.

    Returns:
        tuple: Filtered selected genes DataFrame and filtered selected subtypes DataFrame.
    """

    selected_subtypes = None
    selected_genes = selected_genes.transpose()

    selected_subtypes1 = tcga_clinic[tcga_clinic["bcr_patient_barcode"].isin(selected_genes.index)].astype(
        'object').copy()
    selected_subtypes1.set_index("bcr_patient_barcode", inplace=True)
    selected_subtypes1.dropna(axis=1, inplace=True, how="all")

    try:
        selected_subtypes1 = selected_subtypes1[["gender", "race", "ajcc_pathologic_tumor_stage", "clinical_stage",
                                                 "histological_type", "histological_grade"]]
    except:
        selected_subtypes1 = selected_subtypes1[["gender", "race", "ajcc_pathologic_tumor_stage", "clinical_stage",
                                                 "histological_grade"]]

    selected_subtypes1 = selected_subtypes1.rename(columns={"gender": "Gender", "race": "Race",
                                                            "ajcc_pathologic_tumor_stage": "AJCC Stage",
                                                            "clinical_stage": "Clinical Stage",
                                                            "histological_type": "Hist. Type",
                                                            "histological_grade": "Hist. Grade"})
    selected_subtypes2 = tcga_subtypes[tcga_subtypes["pan.samplesID"].isin(selected_genes.index)].astype(
        'object').copy()

    selected_subtypes2.set_index("pan.samplesID", inplace=True)
    selected_subtypes2.dropna(axis=1, inplace=True, how="all")
    selected_subtypes2 = selected_subtypes2[selected_subtypes2.columns[1:]]
    selected_subtypes2.columns = selected_subtypes2.columns.str.replace("_", "-")
    selected_subtypes2.columns = selected_subtypes2.columns.str.replace("Subtype", "ST")

    selected_subtypes = selected_subtypes1.join(selected_subtypes2, how='left')
    selected_subtypes.fillna('nan', axis=1, inplace=True)

    selected_subtypes = non_unique_values(selected_subtypes)

    columns_of_interest = selected_subtypes.columns.tolist()[0:]
    columns_of_interest = [col for col in columns_of_interest if selected_subtypes[col].nunique() > 1]

    df_multi_subtype = []

    and_or_container = st.container(key="subtype_and_or_tcga")
    anno_and_or = and_or_container.radio(':blue[You Can Choose] :red["AND"] :blue[or] :red["OR"] :blue[Logic for Filtering ]'
                             ' :blue[Subtypes]', options=("OR", "AND"), horizontal=True, index=0,
                             key='or_and_tcga', help=dsc.sub_text())

    with st.container():
        cols = st.columns(4)  # Create 3 columns
        for i, column in enumerate(columns_of_interest):
            with cols[i % 4]:
                subtype_select = create_multiselect(selected_subtypes, column, i)
                df_multi_subtype += subtype_select

    filtered_selected_genes = selected_subtypes[selected_subtypes.isin(df_multi_subtype)]
    filtered_selected_genes.dropna(how='all', inplace=True)

    if len(df_multi_subtype) == 0:
        filtered_selected_genes = selected_subtypes
    else:
        if anno_and_or == 'OR':
            conditions = [selected_subtypes[column].isin(df_multi_subtype) for column in columns_of_interest]
            filtered_selected_genes = selected_subtypes[pd.concat(conditions, axis=1).any(axis=1)]
        elif anno_and_or == 'AND':
            filtered_selected_genes = selected_subtypes[selected_subtypes.isin(df_multi_subtype).sum(axis=1) ==
                                                        len(df_multi_subtype)]

    filtered_selected_genes.to_csv(anno_path)
    name_list = filtered_selected_genes.index
    selected_genes_filtered = selected_genes[selected_genes.index.isin(name_list)]
    selected_genes_filtered = selected_genes_filtered.transpose()

    st.info(f":blue[Number of samples after filtering: ] "
            f" :red[{str(filtered_selected_genes.shape[0])}]  :blue[of ] :red[{str(selected_genes.shape[0])}]")


    return selected_genes_filtered, filtered_selected_genes


def subtypes_selector_tcga_sg(selected_genes, tcga_clinic, tcga_subtypes, anno_path, sg_exp, key):
    """
    Select subtypes for TCGA data based on user input and filter the data accordingly, using a specific expander.

    Args:
        selected_genes (pd.DataFrame): DataFrame of selected genes.
        tcga_clinic (pd.DataFrame): DataFrame of TCGA clinical data.
        tcga_subtypes (pd.DataFrame): DataFrame of TCGA subtypes data.
        anno_path (str): Path to save the filtered annotation file.
        sg_exp (st.expander): Streamlit expander for displaying additional options.
        key (str): Suffix for unique key generation.

    Returns:
        tuple: Filtered selected genes DataFrame and filtered selected subtypes DataFrame.
    """
    selected_subtypes = None
    selected_genes = selected_genes.transpose()

    selected_subtypes1 = tcga_clinic[tcga_clinic["bcr_patient_barcode"].isin(selected_genes.index)].astype(
        'object').copy()
    selected_subtypes1.set_index("bcr_patient_barcode", inplace=True)
    selected_subtypes1.dropna(axis=1, inplace=True, how="all")

    try:
        selected_subtypes1 = selected_subtypes1[["gender", "race", "ajcc_pathologic_tumor_stage", "clinical_stage",
                                                 "histological_type", "histological_grade"]]
    except:
        selected_subtypes1 = selected_subtypes1[["gender", "race", "ajcc_pathologic_tumor_stage", "clinical_stage",
                                                 "histological_grade"]]

    selected_subtypes1 = selected_subtypes1.rename(columns={"gender": "Gender", "race": "Race",
                                                            "ajcc_pathologic_tumor_stage": "AJCC Stage",
                                                            "clinical_stage": "Clinical Stage",
                                                            "histological_type": "Hist. Type",
                                                            "histological_grade": "Hist. Grade"})

    selected_subtypes2 = tcga_subtypes[tcga_subtypes["pan.samplesID"].isin(selected_genes.index)].astype(
        'object').copy()
    selected_subtypes2.set_index("pan.samplesID", inplace=True)
    selected_subtypes2.dropna(axis=1, inplace=True, how="all")
    selected_subtypes2 = selected_subtypes2[selected_subtypes2.columns[1:]]
    selected_subtypes2.columns = selected_subtypes2.columns.str.replace("_", "-")
    selected_subtypes2.columns = selected_subtypes2.columns.str.replace("Subtype", "ST")

    selected_subtypes = selected_subtypes1.join(selected_subtypes2, how='left')
    selected_subtypes.fillna('nan', axis=1, inplace=True)
    selected_subtypes = non_unique_values(selected_subtypes)
    columns_of_interest = selected_subtypes.columns.tolist()[0:]
    columns_of_interest = [col for col in columns_of_interest if selected_subtypes[col].nunique() > 1]
    df_multi_subtype = []

    and_or_container = sg_exp.container(key="subtype_and_or_tcga_sg"+key)
    anno_and_or = and_or_container.radio(
        ':blue[You Can Choose] :red["AND"] :blue[or] :red["OR"] :blue[Logic for Filtering ]'
        ' :blue[Subtypes]', options=("OR", "AND"), horizontal=True, index=0,
        key='or_and_tcga_sg' + key, help=dsc.sub_text())

    with sg_exp.container():
        cols = sg_exp.columns(3)  # Create 3 columns
        for i, column in enumerate(columns_of_interest):
            with cols[i % 3]:
                subtype_select = create_multiselect_sg(selected_subtypes, column, i, key, sg_exp)
                df_multi_subtype += subtype_select

    filtered_selected_genes = selected_subtypes[selected_subtypes.isin(df_multi_subtype)]
    filtered_selected_genes.dropna(how='all', inplace=True)

    if len(df_multi_subtype) == 0:
        filtered_selected_genes = selected_subtypes
    else:
        if anno_and_or == 'OR':
            conditions = [selected_subtypes[column].isin(df_multi_subtype) for column in columns_of_interest]
            filtered_selected_genes = selected_subtypes[pd.concat(conditions, axis=1).any(axis=1)]
        elif anno_and_or == 'AND':
            filtered_selected_genes = selected_subtypes[selected_subtypes.isin(df_multi_subtype).sum(axis=1) ==
                                                        len(df_multi_subtype)]

    # Save the reindex annotation dataframe as a CSV file
    name_list = filtered_selected_genes.index
    selected_genes_filtered = selected_genes[selected_genes.index.isin(name_list)]
    selected_genes_filtered = selected_genes_filtered.transpose()

    sg_exp.info(f":blue[Number of samples after filtering: ]"
                f" :red[{str(filtered_selected_genes.shape[0])}]  :blue[of ] :red[{str(selected_genes.shape[0])}]")

    return selected_genes_filtered, filtered_selected_genes


def create_multiselect_sg(df, column, i, key, sg_exp):
    """
    Create a multiselect widget for a given column in the DataFrame using a specific expander.

    Args:
        df (pd.DataFrame): The input DataFrame.
        column (str): The column name for which to create the multiselect widget.
        i (int): Index for unique key generation.
        key (str): Suffix for unique key generation.
        sg_exp (st.expander): Streamlit expander for displaying additional options.

    Returns:
        list: List of selected values from the multiselect widget.
    """
    value_counts = df[column].value_counts(dropna=False).reset_index()
    value_counts.columns = [column, 'Count']
    value_counts_df = value_counts.set_index(column)
    unique_values = df[column].unique()
    df_multitype = pd.DataFrame(unique_values, columns=[column])
    include_list = [False for _ in range(len(df_multitype))]
    df_multitype.set_index(column, inplace=True)
    df_multitype['Count'] = value_counts_df["Count"].astype('str')
    df_multitype["Include"] = include_list
    df_multitype.reset_index(inplace=True)

    with sg_exp.container():
        sg_exp.subheader(f":blue[{str(column)}]")
        subtype_select = sg_exp.data_editor(df_multitype,
                                            column_config={"Include": st.column_config.CheckboxColumn(
                                                "Select Subtypes",
                                                help="Select Subtypes",
                                                default=True,
                                                width='small'),
                                            },
                                            hide_index=True,
                                            use_container_width=True,
                                            key="subt_de" + str(i) + key,
                                            height=200)

    sg_exp.header(" ", divider="blue")
    subtype_select = subtype_select[(subtype_select == True).any(axis=1)]
    subtype_select = subtype_select[column].values.tolist()

    return subtype_select


def subtype_select_sg(selected_genes, dataset_anno, selected_genes_anno, dataset_option, clinic_df, sg_exp, key):
    """
    Select subtypes based on user input and filter the data accordingly, using a specific expander.

    Args:
        selected_genes (pd.DataFrame): DataFrame of selected genes.
        dataset_anno (pd.DataFrame): DataFrame of dataset annotations.
        selected_genes_anno (pd.DataFrame): DataFrame of selected genes annotations.
        dataset_option (str): Dataset option selected by the user.
        clinic_df (pd.DataFrame): DataFrame of clinical data.
        sg_exp (st.expander): Streamlit expander for displaying additional options.
        key (str): Suffix for unique key generation.

    Returns:
        tuple: Filtered selected genes DataFrame, list of names, and merged DataFrame.
    """
    columns_of_interest = None

    dataset_anno.index.rename("Samples", inplace=True)
    clinic_df.index.rename("Samples", inplace=True)

    if dataset_option == 'Anish-SCLC':
        columns_of_interest = ["NAPY", "NE", 'Sex', 'Race', 'Smoke', 'Tumor_Stage', 'EMT']
        clinic_df.rename(columns={"SCLC staging/initial diagnosis ": "Tumor_Stage",
                                  "Smoking history": "Smoke"}, inplace=True)
    elif dataset_option == 'George-SCLC':
        columns_of_interest = ["NAPY", "NE", 'Sex', 'Smoker', 'UICC', 'Ethnicity Category', 'EMT']
    elif dataset_option == 'Jiang-SCLC':
        columns_of_interest = ["NAPY", "NE", 'Gender', 'Smoke', 'UICC stage', 'EMT', 'Chemotherapy']
    elif dataset_option == 'George-LCNEC':
        columns_of_interest = ["NAPY", "NE", 'Gender', 'Smoke', 'UICC Stage', 'EMT', 'Classification',
                               'Molecular_Subtypes']
    elif dataset_option == 'Fernandez-Carcinoid':
        columns_of_interest = ["NAPY", "NE", 'Gender', 'Smoke', 'UICC Stage', 'EMT', 'cluster_LNET', 'cluster_LNEN',
                               'Histopathology']
    elif dataset_option == 'Liu-SCLC':
        columns_of_interest = ["NAPY", "NE", 'Gender', 'Smoke', 'TNM_Stage', 'EMT', 'Tumor_site', 'Histologic_type']
    elif dataset_option == 'Alcala-Carcinoid':
        columns_of_interest = ["NAPY", "NE", 'Gender', 'Smoke', 'EMT', 'Histopathology_simplified',
                               'Molecular_clusters', 'Histopathology', 'Merged_Mol_clusters']
    elif dataset_option == 'Rousseaux-Mixed(non-NE)':
        columns_of_interest = ["NAPY", "NE", 'Gender', 'EMT', 'Histology', 'T_Stage', 'N', 'M', 'TNM', 'SurvSig_Mol.Subtype']
    elif dataset_option == 'TCGA':
        pass

    if dataset_option in ['Jiang-SCLC', 'George-LCNEC', 'Fernandez-Carcinoid', 'Liu-SCLC', 'Alcala-Carcinoid',
                          'Rousseaux-Mixed(non-NE)']:
        merged_df = clinic_df
    else:
        merged_df = dataset_anno.join(clinic_df)

    df_multi_subtype = []

    and_or_container = sg_exp.container(key="subtype_and_or_sg" + key)
    anno_and_or = and_or_container.radio(
        ':blue[You Can Choose] :red["AND"] :blue[or] :red["OR"] :blue[Logic for Filtering ]'
        ' :blue[Subtypes]', options=("OR", "AND"), horizontal=True, index=0,
        key='and_or_sg' + key, help=dsc.sub_text())

    with sg_exp.container():
        cols = st.columns(4)  # Create 3 columns
        for i, column in enumerate(columns_of_interest):
            with cols[i % 4]:
                subtype_select = create_multiselect_sg(merged_df, column, i, key, sg_exp)
                df_multi_subtype += subtype_select

    filtered_selected_genes = merged_df[merged_df.isin(df_multi_subtype)]
    filtered_selected_genes.dropna(how='all', inplace=True)

    if len(df_multi_subtype) == 0:
        filtered_selected_genes = merged_df
    else:
        if anno_and_or == 'OR':
            conditions = [merged_df[column].isin(df_multi_subtype) for column in columns_of_interest]
            filtered_selected_genes = merged_df[pd.concat(conditions, axis=1).any(axis=1)]
        elif anno_and_or == 'AND':
            filtered_selected_genes = merged_df[merged_df.isin(df_multi_subtype).sum(axis=1) == len(df_multi_subtype)]

    name_list = filtered_selected_genes.index
    selected_genes_filtered = selected_genes_anno[selected_genes_anno.index.isin(name_list)]
    selected_genes_filtered = selected_genes_filtered.transpose()

    sg_exp.info(f":blue[Number of samples after filtering: ] "
                f" :red[{str(selected_genes_filtered.shape[1])} ]  :blue[of ] :red[{str(selected_genes.shape[0])}]")

    return selected_genes_filtered, name_list, merged_df


def subtype_select_gf(selected_genes, dataset_anno, selected_genes_anno, dataset_option, clinic_df, sg_exp, key):

    """
    Select subtypes based on user input and filter the data accordingly, using a specific expander.

    Args:
        selected_genes (pd.DataFrame): DataFrame of selected genes.
        dataset_anno (pd.DataFrame): DataFrame of dataset annotations.
        selected_genes_anno (pd.DataFrame): DataFrame of selected genes annotations.
        dataset_option (str): Dataset option selected by the user.
        clinic_df (pd.DataFrame): DataFrame of clinical data.
        sg_exp (st.expander): Streamlit expander for displaying additional options.
        key (str): Suffix for unique key generation.

    Returns:
        tuple: Filtered selected genes DataFrame, list of names, and merged DataFrame.
    """
    columns_of_interest = None

    dataset_anno.index.rename("Samples", inplace=True)
    clinic_df.index.rename("Samples", inplace=True)

    if dataset_option == 'Anish-SCLC':
        columns_of_interest = ["NAPY", "NE", 'Sex', 'Race', 'Smoke', 'Tumor_Stage', 'EMT']
        clinic_df.rename(columns={"SCLC staging/initial diagnosis ": "Tumor_Stage",
                                  "Smoking history": "Smoke"}, inplace=True)
    elif dataset_option == 'George-SCLC':
        columns_of_interest = ["NAPY", "NE", 'Sex', 'Smoker', 'UICC', 'Ethnicity Category', 'EMT']
    elif dataset_option == 'Jiang-SCLC':
        columns_of_interest = ["NAPY", "NE", 'Gender', 'Smoke', 'UICC stage', 'EMT', 'Chemotherapy']
    elif dataset_option == 'George-LCNEC':
        columns_of_interest = ["NAPY", "NE", 'Gender', 'Smoke', 'UICC Stage', 'EMT', 'Classification',
                               'Molecular_Subtypes']
    elif dataset_option == 'Fernandez-Carcinoid':
        columns_of_interest = ["NAPY", "NE", 'Gender', 'Smoke', 'UICC Stage', 'EMT', 'cluster_LNET', 'cluster_LNEN',
                               'Histopathology']
    elif dataset_option == 'Liu-SCLC':
        columns_of_interest = ["NAPY", "NE", 'Gender', 'Smoke', 'TNM_Stage', 'EMT', 'Tumor_site', 'Histologic_type']
    elif dataset_option == 'Alcala-Carcinoid':
        columns_of_interest = ["NAPY", "NE", 'Gender', 'Smoke', 'EMT', 'Histopathology_simplified',
                               'Molecular_clusters', 'Histopathology', 'Merged_Mol_clusters']
    elif dataset_option == 'Rousseaux-Mixed(non-NE)':
        columns_of_interest = ["NAPY", "NE", 'Gender', 'EMT', 'Histology', 'T_Stage', 'N', 'M', 'TNM', 'SurvSig_Mol.Subtype']
    elif dataset_option == 'TCGA':
        pass

    if dataset_option in ['Jiang-SCLC', 'George-LCNEC', 'Fernandez-Carcinoid', 'Liu-SCLC', 'Alcala-Carcinoid',
                          'Rousseaux-Mixed(non-NE)']:
        merged_df = clinic_df
    else:
        merged_df = dataset_anno.join(clinic_df)
    df_multi_subtype = []

    and_or_container = sg_exp.container(key="subtype_and_or_gf")
    anno_and_or = and_or_container.radio(
        ':blue[You Can Choose] :red["AND"] :blue[or] :red["OR"] :blue[Logic for Filtering ]'
        ' :blue[Subtypes]', options=("OR", "AND"), horizontal=True, index=0,
        key='and_or_gf', help=dsc.sub_text())

    with sg_exp.container():
        cols = st.columns(4)  # Create 3 columns
        for i, column in enumerate(columns_of_interest):
            with cols[i % 4]:
                subtype_select = create_multiselect_sg(merged_df, column, i, key, sg_exp)
                df_multi_subtype += subtype_select

    filtered_selected_genes = merged_df[merged_df.isin(df_multi_subtype)]
    filtered_selected_genes.dropna(how='all', inplace=True)

    if len(df_multi_subtype) == 0:
        filtered_selected_genes = merged_df
    else:
        if anno_and_or == 'OR':
            conditions = [merged_df[column].isin(df_multi_subtype) for column in columns_of_interest]
            filtered_selected_genes = merged_df[pd.concat(conditions, axis=1).any(axis=1)]
        elif anno_and_or == 'AND':
            filtered_selected_genes = merged_df[merged_df.isin(df_multi_subtype).sum(axis=1) == len(df_multi_subtype)]

    name_list = filtered_selected_genes.index
    selected_genes_filtered = selected_genes_anno[selected_genes_anno.index.isin(name_list)]
    selected_genes_filtered = selected_genes_filtered.transpose()

    sg_exp.info(f":blue[Number of samples after filtering: ] "
                f" :red[{str(selected_genes_filtered.shape[1])} ]  :blue[of ] :red[{str(selected_genes.shape[0])}]")

    return selected_genes_filtered, merged_df


def subtypes_selector_tcga_gf(selected_genes, tcga_clinic, tcga_subtypes, sg_exp, key):
    """
    Select subtypes for TCGA data based on user input and filter the data accordingly, using a specific expander.

    Args:
        selected_genes (pd.DataFrame): DataFrame of selected genes.
        tcga_clinic (pd.DataFrame): DataFrame of TCGA clinical data.
        tcga_subtypes (pd.DataFrame): DataFrame of TCGA subtypes data.
        sg_exp (st.expander): Streamlit expander for displaying additional options.
        key (str): Suffix for unique key generation.

    Returns:
        tuple: Filtered selected genes DataFrame and filtered selected subtypes DataFrame.
    """
    selected_subtypes = None
    selected_genes = selected_genes.transpose()

    selected_subtypes1 = tcga_clinic[tcga_clinic["bcr_patient_barcode"].isin(selected_genes.index)].astype(
        'object').copy()
    selected_subtypes1.set_index("bcr_patient_barcode", inplace=True)
    selected_subtypes1.dropna(axis=1, inplace=True, how="all")

    try:
        selected_subtypes1 = selected_subtypes1[["gender", "race", "ajcc_pathologic_tumor_stage", "clinical_stage",
                                                 "histological_type", "histological_grade"]]
    except:
        selected_subtypes1 = selected_subtypes1[["gender", "race", "ajcc_pathologic_tumor_stage", "clinical_stage",
                                                 "histological_grade"]]

    selected_subtypes1 = selected_subtypes1.rename(columns={"gender": "Gender", "race": "Race",
                                                            "ajcc_pathologic_tumor_stage": "AJCC Stage",
                                                            "clinical_stage": "Clinical Stage",
                                                            "histological_type": "Hist. Type",
                                                            "histological_grade": "Hist. Grade"})

    selected_subtypes2 = tcga_subtypes[tcga_subtypes["pan.samplesID"].isin(selected_genes.index)].astype(
        'object').copy()
    selected_subtypes2.set_index("pan.samplesID", inplace=True)
    selected_subtypes2.dropna(axis=1, inplace=True, how="all")
    selected_subtypes2 = selected_subtypes2[selected_subtypes2.columns[1:]]
    selected_subtypes2.columns = selected_subtypes2.columns.str.replace("_", "-")
    selected_subtypes2.columns = selected_subtypes2.columns.str.replace("Subtype", "ST")

    selected_subtypes = selected_subtypes1.join(selected_subtypes2, how='left')
    selected_subtypes.fillna('nan', axis=1, inplace=True)
    columns_of_interest = selected_subtypes.columns.tolist()[0:]
    columns_of_interest = [col for col in columns_of_interest if selected_subtypes[col].nunique() > 1]
    df_multi_subtype = []

    and_or_container = sg_exp.container(key="subtype_and_or_gf_tcga")
    anno_and_or = and_or_container.radio(
        ':blue[You Can Choose] :red["AND"] :blue[or] :red["OR"] :blue[Logic for Filtering ]'
        ' :blue[Subtypes]', options=("OR", "AND"), horizontal=True, index=0,
        key='and_or_gf_tcga', help=dsc.sub_text())

    with sg_exp.container():
        cols = sg_exp.columns(3)  # Create 3 columns
        for i, column in enumerate(columns_of_interest):
            with cols[i % 3]:
                subtype_select = create_multiselect_sg(selected_subtypes, column, i, key, sg_exp)
                df_multi_subtype += subtype_select

    filtered_selected_genes = selected_subtypes[selected_subtypes.isin(df_multi_subtype)]
    filtered_selected_genes.dropna(how='all', inplace=True)

    if len(df_multi_subtype) == 0:
        filtered_selected_genes = selected_subtypes
    else:
        if anno_and_or == 'OR':
            conditions = [selected_subtypes[column].isin(df_multi_subtype) for column in columns_of_interest]
            filtered_selected_genes = selected_subtypes[pd.concat(conditions, axis=1).any(axis=1)]
        elif anno_and_or == 'AND':
            filtered_selected_genes = selected_subtypes[selected_subtypes.isin(df_multi_subtype).sum(axis=1) ==
                                                        len(df_multi_subtype)]

    # Save the reindex annotation dataframe as a CSV file
    name_list = filtered_selected_genes.index
    selected_genes_filtered = selected_genes[selected_genes.index.isin(name_list)]
    selected_genes_filtered = selected_genes_filtered.transpose()

    sg_exp.info(f":blue[Number of samples after filtering: ]"
                f" :red[{str(filtered_selected_genes.shape[0])}]  :blue[of ] :red[{str(selected_genes.shape[0])}]")

    return selected_genes_filtered, filtered_selected_genes
