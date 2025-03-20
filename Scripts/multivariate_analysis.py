import pandas as pd
import streamlit as st
from PIL import Image
from Scripts.rpy2_heatmap_plots import multi_var
from Scripts.r_code_settings import file_copy
from Scripts.data_loader import convert_string_to_float


@st.fragment
def mv_db(cluster, selected_genes, multi_db_name, user_id, dataset_option, multi_info_name, forest_pdf, forest_png,
          anno_file, clinic_df):
    """
    Multivariate database analysis function

    Args:
        cluster (DataFrame): Cluster information.
        selected_genes (DataFrame): Selected genes information.
        multi_db_name (str): Multivariate database name.
        user_id (str): User ID.
        dataset_option (str): Dataset option.
        multi_info_name (str): Multivariate information name.
        forest_pdf (str): Path to save forest plot PDF.
        forest_png (str): Path to save forest plot PNG.
        anno_file (DataFrame): Annotation file information.
        clinic_df (DataFrame): Clinical data frame.

    Returns:
        None
    """
    try:

        # Advanced settings for multivariate analysis
        multi_adv = st.expander(":blue[Advanced Settings]")
        multi_adv_col = multi_adv.columns(2)
        min_num = multi_adv_col[0].number_input(":blue[Select] :red[minimum number of groups]", min_value=1,
                                                max_value=10, value=3, key="min_num")

        mv_db1, cols_to_filter, available_columns = multi_datas(dataset_option, cluster, selected_genes, clinic_df,
                                                                anno_file)

        # Step 1: Identify all the values in the specified columns that do not meet the threshold.
        invalid_values = {}
        for col in available_columns:
            # Check if the column is in df and of object data type
            if col in mv_db1.columns and mv_db1[col].dtype == 'object':
                value_counts = mv_db1[col].value_counts()
                invalid_values[col] = value_counts[value_counts < min_num].index.tolist()
            if col in mv_db1.columns and mv_db1[col].dtype == 'int':
                value_counts = mv_db1[col].value_counts()
                invalid_values[col] = value_counts[value_counts < min_num].index.tolist()

        # Step 2: Set the values in the specified columns to NaN if they are in the invalid_values list.
        for col, values in invalid_values.items():
            mv_db1.loc[mv_db1[col].isin(values), col] = 'nan'

        # Filter columns with only one unique value and remove them from multi-select options
        columns_with_multiple_values = [
            col for col in mv_db1.columns
            if (mv_db1[col].nunique() >= 2) or
               (mv_db1[col].nunique() == 2 and not mv_db1[col].isin(['nan']).any())
        ]

        if dataset_option == "TCGA":
            selected_features_option = [
                col
                for col in mv_db1.columns
                if (mv_db1[col].dtype == 'object') and (mv_db1[col].nunique() > 1)
            ]
            selected_features_option.remove("Status")
            selected_features_option_default = selected_features_option[:2]
            selected_features_option_default.append("Cluster")

        else:
            selected_features_option_default = list(set(columns_with_multiple_values) & set(cols_to_filter))
            selected_features_option = list(set(columns_with_multiple_values) & set(available_columns))

        selected_features = multi_adv_col[1].multiselect(":blue[Select the] :red[Features]",
                                                         options=selected_features_option,
                                                         key="multi_adv_g", default=selected_features_option_default,
                                                         max_selections=5)

        # DataFrame to store selected values from each selectbox
        selectbox_df = pd.DataFrame(columns=['Selected References'])

        if len(selected_features) == 0:
            selected_features = ['Cluster']

        for feature in selected_features:
            # Filter out NaN values and convert to list
            unique_values = [value for value in mv_db1[feature].unique().tolist() if value != 'nan']

            # Getting the most common value for the feature while ensuring it's not NaN
            most_common_value = mv_db1[feature].mode().dropna().iloc[0] if not mv_db1[
                feature].mode().dropna().empty else None

            # If there is no most common value (all values are NaN), skip this feature
            if most_common_value is None:
                continue

            # Get the index of the most common value in the unique_values list, if it exists
            initial_index = unique_values.index(most_common_value) if most_common_value in unique_values else 0

            selected_val = multi_adv.selectbox(f":blue[Select] :red[reference] :blue[value for] :red[{feature}]",
                                               options=unique_values, key=f"select_mv_{feature}", index=initial_index)

            selectbox_df.loc[feature] = selected_val

        if len(selected_features) == mv_db1.shape[1]:
            multi_adv.info(":red[All] :blue[features are selected]")
        else:
            multi_adv.info(f":blue[Selected feature(s): ] :red[{str(selected_features)}]")

        if len(selected_features) != 0:
            selected_features.extend(["Status", "Time"])
            mv_db1 = mv_db1[selected_features]
        mv_db1.replace("NaN", "nan", inplace=True)

        mv_db1.to_csv(multi_db_name)

        multi_col1, multi_col2, multi_col3 = st.columns(3)
        try:
            multi_var(user_id, selectbox_df, multi_info_name, forest_pdf, forest_png, multi_db_name)
            image = Image.open(forest_png)
            multi_col2.subheader(":blue[Multivariate Plot]", divider="blue")
            multi_col2.image(image, use_container_width=True, output_format='PNG')
            multi_col2.toast(':red[Multivariate Analysis ] :blue[is Completed!]', icon='🎉')
        except Exception as e:
            multi_col2.subheader(":blue[Multivariate Plot]", divider="blue")
            multi_col2.image("style_items/error.svg", use_container_width=True)
            if mv_db1.dropna(how='all', inplace=True) is None:
                multi_col2.warning(":blue[No survival data is available!] :red[Please update your] :blue[subtype selection].")
            else:
                multi_col2.warning(":blue[You have ] :red[infinity ] :blue[value(s). Try to remove] :red[features] "
                                   ":blue[ or try to use other] :red[ clustering / dimensionality method]!")
        file_copy(user_id)
    except Exception as e:
           st.error(f"Error in multivariate database analysis: {e}")


def multi_datas(dataset_option, cluster, selected_genes, clinic_df, anno_file):
    mv_db1 = pd.DataFrame()
    # Dataset-specific processing
    if dataset_option == "George-SCLC":
        mv_db = pd.read_csv("source_data/sclc_ucologne_2015_clinical_data.tsv", delimiter="\t", index_col=1)
        mv_db = mv_db.astype(str)
        mv_db = mv_db[mv_db.index.isin(selected_genes.index)]
        mv_db1.index = mv_db.index
        cluster["Cluster"] = cluster["Cluster"].astype('str')

        mv_db1["Cluster"] = cluster["Cluster"]
        mv_db1['Time'] = mv_db["Overall Survival (Months)"]
        mv_db1['Status'] = mv_db["Overall Survival Status"]
        mv_db1['Chemotherapy'] = mv_db['Chemotherapy']
        mv_db1['EMT'] = anno_file['EMT']
        mv_db1['NE'] = anno_file['NE']
        mv_db1['NAPY'] = anno_file['NAPY']
        mv_db1['Ethnicity_Category'] = mv_db['Ethnicity Category']
        mv_db1['Sex'] = mv_db['Sex']
        mv_db1['UICC_stage'] = mv_db['UICC Tumor Stage']
        mv_db1 = mv_db1.replace(['0:LIVING', '1:DECEASED'], ['0', '1'])
        mv_db1 = mv_db1.replace(['IA', 'IB'], ['Ia', 'Ib'])
        mv_db1 = mv_db1.replace(['NaN'], ['nan'])
        mv_db1['Cluster'] = mv_db1['Cluster'].astype(str)
        cols_to_filter = ['UICC_stage', 'Cluster', 'NE', 'EMT', 'NAPY']
        available_columns = ['UICC_stage', 'Cluster', 'Chemotherapy', 'Sex', 'NE', 'EMT', 'NAPY', 'Ethnicity_Category']

    elif dataset_option == 'Anish-SCLC':
        mv_db = pd.read_excel("source_data/Supplementary_table_patient_clinical_characteristics_2021015.xlsx",
                              index_col=2)
        mv_db = mv_db.astype(str)
        mv_db = mv_db[mv_db.index.isin(selected_genes.index)]
        mv_db1.index = mv_db.index
        mv_db1.index.rename("Patient ID", inplace=True)
        cluster["Cluster"] = cluster["Cluster"].astype('str')
        mv_db1["Cluster"] = cluster["Cluster"]
        mv_db1['Time'] = mv_db["Overall survival since diagnosis"]
        mv_db1['Time'] = mv_db1['Time'].astype('float')
        mv_db1['Time'] = mv_db1['Time'].astype('int')
        mv_db1['Time'] = mv_db1['Time'].astype('str')
        mv_db1['Status'] = mv_db["Mortality status"]
        mv_db1 = mv_db1.replace(['Dead', 'Alive'], ['1', '0'])
        mv_db1['Immunotherapy'] = mv_db["Immunotherapy prior to biopsy"]
        mv_db1['Smoking history'] = mv_db["Smoking history"]
        mv_db1['EMT'] = anno_file['EMT']
        mv_db1['NE'] = anno_file['NE']
        mv_db1['NAPY'] = anno_file['NAPY']
        mv_db1['Sex'] = mv_db["Sex"]
        mv_db1['SCLC_staging'] = mv_db['SCLC staging/initial diagnosis ']
        mv_db1['Cluster'] = mv_db1['Cluster'].astype(str)
        cols_to_filter = ['Cluster', 'SCLC_staging', 'Immunotherapy', 'Sex']
        available_columns = ['Cluster', 'SCLC_staging', 'Immunotherapy', 'Sex', 'EMT', 'NAPY', 'NE']

    elif dataset_option == 'Jiang-SCLC':
        mv_db = pd.read_csv("source_data/Jiang_SCLC_annotation.csv", index_col=0)
        mv_db = mv_db.astype(str)
        mv_db = mv_db[mv_db.index.isin(selected_genes.index)]
        mv_db1.index = mv_db.index
        mv_db1.index.rename("Patient ID", inplace=True)
        cluster["Cluster"] = cluster["Cluster"].astype('str')
        mv_db1["Cluster"] = cluster["Cluster"]
        mv_db1['Time'] = mv_db["Months"]
        mv_db1['Time'] = mv_db1['Time'].astype('float')
        mv_db1.dropna(how='any', inplace=True)
        mv_db1['Time'] = mv_db1['Time'].astype('int')
        mv_db1['Time'] = mv_db1['Time'].astype('str')
        mv_db1['Status'] = mv_db["Survival Status"]
        mv_db1['Chemotherapy'] = mv_db["chemotherapy-exposed biopsy"]
        mv_db1['Smoke'] = mv_db["Smoke"]
        mv_db1['Gender'] = mv_db["Gender"]
        mv_db1['SCLC_staging'] = mv_db['UICC stage']
        mv_db1['EMT'] = anno_file['EMT']
        mv_db1['NE'] = anno_file['NE']
        mv_db1['Cluster'] = mv_db1['Cluster'].astype(str)
        cols_to_filter = ['Cluster', 'SCLC_staging', 'NE']
        available_columns = ['Cluster', 'SCLC_staging', 'Gender', 'EMT', 'NE']

    elif dataset_option == 'George-LCNEC':
        mv_db = pd.read_csv("source_data/LCNEC_George_clinical_anno.csv", index_col=0)
        mv_db.replace(
            ['[Not Available]', '[Unknown]', '[Not Evaluated]', 'nan', '[Not Applicable]', '[Discrepancy]', 'NaN'],
            value=None, inplace=True)
        mv_db.dropna(how='all', inplace=True, axis=1)
        mv_db.dropna(how='all', inplace=True)
        mv_db = mv_db.astype(str)
        mv_db = mv_db[mv_db.index.isin(selected_genes.index)]
        mv_db1.index = mv_db.index
        mv_db1.index.rename("Patient ID", inplace=True)
        cluster["Cluster"] = cluster["Cluster"].astype('str')
        mv_db1["Cluster"] = cluster["Cluster"]
        mv_db1['Time'] = mv_db["Survival_months"]
        mv_db1['Time'] = mv_db1['Time'].astype('float')
        mv_db1.dropna(how='any', inplace=True)
        mv_db1['Time'] = mv_db1['Time'].astype('int')
        mv_db1['Time'] = mv_db1['Time'].astype('str')
        mv_db1['Status'] = mv_db["Survival_censor"]
        mv_db1 = mv_db1.replace(['dead', 'alive'], ['1', '0'])
        mv_db1['Molecular_Subtypes'] = mv_db["Molecular_Subtypes"]
        mv_db1['Smoke'] = mv_db["Smoking_Hx"]
        mv_db1['Gender'] = mv_db["Sex"]
        mv_db1['SCLC_staging'] = mv_db['Stage_UICC']
        mv_db1['EMT'] = anno_file['EMT']
        mv_db1['NE'] = anno_file['NE']
        mv_db1['NAPY'] = anno_file['NAPY']
        mv_db1['Classification'] = mv_db['Classification']
        mv_db1['Cluster'] = mv_db1['Cluster'].astype(str)
        cols_to_filter = ['Cluster', 'SCLC_staging', 'Gender', 'Classification']
        available_columns = ['Cluster', 'SCLC_staging', 'Gender', 'Classification', 'EMT', 'NE', 'NAPY']

    elif dataset_option == 'Fernandez-Carcinoid':
        mv_db = pd.read_csv("source_data/carcinoid_clinical_anno.csv", index_col=0)
        mv_db.replace(
            ['[Not Available]', '[Unknown]', '[Not Evaluated]', 'nan', '[Not Applicable]', '[Discrepancy]', 'NaN'],
            value=None, inplace=True)
        mv_db.dropna(how='all', inplace=True, axis=1)
        mv_db.dropna(how='all', inplace=True)
        mv_db = mv_db.astype(str)
        mv_db = mv_db[mv_db.index.isin(selected_genes.index)]
        mv_db1.index = mv_db.index
        mv_db1.index.rename("Patient ID", inplace=True)
        cluster["Cluster"] = cluster["Cluster"].astype('str')
        mv_db1["Cluster"] = cluster["Cluster"]
        mv_db1['Time'] = mv_db["Survival_months"]
        mv_db1['Time'] = mv_db1['Time'].astype('float')
        mv_db1.dropna(how='any', inplace=True)
        mv_db1['Time'] = mv_db1['Time'].astype('int')
        mv_db1['Time'] = mv_db1['Time'].astype('str')
        mv_db1['Status'] = mv_db["Survival_censor"]
        mv_db1 = mv_db1.replace(['dead', 'alive'], ['1', '0'])
        mv_db1['Histopathology'] = mv_db["Histopathology"]
        mv_db1['cluster_LNET'] = mv_db["cluster_LNET"]
        mv_db1['cluster_LNEN'] = mv_db["cluster_LNEN"]
        mv_db1['Smoke'] = mv_db["Smoking_status"]
        mv_db1['Gender'] = mv_db["Sex"]
        mv_db1['EMT'] = anno_file['EMT']
        mv_db1['NE'] = anno_file['NE']
        mv_db1['NAPY'] = anno_file['NAPY']
        mv_db1['SCLC_staging'] = mv_db['Stage_UICC']
        mv_db1['Cluster'] = mv_db1['Cluster'].astype(str)
        cols_to_filter = ['Cluster', 'SCLC_staging', 'Gender', 'Histopathology', 'cluster_LNET']
        available_columns = ['Cluster', 'SCLC_staging', 'Gender', 'Histopathology', 'cluster_LNET', 'EMT', 'NE', 'NAPY']

    elif dataset_option == 'Liu-SCLC':
        mv_db = pd.read_csv("source_data/liu_clinical_anno.csv", index_col=0)
        mv_db.replace(
            ['[Not Available]', '[Unknown]', '[Not Evaluated]', 'nan', '[Not Applicable]', '[Discrepancy]', 'NaN',
             "unknown"], value=None, inplace=True)
        mv_db.dropna(how='all', inplace=True, axis=1)
        mv_db.dropna(how='all', inplace=True)
        mv_db.index = mv_db.index.str.upper()
        mv_db['Survial.(months)'] = mv_db['Survial.(months)'].str.split(r' \(').str[0]
        mv_db = mv_db.applymap(convert_string_to_float)
        mv_db = mv_db.astype(str)
        mv_db = mv_db[mv_db.index.isin(selected_genes.index)]
        mv_db1.index = mv_db.index
        mv_db1.index.rename("Patient ID", inplace=True)
        cluster["Cluster"] = cluster["Cluster"].astype('str')
        mv_db1["Cluster"] = cluster["Cluster"]
        mv_db1['Time'] = mv_db["Survial.(months)"]
        mv_db1['Time'] = mv_db1['Time'].astype('float')
        mv_db1.dropna(how='any', inplace=True)
        mv_db1['Time'] = mv_db1['Time'].astype('int')
        mv_db1['Time'] = mv_db1['Time'].astype('str')
        mv_db1['Status'] = mv_db["Status."]
        mv_db1 = mv_db1.replace(['dead', 'alive'], ['1', '0'])
        mv_db1['Smoke'] = mv_db["Smoking"]
        mv_db1['Tumor_site'] = mv_db["Tumor.site"]
        mv_db1['Histologic_type'] = mv_db["Histologic.type"]
        mv_db1['Gender'] = mv_db["Gender"]
        mv_db1['EMT'] = anno_file['EMT']
        mv_db1['NE'] = anno_file['NE']
        mv_db1['NAPY'] = anno_file['NAPY']
        mv_db1['Stage'] = mv_db['TNM.Stage']
        mv_db1['Cluster'] = mv_db1['Cluster'].astype(str)
        cols_to_filter = ['Cluster', 'Stage', 'Gender', 'Tumor_site', 'Smoke']
        available_columns = ['Cluster', 'Stage', 'Gender', 'Tumor_site', 'Smoke', 'EMT', 'NE', 'NAPY']

    elif dataset_option == 'Alcala-Carcinoid':
        mv_db = pd.read_csv("source_data/alcala_clinical_anno.csv", index_col=0)
        mv_db.replace(
            ['[Not Available]', '[Unknown]', '[Not Evaluated]', 'nan', '[Not Applicable]', '[Discrepancy]', 'NaN'],
            value=None, inplace=True)
        mv_db.dropna(how='all', inplace=True, axis=1)
        mv_db.dropna(how='all', inplace=True)
        mv_db = mv_db.astype(str)
        mv_db = mv_db[mv_db.index.isin(selected_genes.index)]
        mv_db1.index = mv_db.index
        mv_db1.index.rename("Patient ID", inplace=True)
        cluster["Cluster"] = cluster["Cluster"].astype('str')
        mv_db1["Cluster"] = cluster["Cluster"]
        mv_db1['Time'] = mv_db["Survival_months"]
        mv_db1['Time'] = mv_db1['Time'].astype('float')
        mv_db1.dropna(how='any', inplace=True)
        mv_db1['Time'] = mv_db1['Time'].astype('int')
        mv_db1['Time'] = mv_db1['Time'].astype('str')
        mv_db1['Status'] = mv_db["Survival_censor"]
        mv_db1 = mv_db1.replace(['dead', 'alive'], ['1', '0'])
        mv_db1['Histopathology'] = mv_db["Histopathology"]
        mv_db1['Histopathology_simplified'] = mv_db["Histopathology_simplified"]
        mv_db1['Molecular_clusters'] = mv_db["Molecular_clusters"]
        mv_db1['Smoke'] = mv_db["Smoking_status"]
        mv_db1['Gender'] = mv_db["Sex"]
        mv_db1['EMT'] = anno_file['EMT']
        mv_db1['NE'] = anno_file['NE']
        mv_db1['NAPY'] = anno_file['NAPY']
        mv_db1['Cluster'] = mv_db1['Cluster'].astype(str)
        cols_to_filter = ['Cluster', 'Smoke', 'Gender', 'Histopathology', 'Molecular_clusters']
        available_columns = ['Cluster', 'Smoke', 'Gender', 'Histopathology', 'Molecular_clusters', 'EMT', 'NE', 'NAPY',
                             'Histopathology_simplified']

    elif dataset_option == 'Rousseaux-Mixed(non-NE)':
        mv_db = clinic_df
        anno_file = anno_file
        mv_db.replace(
            ['[Not Available]', '[Unknown]', '[Not Evaluated]', 'nan', '[Not Applicable]', '[Discrepancy]', 'NaN'],
            value=None, inplace=True)
        mv_db.dropna(how='all', inplace=True, axis=1)
        mv_db.dropna(how='all', inplace=True)
        mv_db = mv_db.astype(str)
        mv_db = mv_db[mv_db.index.isin(selected_genes.index)]
        mv_db1.index = mv_db.index
        mv_db = pd.merge(mv_db, anno_file, left_index=True, right_index=True)
        mv_db1.index.rename("Patient ID", inplace=True)
        cluster["Cluster"] = cluster["Cluster"].astype('str')
        mv_db1["Cluster"] = cluster["Cluster"]
        mv_db1['Time'] = mv_db["Time"]
        mv_db1['Time'] = mv_db1['Time'].astype('float')
        mv_db1.dropna(how='any', inplace=True)
        mv_db1['Time'] = mv_db1['Time'].astype('int')
        mv_db1['Time'] = mv_db1['Time'].astype('str')
        mv_db1['Status'] = mv_db["Status"]
        mv_db1['EMT'] = mv_db["EMT"]
        mv_db1['NE'] = mv_db["NE"]
        mv_db1['NAPY'] = mv_db["NAPY"]
        mv_db1['Gender'] = mv_db["Gender"]
        mv_db1['SurvSig_Mol.Subtype'] = mv_db["SurvSig_Mol.Subtype"]
        mv_db1['Histology'] = mv_db["Histology"]
        mv_db1['Cluster'] = mv_db1['Cluster'].astype(str)
        cols_to_filter = ['Cluster', 'Gender', 'Histology', 'NE', 'EMT']
        available_columns = ['Cluster', 'Gender', 'Histology', 'NE', 'EMT', 'NAPY', 'SurvSig_Mol.Subtype']

    elif dataset_option == 'TCGA':
        if clinic_df.index.name != 'bcr_patient_barcode':
            clinic_df.set_index('bcr_patient_barcode', inplace=True)
        mv_db = clinic_df
        columns_to_drop_anno = ['ST-Selected', 'ST-other', 'Hist. Grade', 'Hist. Type', 'AJCC Stage', 'Clinical Stage']
        columns_to_drop_anno = [col for col in columns_to_drop_anno if col in anno_file.columns]
        anno_file.drop(columns=columns_to_drop_anno, inplace=True)
        mv_db = mv_db.astype(str)
        mv_db = mv_db[mv_db.index.isin(selected_genes.index)]
        anno_file = anno_file[anno_file.index.isin(selected_genes.index)]

        anno_file.dropna(how='all', axis=1, inplace=True)
        columns_to_drop = [col for col in anno_file.columns if anno_file[col].nunique() < 2]
        anno_file = anno_file.drop(columns=columns_to_drop)
        anno_file.rename(columns=lambda x: x.replace('ST-', 'ST_') if x.startswith('ST-') else x, inplace=True)

        mv_db.drop(
            columns=['type', 'age_at_initial_pathologic_diagnosis', 'gender', 'race', 'initial_pathologic_dx_year',
                     'menopause_status', 'birth_days_to', 'vital_status', 'tumor_status', 'last_contact_days_to',
                     'death_days_to', 'cause_of_death', 'new_tumor_event_type', 'new_tumor_event_site',
                     'new_tumor_event_site_other', 'new_tumor_event_dx_days_to', 'treatment_outcome_first_course',
                     'margin_status', 'residual_tumor', 'Redaction'], inplace=True)

        mv_db = pd.merge(mv_db, anno_file, left_index=True, right_index=True)
        mv_db.replace(
            ['[Not Available]', '[Unknown]', '[Not Evaluated]', 'nan', '[Not Applicable]', '[Discrepancy]', 'NaN'],
            value=None, inplace=True)
        mv_db.index.rename("Patient ID", inplace=True)
        cluster["Cluster"] = cluster["Cluster"].astype('str')
        mv_db["Cluster"] = cluster["Cluster"]
        mv_db.rename(columns={"histological_type": "Hist.type",
                              "histological_grade": "Hist.grade", "ajcc_pathologic_tumor_stage": "AJCC_Stage",
                              'clinical_stage': 'Clinical_Stage'}, inplace=True)

        mv_db.dropna(how='all', inplace=True, axis=1)
        mv_db.dropna(how='all', inplace=True)
        mv_db['Time'] = mv_db['Time'].astype(float)
        mv_db['Status'] = mv_db['Status'].astype(str)
        mv_db.replace('Stage II', '2', inplace=True)
        cols_to_filter = ['Hist.type', 'Cluster', 'AJCC_Stage', 'Hist.grade', 'Clinical_Stage']
        available_columns = ['Hist.type', 'Cluster', 'AJCC_Stage', 'Hist.grade', 'Clinical_Stage']
        cols_to_filter = cols_to_filter + anno_file.columns.to_list()
        mv_db1 = mv_db.copy()

    else:
        mv_db1 = None
        cols_to_filter = None
        available_columns = None
        pass

    return mv_db1, cols_to_filter, available_columns


def sg_ssgsea_multi_data(dataset_option, mv_group_df, surv_df, clinic_df, anno_file):
    mv_db1 = pd.DataFrame()
    if dataset_option == "George-SCLC":
        mv_db = pd.read_csv("source_data/sclc_ucologne_2015_clinical_data.tsv", delimiter="\t", index_col=1)
        anno_file = pd.read_csv("source_data/george_sclc_anno.csv", index_col=0)
        mv_db = mv_db.astype(str)
        mv_db = mv_db[mv_db.index.isin(surv_df.index)]
        mv_db1.index = mv_db.index
        mv_group_df["Cluster"] = mv_group_df["Cluster"].astype('str')
        mv_db1["Cluster"] = mv_group_df["Cluster"]
        mv_db1['Time'] = mv_db["Overall Survival (Months)"]
        mv_db1['Status'] = mv_db["Overall Survival Status"]
        mv_db1['Chemotherapy'] = mv_db['Chemotherapy']
        mv_db1['Ethnicity_Category'] = mv_db['Ethnicity Category']
        mv_db1['Sex'] = mv_db['Sex']
        mv_db1['EMT'] = anno_file['EMT']
        mv_db1['NE'] = anno_file['NE']
        mv_db1['NAPY'] = anno_file['NAPY']
        mv_db1['UICC_stage'] = mv_db['UICC Tumor Stage']
        mv_db1 = mv_db1.replace(['0:LIVING', '1:DECEASED'], ['0', '1'])
        mv_db1 = mv_db1.replace(['IA', 'IB'], ['Ia', 'Ib'])
        mv_db1 = mv_db1.replace(['NaN'], ['nan'])
        mv_db1['Cluster'] = mv_db1['Cluster'].astype(str)
        cols_to_filter = ['UICC_stage', 'Cluster', 'NE', 'EMT', 'NAPY']
        available_columns = ['UICC_stage', 'Cluster', 'Chemotherapy', 'Sex', 'NE', 'EMT', 'NAPY', 'Ethnicity_Category']

    elif dataset_option == 'Rousseaux-Mixed(non-NE)':
        clinic_df.set_index('sample', inplace=True)
        mv_db = clinic_df
        mv_db.replace(
            ['[Not Available]', '[Unknown]', '[Not Evaluated]', 'nan', '[Not Applicable]', '[Discrepancy]', 'NaN'],
            value=None, inplace=True)
        mv_db.dropna(how='all', inplace=True, axis=1)
        mv_db.dropna(how='all', inplace=True)
        mv_db = mv_db.astype(str)
        mv_db = mv_db[mv_db.index.isin(surv_df.index)]
        mv_db1.index = mv_db.index
        mv_db = pd.merge(mv_db, anno_file, left_index=True, right_index=True)
        mv_db1.index.rename("Patient ID", inplace=True)
        mv_db1['Time'] = mv_db["Time"]
        mv_db1['Time'] = mv_db1['Time'].astype('float')
        mv_db1.dropna(how='any', inplace=True)
        mv_db1['Time'] = mv_db1['Time'].astype('int')
        mv_db1['Time'] = mv_db1['Time'].astype('str')
        mv_db1['Status'] = mv_db["Status"]
        mv_db1['Histology'] = mv_db["Histology"]
        mv_db1['EMT'] = mv_db["EMT"]
        mv_db1['NE'] = mv_db["NE"]
        mv_db1['NAPY'] = mv_db["NAPY"]
        mv_db1['SurvSig_Mol.Subtype'] = mv_db["SurvSig_Mol.Subtype"]
        mv_db1['Gender'] = mv_db["Gender"]
        mv_group_df["Cluster"] = mv_group_df["Cluster"].astype('str')
        mv_db1["Cluster"] = mv_group_df["Cluster"]
        mv_db1['Cluster'] = mv_db1['Cluster'].astype(str)
        cols_to_filter = ['Cluster', 'Gender', 'Histology', 'NE', 'EMT']
        available_columns = ['Cluster', 'Gender', 'Histology', 'NE', 'EMT', 'NAPY', 'SurvSig_Mol.Subtype']

    elif dataset_option == 'Anish-SCLC':
        mv_db = pd.read_excel("source_data/Supplementary_table_patient_clinical_characteristics_2021015.xlsx",
                              index_col=2)
        anno_file = pd.read_csv("source_data/anish_anno.csv", index_col=0)
        mv_db = mv_db.astype(str)
        mv_db = mv_db[mv_db.index.isin(surv_df.index)]
        mv_db1.index = mv_db.index
        mv_db1.index.rename("Patient ID", inplace=True)
        mv_group_df["Cluster"] = mv_group_df["Cluster"].astype('str')
        mv_db1["Cluster"] = mv_group_df["Cluster"]
        mv_db1['Time'] = mv_db["Overall survival since diagnosis"]
        mv_db1['Time'] = mv_db1['Time'].astype('float')
        mv_db1['Time'] = mv_db1['Time'].astype('int')
        mv_db1['Time'] = mv_db1['Time'].astype('str')
        mv_db1['Status'] = mv_db["Mortality status"]
        mv_db1 = mv_db1.replace(['Dead', 'Alive'], ['1', '0'])
        mv_db1['Immunotherapy'] = mv_db["Immunotherapy prior to biopsy"]
        mv_db1['Smoking history'] = mv_db["Smoking history"]
        mv_db1['EMT'] = anno_file['EMT']
        mv_db1['NE'] = anno_file['NE']
        mv_db1['NAPY'] = anno_file['NAPY']
        mv_db1['Sex'] = mv_db["Sex"]
        mv_db1['SCLC_staging'] = mv_db['SCLC staging/initial diagnosis ']
        mv_db1['Cluster'] = mv_db1['Cluster'].astype(str)
        cols_to_filter = ['Cluster', 'SCLC_staging', 'Immunotherapy', 'Sex']
        available_columns = ['Cluster', 'SCLC_staging', 'Immunotherapy', 'Sex', 'EMT', 'NAPY', 'NE']

    elif dataset_option == 'Jiang-SCLC':
        mv_db = pd.read_csv("source_data/Jiang_SCLC_annotation.csv", index_col=0)
        mv_db = mv_db.astype(str)
        mv_db = mv_db[mv_db.index.isin(surv_df.index)]
        mv_db1.index = mv_db.index
        mv_db1.index.rename("Patient ID", inplace=True)
        mv_group_df["Cluster"] = mv_group_df["Cluster"].astype('str')
        mv_db1["Cluster"] = mv_group_df["Cluster"]
        mv_db1['Time'] = mv_db["Months"]
        mv_db1['Time'] = mv_db1['Time'].astype('float')
        mv_db1.dropna(how='any', inplace=True)
        mv_db1['Time'] = mv_db1['Time'].astype('int')
        mv_db1['Time'] = mv_db1['Time'].astype('str')
        mv_db1['Status'] = mv_db["Survival Status"]
        mv_db1['Chemotherapy'] = mv_db["chemotherapy-exposed biopsy"]
        mv_db1['Smoke'] = mv_db["Smoke"]
        mv_db1['Gender'] = mv_db["Gender"]
        mv_db1['SCLC_staging'] = mv_db['UICC stage']
        mv_db1['Cluster'] = mv_db1['Cluster'].astype(str)
        cols_to_filter = ['Cluster', 'SCLC_staging', 'Chemotherapy', 'Gender']
        available_columns = ['Cluster', 'SCLC_staging', 'Chemotherapy', 'Gender']

    elif dataset_option == 'George-LCNEC':
        mv_db = pd.read_csv("source_data/LCNEC_George_clinical_anno.csv", index_col=0)
        mv_db.replace(
            ['[Not Available]', '[Unknown]', '[Not Evaluated]', 'nan', '[Not Applicable]', '[Discrepancy]', 'NaN'],
            value=None, inplace=True)
        mv_db.dropna(how='all', inplace=True, axis=1)
        mv_db.dropna(how='all', inplace=True)
        mv_db = mv_db.astype(str)
        mv_db = mv_db[mv_db.index.isin(surv_df.index)]
        mv_db1.index = mv_db.index
        mv_db1.index.rename("Patient ID", inplace=True)
        mv_group_df["Cluster"] = mv_group_df["Cluster"].astype('str')
        mv_db1["Cluster"] = mv_group_df["Cluster"]
        mv_db1['Time'] = mv_db["Survival_months"]
        mv_db1['Time'] = mv_db1['Time'].astype('float')
        mv_db1.dropna(how='any', inplace=True)
        mv_db1['Time'] = mv_db1['Time'].astype('int')
        mv_db1['Time'] = mv_db1['Time'].astype('str')
        mv_db1['Status'] = mv_db["Survival_censor"]
        mv_db1 = mv_db1.replace(['dead', 'alive'], ['1', '0'])
        mv_db1['Molecular_Subtypes'] = mv_db["Molecular_Subtypes"]
        mv_db1['Smoke'] = mv_db["Smoking_Hx"]
        mv_db1['Gender'] = mv_db["Sex"]
        mv_db1['SCLC_staging'] = mv_db['Stage_UICC']
        mv_db1['Classification'] = mv_db['Classification']
        mv_db1['Cluster'] = mv_db1['Cluster'].astype(str)
        cols_to_filter = ['Cluster', 'SCLC_staging', 'Gender', 'Classification']
        available_columns = ['Cluster', 'SCLC_staging', 'Gender', 'Classification']

    elif dataset_option == 'Fernandez-Carcinoid':
        mv_db = pd.read_csv("source_data/carcinoid_clinical_anno.csv", index_col=0)
        mv_db.replace(
            ['[Not Available]', '[Unknown]', '[Not Evaluated]', 'nan', '[Not Applicable]', '[Discrepancy]', 'NaN'],
            value=None, inplace=True)
        mv_db.dropna(how='all', inplace=True, axis=1)
        mv_db.dropna(how='all', inplace=True)
        mv_db = mv_db.astype(str)
        mv_db = mv_db[mv_db.index.isin(surv_df.index)]
        mv_db1.index = mv_db.index
        mv_db1.index.rename("Patient ID", inplace=True)
        mv_group_df["Cluster"] = mv_group_df["Cluster"].astype('str')
        mv_db1["Cluster"] = mv_group_df["Cluster"]
        mv_db1['Time'] = mv_db["Survival_months"]
        mv_db1['Time'] = mv_db["Survival_months"]
        mv_db1['Time'] = mv_db1['Time'].astype('float')
        mv_db1.dropna(how='any', inplace=True)
        mv_db1['Time'] = mv_db1['Time'].astype('int')
        mv_db1['Time'] = mv_db1['Time'].astype('str')
        mv_db1['Status'] = mv_db["Survival_censor"]
        mv_db1 = mv_db1.replace(['dead', 'alive'], ['1', '0'])
        mv_db1['Histopathology'] = mv_db["Histopathology"]
        mv_db1['cluster_LNET'] = mv_db["cluster_LNET"]
        mv_db1['cluster_LNEN'] = mv_db["cluster_LNEN"]
        mv_db1['Smoke'] = mv_db["Smoking_status"]
        mv_db1['Gender'] = mv_db["Sex"]
        mv_db1['SCLC_staging'] = mv_db['Stage_UICC']
        mv_db1['Cluster'] = mv_db1['Cluster'].astype(str)
        cols_to_filter = ['Cluster', 'SCLC_staging', 'Gender', 'Histopathology', 'cluster_LNET']
        available_columns = ['Cluster', 'SCLC_staging', 'Gender', 'Histopathology', 'cluster_LNET']

    elif dataset_option == 'Liu-SCLC':
        mv_db = pd.read_csv("source_data/liu_clinical_anno.csv", index_col=0)
        mv_db.replace(
            ['[Not Available]', '[Unknown]', '[Not Evaluated]', 'nan', '[Not Applicable]', '[Discrepancy]', 'NaN',
             "unknown"], value=None, inplace=True)
        mv_db.dropna(how='all', inplace=True, axis=1)
        mv_db.dropna(how='all', inplace=True)
        mv_db.index = mv_db.index.str.upper()
        mv_db['Survial.(months)'] = mv_db['Survial.(months)'].str.split(r' \(').str[0]
        mv_db = mv_db.applymap(convert_string_to_float)
        mv_db = mv_db.astype(str)
        mv_db = mv_db[mv_db.index.isin(surv_df.index)]
        mv_db1.index = mv_db.index
        mv_db1.index.rename("Patient ID", inplace=True)
        mv_group_df["Cluster"] = mv_group_df["Cluster"].astype('str')
        mv_db1["Cluster"] = mv_group_df["Cluster"]
        mv_db1['Time'] = mv_db["Survial.(months)"]
        mv_db1['Time'] = mv_db1['Time'].astype('float')
        mv_db1.dropna(how='any', inplace=True)
        mv_db1['Time'] = mv_db1['Time'].astype('int')
        mv_db1['Time'] = mv_db1['Time'].astype('str')
        mv_db1['Status'] = mv_db["Status."]
        mv_db1 = mv_db1.replace(['dead', 'alive'], ['1', '0'])
        mv_db1['Smoke'] = mv_db["Smoking"]
        mv_db1['Tumor_site'] = mv_db["Tumor.site"]
        mv_db1['Histologic_type'] = mv_db["Histologic.type"]
        mv_db1['Gender'] = mv_db["Gender"]
        mv_db1['Stage'] = mv_db['TNM.Stage']
        mv_db1['Cluster'] = mv_db1['Cluster'].astype(str)
        cols_to_filter = ['Cluster', 'Stage', 'Gender', 'Tumor_site', 'Smoke']
        available_columns = ['Cluster', 'Stage', 'Gender', 'Tumor_site', 'Smoke']

    elif dataset_option == 'Alcala-Carcinoid':
        mv_db = pd.read_csv("source_data/alcala_clinical_anno.csv", index_col=0)
        mv_db.replace(
            ['[Not Available]', '[Unknown]', '[Not Evaluated]', 'nan', '[Not Applicable]', '[Discrepancy]', 'NaN'],
            value=None, inplace=True)
        mv_db.dropna(how='all', inplace=True, axis=1)
        mv_db.dropna(how='all', inplace=True)
        mv_db = mv_db.astype(str)
        mv_db = mv_db[mv_db.index.isin(surv_df.index)]
        mv_db1.index = mv_db.index
        mv_db1.index.rename("Patient ID", inplace=True)
        mv_group_df["Cluster"] = mv_group_df["Cluster"].astype('str')
        mv_db1["Cluster"] = mv_group_df["Cluster"]
        mv_db1['Time'] = mv_db["Survival_months"]
        mv_db1['Time'] = mv_db1['Time'].astype('float')
        mv_db1.dropna(how='any', inplace=True)
        mv_db1['Time'] = mv_db1['Time'].astype('int')
        mv_db1['Time'] = mv_db1['Time'].astype('str')
        mv_db1['Status'] = mv_db["Survival_censor"]
        mv_db1 = mv_db1.replace(['dead', 'alive'], ['1', '0'])
        mv_db1['Histopathology'] = mv_db["Histopathology"]
        mv_db1['Histopathology_simplified'] = mv_db["Histopathology_simplified"]
        mv_db1['Molecular_clusters'] = mv_db["Molecular_clusters"]
        mv_db1['Smoke'] = mv_db["Smoking_status"]
        mv_db1['Gender'] = mv_db["Sex"]
        mv_db1['Cluster'] = mv_db1['Cluster'].astype(str)
        cols_to_filter = ['Cluster', 'Smoke', 'Gender', 'Histopathology', 'Molecular_clusters']
        available_columns = ['Cluster', 'Smoke', 'Gender', 'Histopathology', 'Molecular_clusters',
                             'Histopathology_simplified']

    elif dataset_option == 'TCGA':
        mv_db = clinic_df
        columns_to_drop_anno = ['ST-Selected', 'ST-other', 'Hist. Grade', 'Hist. Type', 'AJCC Stage', 'Clinical Stage']
        columns_to_drop_anno = [col for col in columns_to_drop_anno if col in anno_file.columns]
        anno_file.drop(columns=columns_to_drop_anno, inplace=True)
        clinic_df.set_index('sample', inplace=True)
        mv_db = mv_db.astype(str)
        mv_db = mv_db[mv_db.index.isin(surv_df.index)]
        anno_file.index.rename('sample', inplace=True)
        anno_file = anno_file[anno_file.index.isin(surv_df.index)]
        anno_file.dropna(how='all', axis=1, inplace=True)
        columns_to_drop = [col for col in anno_file.columns if anno_file[col].nunique() < 2]
        anno_file = anno_file.drop(columns=columns_to_drop)
        anno_file.rename(columns=lambda x: x.replace('ST-', 'ST_') if x.startswith('ST-') else x, inplace=True)

        mv_db.replace(
            ['[Not Available]', '[Unknown]', '[Not Evaluated]', 'nan', '[Not Applicable]', '[Discrepancy]', 'NaN'],
            value=None, inplace=True)

        mv_db.drop(
            columns=['type', 'age_at_initial_pathologic_diagnosis', 'gender', 'race', 'initial_pathologic_dx_year',
                     'menopause_status', 'birth_days_to', 'vital_status', 'tumor_status', 'last_contact_days_to',
                     'death_days_to', 'cause_of_death', 'new_tumor_event_type', 'new_tumor_event_site',
                     'new_tumor_event_site_other', 'new_tumor_event_dx_days_to', 'treatment_outcome_first_course',
                     'margin_status', 'residual_tumor', 'Redaction'], inplace=True)

        mv_db = pd.merge(mv_db, anno_file, left_index=True, right_index=True)

        mv_db.replace(
            ['[Not Available]', '[Unknown]', '[Not Evaluated]', 'nan', '[Not Applicable]', '[Discrepancy]', 'NaN'],
            value=None, inplace=True)
        mv_db.index.rename("Patient ID", inplace=True)
        mv_group_df["Cluster"] = mv_group_df["Cluster"].astype('str')
        mv_db["Cluster"] = mv_group_df["Cluster"]

        mv_db.rename(columns={"histological_type": "Hist.type",
                              "histological_grade": "Hist.grade", "ajcc_pathologic_tumor_stage": "AJCC_Stage",
                              'clinical_stage': 'Clinical_Stage'}, inplace=True)

        mv_db.dropna(how='all', inplace=True, axis=1)
        mv_db.dropna(how='all', inplace=True)
        mv_db['Time'] = mv_db['Time'].astype(float)
        mv_db['Status'] = mv_db['Status'].astype(str)
        mv_db.replace('Stage II', '2', inplace=True)
        cols_to_filter = ['Hist.type', 'Cluster', 'AJCC_Stage', 'Hist.grade', 'Clinical_Stage']
        available_columns = ['Hist.type', 'Cluster', 'AJCC_Stage', 'Hist.grade', 'Clinical_Stage']
        cols_to_filter = cols_to_filter + anno_file.columns.to_list()
        mv_db1 = mv_db

    else:
        mv_db1 = None
        cols_to_filter = None
        available_columns = None
        pass

    return mv_db1, cols_to_filter, available_columns
