import pandas as pd
import streamlit as st
import shutil
import os
import json
import plotly.express as px
from streamlit_lottie import st_lottie_spinner
from Scripts.rpy2_heatmap_plots import py_hm
from Scripts.user_interaction import naming
import numpy as np
from multiprocessing import Process, Queue
import Scripts.descriptions as dsc


# Wrapper function to handle the R task and send status to the queue
def run_py_hm(queue, *args):
    try:
        py_hm(*args)  # Run the R-related function
        queue.put(True)  # Indicate success
    except Exception as e:
        queue.put(False)  # Indicate failure
        print(f"Error in py_hm: {e}")


def file_copy(user_id):
    """
    Copy user guide to the user's result directory and create a zip archive.

    Arguments:
    user_id -- User ID for file naming and tracking.
    """
    output_filename = f"result_data/result_data_{user_id}"
    dir_name = f"result_data/result_data_{user_id}"

    source_path = "source_data/user_guide.txt"
    destination_path = f"{dir_name}/user_guide.txt"

    try:
        # Copy the file from the source to the destination
        shutil.copyfile(source_path, destination_path)

        # Create the zip archive
        shutil.make_archive(output_filename, 'zip', dir_name)
    except Exception as e:
        st.error(f"An error occurred while copying the file: {e}")


def downloader(user_id, dl_name, download_exp):
    """
    Provide a download button for the results.

    Arguments:
    user_id -- User ID for file naming and tracking.
    dl_name -- The name provided by the user for the results.
    download_exp -- The expander in which to place the download button.
    """
    try:
        # Paths
        output_filename = f"result_data/result_data_{user_id}"
        zip_file_path = f"{output_filename}.zip"

        dir_name = f"result_data/result_data_{user_id}"

        source_path = "source_data/user_guide.txt"
        destination_path = f"{dir_name}/user_guide.txt"

        # Copy the file from the source to the destination
        shutil.copyfile(source_path, destination_path)

        # Create the zip archive
        shutil.make_archive(output_filename, 'zip', dir_name)

        # Read the zip file contents
        if os.path.isfile(zip_file_path):
            with open(zip_file_path, "rb") as fp:
                zip_bytes = fp.read()

            # Provide the zip file for download
            btn = download_exp.download_button(
                label="Download Results",
                data=zip_bytes,
                file_name=dl_name + ".zip",
                mime="application/zip",
                key=f"dwnld_btn_{user_id}"
            )
            if btn:
                st.toast('Download Successful!', icon='🎉')
                # Optionally reset the flag to allow re-preparation
                # st.session_state['download_ready'] = False
        else:
            download_exp.error('The zip file could not be found. Please ensure you have generated the results.')
    except Exception as e:
        download_exp.error(f"An error occurred during the download process: {e}")


def rgb_to_hex(rgb):
    """
    Convert RGB color to HEX format.

    Arguments:
    rgb -- RGB color string.

    Returns:
    HEX color string.
    """
    try:
        if rgb.startswith('rgb'):
            # Split the RGB color values
            r, g, b = map(int, rgb.split('(')[1].split(')')[0].split(','))
            # Convert the RGB color values to hexadecimal format
            return f'#{r:02x}{g:02x}{b:02x}'
        else:
            return rgb
    except Exception as e:
        st.error(f"An error occurred while converting RGB to HEX: {e}")
        return rgb


def ht_adv_set(column_cluster_path, ht_row_cluster_path, user_session_id, ht_top_annotation_path, option_dataset,
               surv_dataframe_path, ht_png, ht_pdf, surv_png, surv_pdf, ht_expression_path, cancer_opt):
    """
    Generate advanced heatmap settings in Streamlit.

    Arguments:
    column_cluster_path -- Path to the CSV file for column clustering information.
    ht_row_cluster_path -- Path to the CSV file for row clustering information.
    user_session_id -- User session ID for file naming and tracking.
    ht_top_annotation_path -- Path to the CSV file for top annotation data.
    option_dataset -- Dataset option for the analysis.
    surv_dataframe_path -- Path to the CSV file containing survival data.
    ht_png -- Path to save the PNG output of the heatmap.
    ht_pdf -- Path to save the PDF output of the heatmap.
    surv_png -- Path to save the PNG output of the survival plot.
    surv_pdf -- Path to save the PDF output of the survival plot.
    ht_expression_path -- Path to the CSV file containing expression data.
    cancer_opt -- Option for specific cancer type.
    """
    try:
        df_cluster_number = pd.read_csv(column_cluster_path)
        num_of_cluster = df_cluster_number["Cluster"].nunique()

        df_cluster_number_gene = pd.read_csv(ht_row_cluster_path)
        clust_num_gene = df_cluster_number_gene["Cluster"].nunique()

        if num_of_cluster < 20 and clust_num_gene < 20:
            if (num_of_cluster > 1 and clust_num_gene > 0) or (
                    df_cluster_number_gene["Cluster"].unique()[0] == "Genes" and
                    num_of_cluster > 1):
                if clust_num_gene == 'Genes':
                    clust_num_gene = 1

                cluster_color = px.colors.qualitative.Set1 if num_of_cluster <= 9 else px.colors.qualitative.Light24
                cluster_color = cluster_color[:num_of_cluster]

                adv_ht_form = st.form(key="adv_ht_form", border=False)
                adv_ht_expander = adv_ht_form.expander(":blue[Advanced Setting of Plots]")

                adv_ht_expander.subheader(":blue[Heatmap's Settings]")

                perc_slider = adv_ht_expander.slider(':blue[Select a ] :red[Range of Values]', 0, 100, (1, 99),
                                                     key="perc_setting")

                exp_quantile = pd.read_csv(ht_expression_path, index_col=0).values
                z_score_check = exp_quantile.min()
                min_value_perc = np.quantile(exp_quantile, (perc_slider[0] / 100))
                max_value_perc = np.quantile(exp_quantile, (perc_slider[1] / 100))

                perc_mid = int(((perc_slider[1] - perc_slider[0]) / 2) + perc_slider[0])

                if z_score_check < 0:
                    if min_value_perc > 0:
                        adv_ht_expander.warning(
                            ":blue[The minimum value is ] :red[ positive number.] :blue[It can cause ] "
                            ":red[ misunderstandings ] :blue[in the result.] ")
                    elif max_value_perc < 0:
                        adv_ht_expander.warning(
                            ":blue[The maximum value is ] :red[ negative number.] :blue[It can cause ] "
                            ":red[ misunderstandings ] :blue[in the result.] ")

                perc_max = str(perc_slider[1])
                perc_min = str(perc_slider[0])
                perc_mid = str(perc_mid)
                perc_col1, perc_col2, perc_col3, perc_col4 = adv_ht_expander.columns(4)
                perc_col4.write(f':blue[Max: ] :red[{perc_slider[1]}]')
                perc_col3.write(f':blue[Mid: ] :red[{perc_mid}]')
                perc_col2.write(f':blue[Min: ] :red[{perc_slider[0]}]')

                hm_name = "Z-Score" if z_score_check < 0 else "Exp(log2)"
                hm_name = adv_ht_expander.text_input(":red[Name ] :blue[of your Heatmap]", value=hm_name, max_chars=20,
                                                     key="hm_name")

                adv_ht_expander.subheader(" ", divider="blue")
                adv_ht_expander.subheader(":blue[Heatmap's Colors]")
                color_col1, color_col2, color_col3 = adv_ht_expander.columns(3)
                color_max = color_col1.color_picker(":red[Max] :blue[ Value's Color]", value="#FF0000",
                                                    key="top_color_hm")
                color_mid = color_col2.color_picker(":red[Mid] :blue[ Value's Color]", value="#FFFFFF",
                                                    key="middle_color_hm")
                color_min = color_col3.color_picker(":red[Min] :blue[ Value's Color]", value="#0000FF",
                                                    key="bot_color_hm")

                clust_color_df = pd.DataFrame()
                clust_color_df["cluster"] = df_cluster_number["Cluster"].unique()
                clust_color_df = clust_color_df.sort_values(by="cluster", ascending=True)
                clust_color_df["color"] = cluster_color[:num_of_cluster]
                clust_color_df['color'] = clust_color_df['color'].apply(rgb_to_hex)

                cluster_color_path_streamlit = naming(user_session_id)[41]

                surv_clusters = pd.read_csv(surv_dataframe_path, index_col=1)
                surv_clusters = surv_clusters["Cluster"].unique()
                clust_color_df = clust_color_df[clust_color_df['cluster'].isin(surv_clusters)]

                adv_ht_expander.subheader(" ", divider="blue")
                adv_ht_expander.subheader(":blue[Display Heatmap's Dendogram]")

                dendrogram_cols = adv_ht_expander.columns(2)

                if df_cluster_number.shape[0] > 200:
                    default_dendro_sample = False
                else:
                    default_dendro_sample = True
                if df_cluster_number_gene.shape[0] > 400:
                    default_dendro_gene = False
                else:
                    default_dendro_gene = True

                dendrogram = dendrogram_cols[0].toggle(":blue[Use ] :red[Dendrogram ] :blue[for Heatmap's Samples]",
                                                      key="dendrogram_sample", value=default_dendro_sample)
                dendrogram_gene = dendrogram_cols[1].toggle(":blue[Use ] :red[Dendrogram ] :blue[for Heatmap's Genes]",
                                                           key="dendrogram_gene", value=default_dendro_gene)

                adv_ht_expander.subheader(" ", divider="blue")
                adv_ht_expander.subheader(":blue[Annotation's Colors]")
                annotations = pd.read_csv(ht_top_annotation_path, delimiter=',', index_col=0, dtype='object')
                df_cluster_number.set_index('Samples', inplace=True)
                annotations['Cluster'] = df_cluster_number['Cluster']
                df_cluster_number.reset_index(inplace=True)
                annotations.fillna('nan', axis=1, inplace=True)
                cluster_column = annotations.pop('Cluster')
                annotations.insert(0, 'Cluster', cluster_column)


                if option_dataset == 'TCGA':
                    annos_list = [
                        col
                        for col in annotations.columns.tolist()[1:]
                        if annotations[col].nunique() > 1
                    ]
                else:
                    annos_list = [
                        col
                        for col in annotations.columns.tolist()[1:]
                    ]

                if option_dataset == 'George-SCLC':
                    default_anno = ['NAPY', 'NE', 'EMT']
                elif option_dataset == 'Anish-SCLC':
                    default_anno = ['NAPY', 'NE', 'EMT']
                elif option_dataset == 'Liu-SCLC':
                    default_anno = ['NAPY', 'NE', 'EMT']
                elif option_dataset == 'Jiang-SCLC':
                    default_anno = ['NAPY', 'NE', 'EMT']
                elif option_dataset == 'George-LCNEC':
                    default_anno = ["NAPY", "NE", 'Molecular_Subtypes', 'EMT']
                elif option_dataset == 'Fernandez-Carcinoid':
                    default_anno = ["NAPY", 'cluster_LNET', 'EMT']
                elif option_dataset == 'Alcala-Carcinoid':
                    default_anno = ["NAPY", 'Histopathology', 'EMT', 'NE', 'Molecular_clusters']
                elif option_dataset == 'Rousseaux-Mixed(non-NE)':
                    default_anno = ["NAPY", "NE", 'Histology', 'EMT']
                else:
                    if cancer_opt == 'ACC':
                        default_anno = ['ST-DNAmeth', 'ST-mRNA', 'ST-other']
                    elif cancer_opt == 'BLCA':
                        default_anno = ['ST-Selected']
                    elif cancer_opt == 'BRCA':
                        default_anno = ['ST-mRNA', 'Hist. Type']
                    elif cancer_opt == 'CESC':
                        default_anno = ['Hist. Type']
                    elif cancer_opt == 'CHOL':
                        default_anno = ['AJCC Stage', 'Hist. Type', 'Hist. Grade']
                    elif cancer_opt == 'COAD':
                        default_anno = ['Hist. Type', 'ST-Selected']
                    elif cancer_opt == 'DLBC':
                        default_anno = ['Hist. Type']
                    elif cancer_opt == 'ESCA':
                        default_anno = ['Hist. Type', 'ST-other']
                    elif cancer_opt == 'GBM':
                        default_anno = ['ST-mRNA', 'ST-other', 'ST-DNAmeth']
                    elif cancer_opt == 'HNSC':
                        default_anno = ['Hist. Grade', 'ST-other', 'ST-DNAmeth']
                    elif cancer_opt == 'KICH':
                        default_anno = ['ST-other']
                    elif cancer_opt == 'KIRC':
                        default_anno = ['ST-miRNA', 'ST-mRNA']
                    elif cancer_opt == 'KIRP':
                        default_anno = ['ST-Integrative', 'ST-miRNA', 'ST-CNA']
                    elif cancer_opt == 'LAML':
                        default_anno = ['ST-Selected', 'ST-miRNA']
                    elif cancer_opt == 'LGG':
                        default_anno = ['ST-mRNA', 'ST-DNAmeth', 'ST-other']
                    elif cancer_opt == 'LIHC':
                        default_anno = ['ST-mRNA', 'ST-DNAmeth', 'ST-Selected']
                    elif cancer_opt == 'LUAD':
                        default_anno = ['ST-Integrative', 'ST-DNAmeth']
                    elif cancer_opt == 'LUSC':
                        default_anno = ['ST-mRNA']
                    elif cancer_opt == 'MESO':
                        default_anno = ['Hist. Type']
                    elif cancer_opt == 'OV':
                        default_anno = ['ST-mRNA', ]
                    elif cancer_opt == 'PAAD':
                        default_anno = ['Hist. Type']
                    elif cancer_opt == 'PCPG':
                        default_anno = ['ST-mRNA', 'ST-DNAmeth', 'ST-Selected',]
                    elif cancer_opt == 'PRAD':
                        default_anno = ['ST-CNA', 'ST-DNAmeth', 'ST-other', 'ST-mRNA', 'ST-miRNA']
                    elif cancer_opt == 'READ':
                        default_anno = ['Hist. Type', 'ST-other']
                    elif cancer_opt == 'SARC':
                        default_anno = ['Hist. Type']
                    elif cancer_opt == 'SKCM':
                        default_anno = ['ST-mRNA', 'ST-DNAmeth', 'ST-miRNA', 'ST-other']
                    elif cancer_opt == 'STAD':
                        default_anno = ['Hist. Type', 'ST-other']
                    elif cancer_opt == 'TGCT':
                        default_anno = ['Clinical Stage']
                    elif cancer_opt == 'THCA':
                        default_anno = ['Hist. Type', 'ST-DNAmeth', 'ST-miRNA', 'ST-Selected']
                    elif cancer_opt == 'THYM':
                        default_anno = ['Hist. Type']
                    elif cancer_opt == 'UCEC':
                        default_anno = ['Hist. Type', 'ST-mRNA', 'ST-CNA', 'ST-Integrative', 'ST-other']
                    elif cancer_opt == 'UCS':
                        default_anno = ['Hist. Type', 'ST-Selected']
                    elif cancer_opt == 'UVM':
                        default_anno = ['Hist. Type']
                    else:
                        default_anno = ['Gender']

                anno_options = adv_ht_expander.multiselect(":blue[Select ] :red[Annotations]", options=annos_list,
                                                           max_selections=100, key='ms_tcga', default=default_anno,
                                                           placeholder="Select Annotations for Heatmap")
                adv_ht_expander.info(
                    'If you add or remove annotations please click the :red["Refresh the Annotations"] '
                    'button!')
                refresh_btn = adv_ht_expander.form_submit_button("Refresh the Annotations", help=dsc.act_btn(),
                                                                 icon="🔃")

                if refresh_btn:
                    adv_ht_expander.info("The changes have been :red[SUCCESSFULLY] refreshed")

                if len(anno_options) == 0:
                    annotations = annotations[["Cluster"]]
                else:
                    anno_options.extend(['Cluster'])
                    annotations = annotations[anno_options]
                    annotations.transpose()
                    annotations.drop_duplicates(keep='first')
                    annotations.transpose()

                annotations = annotations.astype(str)
                cluster_column = annotations.pop('Cluster')
                annotations.insert(0, 'Cluster', cluster_column)
                annotations = annotations.select_dtypes(include=['object'])
                annotations.to_csv(ht_top_annotation_path)
                annotations = annotations.reset_index(drop=True)
                annotations_elements = annotations.values.flatten()
                annotations_elements = set(annotations_elements)
                annotations_elements = sorted(annotations_elements)

                object_columns = [col for col in annotations.columns if annotations[col].dtype == 'object']

                results = {'type': [], 'column': []}

                for elem in annotations_elements:
                    for col in object_columns:
                        if elem in annotations[col].values:
                            results['type'].append(elem)
                            results['column'].append(col)

                anno_colors = pd.DataFrame(results)
                anno_colors.sort_values(by=["column", "type"], ascending=[True, True], inplace=True)
                anno_colors.reset_index(inplace=True)
                anno_colors.drop('index', axis=1, inplace=True)

                color_pal = ['#0099FF', '#FF9900', '#66CC00', '#FF66CC', '#00CC66',
                             '#FF3300', '#CC00CC', '#33CCCC', '#FFCC00', '#99FF33',
                             '#FF0066', '#33CC00', '#FF99CC', '#00CCFF', '#CC3300',
                             '#66CCFF', "#CDEEAA", "#EFB960", "#C3AED6", "#874C62",
                             "#EA907A", "#2192FF", "#557571", "#F9FD50", "#B2EC5D"]

                all_color_palettes = []
                color_palette = []

                cluster_color_dict = {'Cluster': dict(zip(clust_color_df['cluster'], clust_color_df['color']))}

                type_color_mapping_groups = {
                    'Histology': {"ADC": "#1f77b4", "LCNEC": "#9beaf2", "CARCI": "#ff7f0e", "SQC": "#ffb5e9",
                                  "BAS": "#2ca02c", "LCC": "#7a37b8", "Other": "#87888a", "SCLC": "#d62728",
                                  "NTL": "#000000"},
                    'NAPY': {'NEUROD1': '#8FEBC8', 'ASCL1': '#FD8A8A', 'YAP1': '#EFFA7A', 'POU2F3': '#8E60DB'},
                    'NE': {'NE': '#FCCA53', 'non-NE': '#84CBD9'},
                    'Sex': {'Female': '#c79558', 'Male': '#7df8ff', 'F': '#c79558', 'M': '#7df8ff', 'nan': '#a9a9ab'},
                    'EMT': {'Epithelial': '#9DE481', 'Epithelial-Mesenchymal': '#D3B7D9', 'Mesenchymal': '#674E62'},
                    'Smoker': {'Current': '#F2A359', 'Former': '#64A5DE', 'Never': '#98D19C', 'nan': '#a9a9ab',
                               'NA_Smoke': '#a9a9ab', 'NA_Smoker': '#a9a9ab'},
                    'Smoke': {'No': '#98D19C', 'Yes': '#F2A359', 'nan': '#a9a9ab', 'Current': '#F2A359',
                              'Former': '#64A5DE', 'Never': '#98D19C', 'Passive': '#ae95b8', 'no': '#98D19C',
                              'yes': '#F2A359', 'unknown': '#a9a9ab', 'NA_Smoke': '#a9a9ab', 'NA_Smoker': '#a9a9ab'},
                    'Tumor_Stage': {'Extensive stage': '#f74e05', 'Limited stage': '#fcacf0'},
                    'UICC': {'I': '#7AB1E7', 'IA': '#8FBCEC', 'IB': '#A4C7F1', 'II': '#B9D2F6', 'IIA': '#CEDDFB',
                             'IIB': '#E3E8FF', 'IIIA': '#D1C2D2', 'IIIB': '#BF9CAB', 'IV': '#AD769E', 'nan': '#a9a9ab',
                             'NA_UICC': '#a9a9ab'},
                    'UICC stage': {'I': '#7AB1E7', 'IA': '#8FBCEC', 'IB': '#A4C7F1', 'II': '#B9D2F6', 'IIA': '#CEDDFB',
                                   'IIB': '#E3E8FF', 'IIIA': '#D1C2D2', 'IIIB': '#BF9CAB', 'IV': '#AD769E',
                                   'IIIA/IIIB': "#c0a7c2", 'NA_UICC': '#a9a9ab'},
                    'UICC Stage': {'Ia': '#8FBCEC', 'Ib': '#A4C7F1', 'IIa': '#CEDDFB',
                                   'IIb': '#E3E8FF', 'III': "#D1C2D2", 'IIIa': '#c0a7c2', 'IIIb': '#BF9CAB',
                                   'IV': '#AD769E',
                                   'nan': '#a9a9ab', 'I': '#7AB1E7', 'IA': '#8FBCEC', 'IA1': '#A4C7F1',
                                   'IA2': '#B9D2F6', 'NA_UICC': '#a9a9ab',
                                   'IA3': '#CEDDFB', 'IB': '#E3E8FF', 'IIA': '#D1C2D2', 'IIB': '#BF9CAB',
                                   'IIIA': '#AD769E'},
                    'Gender': {'MALE': '#7df8ff', 'FEMALE': '#c79558', 'M': '#7df8ff', 'F': '#c79558',
                               'Male': '#7df8ff',
                               'Female': '#c79558'},
                    'Race': {'[Not Available]': '#a9a9ab', '[Not Evaluated]': '#a9a9ab', '[Unknown]': '#a9a9ab',
                             'AMER. IND. / ALASKA NATIVE': '#B3FFD9', 'ASIAN': '#CCE6FF',
                             'BLACK / AFR. AMER.': '#D9B3FF', 'NA_Race': '#a9a9ab',
                             'NAT. HAWAIIAN / PACIFIC ISLANDER': '#FFCCE6', 'WHITE': '#FFD9B3', 'Asian': '#CCE6FF',
                             'Black or African American': '#D9B3FF', 'White': '#FFD9B3', 'nan': '#a9a9ab'},
                    'Chemotherapy': {'Naive': '#98D19C', 'Chemo': '#F2A359'},
                    'TNM_Stage': {'Ia2': '#8FBCEC', 'Ia3': '#A4C7F1', 'Ib': '#CEDDFB', 'IIa': '#E3E8FF',
                                  'IIb': '#D1C2D2',
                                  'IIIa': '#c0a7c2', 'IIIb': '#BF9CAB', 'IV': '#AD769E'},
                    'Histopathology_simplified': {"Typical_HS": "#FF9900", "Atypical_HS": "#0099FF",
                                                  'Supra_carcinoid_HS': '#CC00CC'},
                    'Histopathology': {"Typical": "#FF9900", "Atypical": "#0099FF", 'Carcinoid': '#CC00CC'},
                    'Ethnicity Category': {'Asian': '#CCE6FF', 'Caucasian': '#FFD9B3', 'NA_EC': '#a9a9ab'},
                    'Molecular_Subtypes': {'NA_Mol_ST': '#a9a9ab', 'Type 1 LCNEC': '#00A20B', 'Type 2 LCNEC': '#AF5600'},
                    'Merged_Mol_clusters': {'A1_mol_sclust': '#0099FF', 'A2_mol_clust': '#FF9900',
                                            'B_mol_clust': '#66CC00', 'Supra_carcinoid_mol_clust': '#CC00CC'},
                    'Molecular_clusters': {'Carcinoid-A1': '#0099FF', 'Carcinoid-A2': '#FF9900',
                                            'Carcinoid-B': '#66CC00', 'Supra_carcinoid_MC': '#CC00CC',
                                           'LC1': '#FF0009', 'LC2': '#45029C', 'LC3': '#046220'}

                }

                type_color_mapping = {**type_color_mapping_groups, **cluster_color_dict}

                unknown_color = "#a9a9ab"

                for i in range(anno_colors['column'].nunique()):
                    column = anno_colors['column'].unique()[i]
                    unique_values = anno_colors.loc[anno_colors['column'] == column, 'type'].unique()

                    if column in type_color_mapping:
                        if isinstance(type_color_mapping[column], dict):
                            # Update existing dictionary-based mapping
                            value_color_mapping = type_color_mapping[column]
                            # For keys that start with certain patterns, force them to unknown_color
                            for key in list(value_color_mapping.keys()):
                                if any(key.lower().startswith(pattern)
                                       for pattern in ["[not available]",
                                                       "[unknown]",
                                                       "[not evaluated]",
                                                       "nan",
                                                       "[not applicable]",
                                                       "[discrepancy]"]):
                                    value_color_mapping[key] = unknown_color

                            # Build the color_palette using our updated mappings
                            color_palette = [value_color_mapping.get(val, unknown_color)
                                             for val in unique_values]
                        else:
                            # If the mapping isn’t a dict, just use the existing palette
                            color_palette = type_color_mapping[column]
                    else:
                        # Create a new mapping if the column is not in type_color_mapping
                        # For values that *start with* certain patterns, use unknown_color
                        value_color_mapping = {
                            val: (color_pal[i % len(color_pal)]
                                  if not any(val.lower().startswith(pattern)
                                             for pattern in ["[not available]",
                                                             "[unknown]",
                                                             "[not evaluated]",
                                                             "[not applicable]",
                                                             "[discrepancy]",
                                                             "nan",
                                                             "ntl_m",
                                                             "ntl_n",
                                                             "ntl_t_stage",
                                                             "ntl_gender",
                                                             "ntl_tnm",
                                                             "non_CARCI"])
                                  else unknown_color)
                            for i, val in enumerate(unique_values)
                        }

                        color_palette = [value_color_mapping[val] for val in unique_values]

                    color_palette = [rgb_to_hex(rgb) for rgb in color_palette]
                    all_color_palettes.extend(color_palette)

                anno_color_col1, anno_color_col2 = adv_ht_expander.columns(2)

                default_color = all_color_palettes

                unique_columns = anno_colors['column'].nunique()
                adv_ht_expander.header(" ", divider="blue")
                cols = adv_ht_expander.columns(unique_columns)
                color_list = []
                total_count = 0
                for i, group in enumerate(anno_colors.groupby('column')):
                    count_j = 0
                    cols[i].subheader(f':blue[{group[0]}]')
                    cols[i].header(" ", divider="blue")
                    for j, (_, row) in enumerate(group[1].iterrows()):
                        current_color = default_color[j + total_count]
                        count_j += 1
                        color_anno = cols[i].color_picker(row['type'], value=current_color,
                                                          key=f"color_picker_{group[0]}_{j}")
                        color_list.append(color_anno)

                    total_count += count_j

                adv_ht_expander.header(" ", divider="blue")

                anno_colors.set_index("column", inplace=True)

                anno_colors['color'] = color_list
                anno_colors['color'] = anno_colors['color'].apply(str.upper)
                anno_colors_path = naming(user_session_id)[42]

                sort_order = {'Cluster': 1}
                anno_colors['Sort_Column'] = anno_colors.index.map(sort_order)
                anno_colors.loc[:, 'Sort_Column'] = anno_colors['Sort_Column'].fillna(999)

                anno_colors.sort_values(by=['Sort_Column', 'column'], inplace=True)
                anno_colors.drop('Sort_Column', axis=1, inplace=True)
                anno_colors.to_csv(anno_colors_path)

                color_mapping = anno_colors.set_index('type')['color'].to_dict()

                clust_color_df['color'] = clust_color_df['cluster'].map(color_mapping)

                clust_color_df.to_csv(cluster_color_path_streamlit, index=False)

                adv_ht_expander.info('If you have entered the heatmap settings, please click on the '
                                     ':red["Accept All Changes"] button to accept the changes. ')

                acpt_btn = adv_ht_expander.form_submit_button("Accept All Changes and Draw Plots", help=dsc.act_btn(),
                                                              icon="✏️")

                if acpt_btn:
                    adv_ht_expander.info("The changes have been :red[SUCCESSFULLY] updated")

                hm_button = st.sidebar.button("Heatmap and Survival Plot", icon="✏️", key="hm_drawer")

                if hm_button or acpt_btn:
                    lottie_url = "style_items/surv_sig_lottie.json"
                    with open(lottie_url, "r") as f:
                        animation = json.load(f)

                    # Use the Queue in your Streamlit app
                    with st_lottie_spinner(animation, height=200, width=200, key="plot_spinner"):
                        # Create a Queue for communication
                        queue = Queue()

                        # Create and start the process
                        r_process = Process(
                            target=run_py_hm,
                            args=(
                                queue,  # Pass the Queue as the first argument
                                perc_max, perc_mid, perc_min, color_max, color_mid, color_min, hm_name,
                                cluster_color_path_streamlit, anno_colors_path, ht_expression_path,
                                column_cluster_path, ht_top_annotation_path, surv_dataframe_path, ht_png,
                                ht_pdf, surv_png, surv_pdf, user_session_id, ht_row_cluster_path,
                                option_dataset, dendrogram, dendrogram_gene
                            )
                        )
                        r_process.start()
                        r_process.join()  # Wait for the process to finish

                        # Get the result from the Queue
                        success = queue.get()

                        if success:
                            st.toast(':red[Heatmap and Survival Plot ] :blue[are Completed!]', icon='🎉')
                        else:
                            st.error(':red[Heatmap and Survival Plot ] :blue[encountered an error.]')

            else:
                st.error(
                    ":blue[Only ] :red[ONE ] :red[cluster. ] :blue[Please use advanced settings to ] :red[increase the number of clusters.]")
        else:
            st.error(
                ":blue[Too many clusters ] :red[(>20). ] :blue[Please use advanced settings to ] :red[increase the number of clusters.]")
    except Exception as e:
        st.error(f"An error occurred while generating advanced heatmap settings: {e}")
