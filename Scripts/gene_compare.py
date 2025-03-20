import pandas as pd
import streamlit as st
from Scripts.r_code_settings import file_copy
from scipy.stats import zscore
from scipy import stats
from scipy.stats import spearmanr, pearsonr
from statsmodels.stats.multitest import multipletests
from Scripts.user_interaction import naming
from Scripts.plotly_plots import plot_custom_fusion, plot_gene_correlation

@st.fragment
def gene_comp(user_id, df, anno_file, dataset_option, clinic_df, sample_cluster):
    """
    Performs gene comparison analysis, including z-score normalization, correlation calculation,
    and various plotting functionalities for gene data. It integrates clinical and annotation data
    for comprehensive analysis.

    Parameters:
    - user_id (str): User identifier for session-specific operations.
    - df (DataFrame): DataFrame containing gene expression data.
    - anno_file (DataFrame): DataFrame with annotation data.
    - dataset_option (str): Dataset selection option.
    - clinic_df (DataFrame): DataFrame containing clinical data.

    Returns:
    - None. The function performs analysis and visualization, displaying results in the Streamlit application.
    """
    try:
        settings_col = st.columns(4)
        axis_gene_col = st.columns(2)
        # Identify overlapping columns
        default_index = None
        corr_comp_df = None
        overlap_columns = anno_file.columns.intersection(clinic_df.columns)

        # Rename overlapping columns in clinic_df before merging
        clinic_df_renamed = clinic_df.rename(columns={col: f"{col}_from_clinic_df" for col in overlap_columns})

        # Merge anno_file and clinic_df_renamed
        merged_df = pd.merge(anno_file, clinic_df_renamed, left_index=True, right_index=True, how='outer')

        if sample_cluster is not None:
            merged_df['Sample_Cluster'] = sample_cluster['Cluster'].reindex(merged_df.index)

        for col in overlap_columns:
            merged_df.drop(col, axis=1, inplace=True)  # Drop anno_file's original overlapping columns
            merged_df.rename(columns={f"{col}_from_clinic_df": col}, inplace=True)  # Rename clinic_df's columns back

        settings_col[2].write(":blue[Use ] :red[**Z-Score:**]")
        z_score_toggle = settings_col[3].toggle(':blue[Activate feature]', key="ztoggc")
        if z_score_toggle:
            df_z_gc = zscore(df.copy(), axis=1, nan_policy='omit')
            df_z_gc = df_z_gc.dropna(axis=0)
        else:
            df_z_gc = df.copy()

        df = df_z_gc

        # For selectbox
        # Ensure your DataFrame index is a list or convert it to a list
        df_index = df.index.tolist()

        # Find the index of 'TP53' in the list
        try:
            selected_index = df_index.index('TP53')
        except ValueError:
            selected_index = 110  # Default to the first item if TP53 is not found

        selected_gene_x = axis_gene_col[0].selectbox(":blue[Select a ] :red[X Axis Gene]", options=df_index,
                                                     key="gene_selx", index=selected_index)
        axis_gene_col[0].caption(":blue[Check your ] :red[X Axis Gene]",
                                 help=f"https://www.genecards.org/Search/Keyword?queryString={selected_gene_x}")

        corr_comp_chb = settings_col[0].toggle(
            f":blue[Compare ] :red[All ] :blue[Genes with: ] :red[{selected_gene_x}]",
            key='corr_comp_chb')
        if corr_comp_chb:
            corr_meth = st.selectbox(":blue[Select the ] :red[Correlation Method]", options=["Spearman", "Pearson"],
                                     index=0, key='corr_meth')
            corr_comp_df = all_vs_one_correlation(df, selected_gene_x, corr_meth)
            corr_comp_df.set_index("Gene", inplace=True)

        if corr_comp_chb:
            selected_gene_y = axis_gene_col[1].selectbox(":blue[Select a ] :red[Y Axis Gene]",
                                                         options=corr_comp_df.index, index=0,
                                                         key="gene_sely")
        else:
            selected_gene_y = axis_gene_col[1].selectbox(":blue[Select a ] :red[Y Axis Gene]", options=df_index,
                                                         key="gene_sely",
                                                         index=500)

        axis_gene_col[1].caption(":blue[Check your ] :red[Y Axis Gene]",
                                 help=f"https://www.genecards.org/Search/Keyword?queryString={selected_gene_y}")

        st.subheader(" ", divider="blue")

        df_x_gene = df.loc[[selected_gene_x]]
        df_y_gene = df.loc[[selected_gene_y]]

        advanced_set_col = st.columns(2)

        gc_corr_meth = advanced_set_col[0].selectbox(":blue[Select the ] :red[Correlation Method]",
                                                     options=["Spearman", "Pearson"],
                                                     index=0, key='gc_corr_method')

        color_options = list(merged_df.columns)
        color_options.append("Monochrome")

        if dataset_option != 'TCGA':
            if 'Histology' in color_options:
                default_index = color_options.index('Histology')
            elif 'NAPY' in color_options:
                default_index = color_options.index('NAPY')
            elif 'NE' in color_options:
                default_index = color_options.index('NE')
            elif 'EMT' in color_options:
                default_index = color_options.index('EMT')
            elif 'TNM' in color_options:
                default_index = color_options.index('TNM')
            elif 'Histopathology' in color_options:
                default_index = color_options.index('Histopathology')
            elif 'Molecular_clusters' in color_options:
                default_index = color_options.index('Molecular_clusters')
            elif 'UICC Stage' in color_options:
                default_index = color_options.index('UICC Stage')
            elif 'Tumor_Stage' in color_options:
                default_index = color_options.index('Tumor_Stage')
        else:
            if 'ST-Selected' in color_options:
                default_index = color_options.index('ST-Selected')
            elif 'ST-Other' in color_options:
                default_index = color_options.index('ST-Other')
            elif 'ST-other' in color_options:
                default_index = color_options.index('ST-other')
            elif 'Hist.Grade' in color_options:
                default_index = color_options.index('Hist.Grade')
            elif 'Hist.Type' in color_options:
                default_index = color_options.index('Hist.Type')
            elif 'Clinical Stage' in color_options:
                default_index = color_options.index('Clinical Stage')
            elif 'AJCC Stage' in color_options:
                default_index = color_options.index('AJCC Stage')
            elif 'ST-mRNA' in color_options:
                default_index = color_options.index('ST-mRNA')
            elif 'ST-protein' in color_options:
                default_index = color_options.index('ST-protein')
            else:
                default_index = color_options.index('Gender')

        gc_dot_color = advanced_set_col[1].selectbox(":blue[Select the ] :red[Coloring Palette]",
                                                     options=color_options, index=default_index, key='gc_dot_color')
        gc_mc_cp = None

        if gc_dot_color == 'Monochrome':
            gc_mc_cp = advanced_set_col[1].color_picker(":blue[Pick the ] :red[Color of Dots ] :blue[in the Plot]",
                                                        value='#00AAF0',
                                                        key='gc_mc_cp')



        # Ensure that the indices align and only the common samples are used for correlation
        common_samples = df.loc[selected_gene_x].index.intersection(df.loc[selected_gene_y].index)

        df_x_gene = df.loc[selected_gene_x, common_samples]
        df_y_gene = df.loc[selected_gene_y, common_samples]

        corr, p_value = gc_corr(df_x_gene, df_y_gene, gc_corr_meth, selected_gene_y)

        st.write(f':blue[{gc_corr_meth} Correlation Coefficient: ] :red[{corr}]')
        st.write(f':blue[{gc_corr_meth} P Value: ] :red[{p_value}]')

        fig_scatter = plot_gene_correlation(df_x_gene, df_y_gene, selected_gene_x, selected_gene_y, corr, p_value,
                                            gc_dot_color, gc_mc_cp, merged_df, comp_type="gc", user_id=user_id)

        order_col = st.columns(2)
        fig_violin_y, p_value_violin_y = plot_custom_fusion(df_y_gene, selected_gene_y, gc_dot_color, gc_mc_cp,
                                                            merged_df, comp_type="gc", user_id=user_id,
                                                            axis_id="y_gene", key='y_gene_gc', multi_col=order_col[1])
        fig_violin_x, p_value_violin_x = plot_custom_fusion(df_x_gene, selected_gene_x, gc_dot_color, gc_mc_cp,
                                                            merged_df, comp_type="gc", user_id=user_id,
                                                            axis_id="x_gene", key='x_gene_gc', multi_col=order_col[0])

        # Show plot
        if corr_comp_chb:
            plot_col = st.columns(2)
            plot_col[1].subheader(":blue[Top 1000 Lowest P Value]", divider="blue")
            plot_col[1].dataframe(corr_comp_df.head(999), use_container_width=True, hide_index=False)
            file_path_ssgsea_p = naming(user_id)[37]
            corr_comp_df.to_csv(file_path_ssgsea_p)
            plot_col[0].plotly_chart(fig_scatter, use_container_width=False)
        else:
            plot_col_no_compare = st.columns(3)
            plot_col_no_compare[1].plotly_chart(fig_scatter, use_container_width=False)

        plot_col_2 = st.columns(2)
        plot_col_2[0].plotly_chart(fig_violin_x, use_container_width=False)
        plot_col_2[1].plotly_chart(fig_violin_y, use_container_width=False)

        if gc_dot_color != 'Monochrome':
            plot_col_2[0].dataframe(p_value_violin_x, use_container_width=True)
            plot_col_2[1].dataframe(p_value_violin_y, use_container_width=True)

        st.toast(":red[Gene Compare ] :blue[is Completed!]", icon='🎉')

    except Exception as e:
        st.warning(f":blue[SVD did ] :red[NOT ] :blue[converge in Linear Least Squares. ] :blue[Please choose ] :red[another gene!] {e}")

    file_copy(user_id)


def gc_corr(x_gene, y_gene, corr_meth, auto):
    """
    Calculates the correlation coefficient and p-value between two genes using the specified
    correlation method.

    Parameters:
    - x_gene (Series): Expression data for the first gene.
    - y_gene (Series): Expression data for the second gene.
    - corr_meth (str): Method of correlation ('Spearman' or 'Pearson').
    - auto (str): Automatic setting indicator (currently unused).

    Returns:
    - corr (float): Calculated correlation coefficient.
    - p_value (float): P-value for the correlation.
    """
    try:
        if auto != 'Automatic':
            if corr_meth == "Spearman":
                corr, p_value = stats.spearmanr(x_gene, y_gene)
            else:
                corr, p_value = stats.pearsonr(x_gene, y_gene)
        else:
            p_value = 'NaN'
            corr = 'NaN'
        return corr, p_value
    except Exception as e:
        st.error(f"Error calculating correlation: {e}")
        return None, None


@st.cache_data(show_spinner="Correlation Calculating...", max_entries=10, ttl=1800)
def all_vs_one_correlation(df, x_gene, corr_meth):
    """
    Calculates the correlation and p-value between a selected gene and all other genes in the dataset,
    using the specified correlation method.

    Parameters:
    - df (DataFrame): DataFrame containing gene expression data.
    - x_gene (str): Expression data for the reference gene.
    - corr_meth (str): Method of correlation ('Spearman' or 'Pearson').

    Returns:
    - results_df (DataFrame): DataFrame containing correlation results and corrected p-values for all gene comparisons.
    """
    try:
        x_gene_expression = df.loc[x_gene]
        results = []

        for gene in df.index:
            if gene != x_gene:
                common_samples = df.columns.intersection(df.columns)
                y_gene_expression = df.loc[gene, common_samples]
                x_gene_common_expression = x_gene_expression[common_samples]
                corr, p_value = gc_corr_for_all(x_gene_common_expression, y_gene_expression, corr_meth)
                results.append({'Gene': gene, 'Correlation': corr, 'P-value': p_value})

        results_df = pd.DataFrame(results)
        results_df.sort_values(by=['P-value'], inplace=True)

        # Calculating FDR
        results_df.dropna(how='any', inplace=True)
        p_values = results_df['P-value'].values
        fdr_corrected_p_values = multipletests(p_values, method='fdr_bh')[1]
        results_df['FDR'] = fdr_corrected_p_values

        # Set pandas option to format numbers in scientific notation
        results_df['FDR'] = results_df['FDR'].apply(lambda x: f'{x:.2e}')
        results_df['P-value'] = results_df['P-value'].apply(lambda x: f'{x:.2e}')

        return results_df
    except Exception as e:
        st.error(f"Error calculating all vs one correlation: {e}")
        return pd.DataFrame()


# Define a helper function for calculating correlation
def gc_corr_for_all(x, y, method):
    """
    Calculates the correlation coefficient and p-value between two genes using the specified
    correlation method.

    Parameters:
    - x (Series): Expression data for the first gene.
    - y (Series): Expression data for the second gene.
    - method (str): Method of correlation ('Spearman' or 'Pearson').

    Returns:
    - corr (float): Calculated correlation coefficient.
    - p_value (float): P-value for the correlation.
    """
    try:
        if method == "Spearman":
            corr, p_value = spearmanr(x, y)
        else:  # Pearson
            corr, p_value = pearsonr(x, y)
        return corr, p_value
    except Exception as e:
        st.error(f"Error calculating correlation for all: {e}")
        return None, None
