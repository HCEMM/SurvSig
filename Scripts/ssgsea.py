import pandas as pd
from gseapy import ssgsea
import streamlit as st
import numpy as np
from Scripts.multivariate_analysis import sg_ssgsea_multi_data
from Scripts.rpy2_heatmap_plots import py_hm_ssgsea_surv
from Scripts.user_interaction import get_session_id, naming
from lifelines.statistics import logrank_test
from lifelines import CoxPHFitter
from statsmodels.stats.multitest import multipletests
from Scripts.rpy2_heatmap_plots import surv_plot
from scipy.stats import zscore
from Scripts.rpy2_heatmap_plots import multi_var
from Scripts.plotly_plots import plot_custom_strip, dot_scatter_ssgsea, fdr_p_plot
from PIL import Image
from Scripts.gene_compare import plot_gene_correlation, plot_custom_fusion, all_vs_one_correlation, gc_corr
import Scripts.descriptions as dsc
from Scripts.r_code_settings import file_copy
from Scripts.survival_info import survival_info
import concurrent.futures
from pathlib import Path


def ssgsea_surv(df_melt, clinic_df, dataset_option, anno_file, gene_cluster, gene_clust):
    """
    Perform survival analysis based on ssGSEA results.

    Args:
        df_melt (pd.DataFrame): Melted DataFrame containing ssGSEA results.
        clinic_df (pd.DataFrame): Clinical data DataFrame.
        dataset_option (str): Dataset option selected by the user.
        anno_file (str): Annotation file path.
        gene_cluster (bool): Whether to use gene clustering.
        gene_clust (pd.DataFrame): DataFrame containing gene clusters.

    Returns:
        None
    """
    step_size = None
    surv_df = None
    results_df = None
    th_num1 = None
    th_num2 = None
    th_num = None
    user_id = get_session_id()
    auto_correct = False

    st.subheader(" ", divider="blue")
    if gene_cluster:
        clust_opt = st.selectbox(":blue[Select a ] :red[cluster]", key="clust_opt_surv",
                                 options=list(df_melt['Cluster'].unique()), index=0)
    else:
        clust_opt = "All Genes"

    user_dir = f"/result_data_{user_id}"
    multi_db_name = f"result_data{user_dir}/multi_db_ssgsea_{user_id}.csv"
    multi_info_name = f"result_data{user_dir}/multi_var_info_ssgsea_{user_id}.tsv"
    forest_pdf = f"result_data{user_dir}/forest_ssgsea_gene_{user_id}.pdf"
    forest_png = f"result_data{user_dir}/forest_ssgsea_gene_{user_id}.png"
    surv_colors_path = f"result_data{user_dir}/surv_ssgsea_gene_colors{user_id}.csv"
    surv_dataframe_path = f"result_data{user_dir}/surv_ssgsea{user_id}.csv"
    surv_pdf = f"result_data{user_dir}/surv_plot_ssgsea{user_id}.pdf"
    surv_png = f"result_data{user_dir}/surv_plot_ssgsea{user_id}.png"
    pv_dot_name = f"result_data{user_dir}/p_value_ssgsea_{user_id}.pdf"

    df_melt = df_melt[df_melt['Cluster'] == clust_opt]

    mv_db1 = pd.DataFrame()

    if dataset_option == 'George-SCLC':
        clinic_df.index.rename("sample", inplace=True)
        clinic_df = pd.DataFrame(clinic_df)
        clinic_df.reset_index(inplace=True)
    elif dataset_option == 'Jiang-SCLC':
        clinic_df.index.rename("sample", inplace=True)
        clinic_df = pd.DataFrame(clinic_df)
        clinic_df.reset_index(inplace=True)
    elif dataset_option == 'George-LCNEC':
        clinic_df.index.rename("sample", inplace=True)
        clinic_df = pd.DataFrame(clinic_df)
        clinic_df.reset_index(inplace=True)
    elif dataset_option == 'Fernandez-Carcinoid':
        clinic_df.index.rename("sample", inplace=True)
        clinic_df = pd.DataFrame(clinic_df)
        clinic_df.reset_index(inplace=True)
    elif dataset_option == 'Liu-SCLC':
        clinic_df.index.rename("sample", inplace=True)
        clinic_df = pd.DataFrame(clinic_df)
        clinic_df.reset_index(inplace=True)
    elif dataset_option == 'Alcala-Carcinoid':
        clinic_df.index.rename("sample", inplace=True)
        clinic_df = pd.DataFrame(clinic_df)
        clinic_df.reset_index(inplace=True)
    elif dataset_option == 'Anish-SCLC':
        clinic_df.index.rename("sample", inplace=True)
        clinic_df = pd.DataFrame(clinic_df)
        clinic_df.reset_index(inplace=True)
    elif dataset_option == 'Rousseaux-Mixed(non-NE)':
        clinic_df.index.rename("sample", inplace=True)
        clinic_df = pd.DataFrame(clinic_df)
        clinic_df.reset_index(inplace=True)

    elif dataset_option == 'TCGA':
        clinic_df.rename(columns={"bcr_patient_barcode": "sample"}, inplace=True)
        clinic_df.dropna(subset=['Status', 'Time'])
        clinic_df["Time"] = clinic_df["Time"].astype('float')
        clinic_df["Status"] = clinic_df["Status"].astype('float')
        if "sample" in clinic_df.columns:
            clinic_df.set_index('sample', inplace=True)
        clinic_df = pd.DataFrame(clinic_df)
        clinic_df.reset_index(inplace=True)
    else:
        pass

    merged = pd.merge(df_melt, clinic_df, on="sample", how="inner")

    if merged.shape[0] > 1:

        separation_opt = ['Median', 'Percentage decomposition', 'Percentage decomposition (lower and upper limit)',
                          'Automatic', 'ssGSEA values cutoff']
        sep_opt = st.selectbox(":blue[Please Select the ] :red[Method of Sample Separation]", options=separation_opt,
                               key='sep_opt', index=3, help=dsc.sep_method())

        ssgsea_form = st.form("ssgsea_form")

        if sep_opt == "Median":
            median_value = merged['ssGSEA'].median()
            st.write(f":blue[Median value: ] :red[{median_value}]")
            merged['Group'] = np.where(merged['ssGSEA'] <= median_value, "Low", "High")
            surv_df = merged[['sample', 'Time', 'Status', 'Group']]
            surv_df.rename(columns={'Group': 'Cluster'}, inplace=True)
            surv_df.set_index('sample', inplace=True)
            th_num = median_value

        elif sep_opt == 'ssGSEA values cutoff':
            threshold = ssgsea_form.slider(":blue[Set the ] :red[threshold]", max_value=merged['ssGSEA'].max() - 0.01,
                                           min_value=merged['ssGSEA'].min() + 0.01,
                                           value=merged['ssGSEA'].median(), key='ssgsea_th')
            merged['Group'] = np.where(merged['ssGSEA'] < threshold, "Low", "High")
            st.write(f":blue[Threshold value: ] :red[{threshold}]")
            surv_df = merged[['sample', 'Time', 'Status', 'Group']]
            surv_df.rename(columns={'Group': 'Cluster'}, inplace=True)
            surv_df.set_index('sample', inplace=True)
            th_num = threshold

        elif sep_opt == 'Percentage decomposition':
            threshold_per = ssgsea_form.slider(":blue[Set the ] :red[threshold]", max_value=90, min_value=10, value=50,
                                               key='ssgsea_per')
            threshold = merged['ssGSEA'].quantile(threshold_per / 100)
            merged['Group'] = np.where(merged['ssGSEA'] < threshold, "Low", "High")
            st.write(f":blue[Threshold value: ] :red[{threshold}]")
            if int(merged['Group'].nunique()) < 2:
                st.error(':blue[You have only ] :red[ONE ] :blue[group. Please choose ] :red[lower/higher limit!] ')
            else:
                surv_df = merged[['sample', 'Time', 'Status', 'Group']]
                surv_df.rename(columns={'Group': 'Cluster'}, inplace=True)
                surv_df.set_index('sample', inplace=True)
            th_num = threshold

        elif sep_opt == 'Percentage decomposition (lower and upper limit)':
            threshold_per = ssgsea_form.slider(":blue[Set the ] :red[threshold]", 10, 90, (33, 66),
                                               key="perc_ssgsea_up_down")

            lower_threshold = merged['ssGSEA'].quantile(threshold_per[0] / 100)
            upper_threshold = merged['ssGSEA'].quantile(threshold_per[1] / 100)
            perc_cut_low_hight = st.columns(5)
            perc_cut_low_hight[0].write(f":blue[Lower threshold: ] :red[{lower_threshold}]")
            perc_cut_low_hight[4].write(f":blue[Upper threshold: ] :red[{upper_threshold}]")

            conditions = [
                (merged['ssGSEA'] < lower_threshold),
                (merged['ssGSEA'] >= lower_threshold) & (merged['ssGSEA'] <= upper_threshold),
                (merged['ssGSEA'] > upper_threshold)
            ]
            values = ["Low", "Intermediate", "High"]

            merged['Group'] = np.select(conditions, values)

            th_num1 = lower_threshold
            th_num2 = upper_threshold

            if int(merged['Group'].nunique()) < 3:
                st.error(':blue[You have only ] :red[ONE ] :blue[group. Please choose ] :red[lower/higher limit4]')
            else:
                surv_df = merged[['sample', 'Time', 'Status', 'Group']]
                surv_df.rename(columns={'Group': 'Cluster'}, inplace=True)
                surv_df.set_index('sample', inplace=True)

        elif sep_opt == 'Automatic':
            try:
                surv_df = clinic_df[['sample', 'Time', 'Status']]
                surv_df.set_index('sample', inplace=True)
                auto_col = ssgsea_form.columns(2)
                threshold_interval = auto_col[0].slider(":blue[Set the ] :red[Interval]", 5, 95, (10, 90),
                                                        key="perc_automatic_ssgsea", help=dsc.perc_int())
                step_size = int(auto_col[1].select_slider(":blue[Set the ] :red[Step]", options=['1', '5', '10'], value='1',
                                                          key="step_ssgsea", help=dsc.sep_step()))

                data_list = []

                min_p_value = float('inf')
                optimal_threshold = None
                for threshold_per in range(threshold_interval[0], threshold_interval[1] + 1, step_size):
                    threshold = merged['ssGSEA'].quantile(threshold_per / 100)
                    merged['Group'] = np.where(merged['ssGSEA'] < threshold, "Low", "High")

                    merged['High_Group'] = np.where(merged['Group'] == "High", 1, 0)
                    group1_time = merged[merged['Group'] == "Low"]['Time']
                    group2_time = merged[merged['Group'] == "High"]['Time']
                    event_observed_A = merged[merged['Group'] == "Low"]['Status']
                    event_observed_B = merged[merged['Group'] == "High"]['Status']

                    result = logrank_test(group1_time, group2_time, event_observed_A=event_observed_A,
                                          event_observed_B=event_observed_B)
                    p_value = result.p_value
                    cox_df = merged[['Time', 'Status', 'High_Group']].dropna(subset=['Time', 'Status',
                                                                                     'High_Group'])

                    if not cox_df.empty:
                        cph = CoxPHFitter(penalizer=0.0001)
                        cph.fit(cox_df, 'Time', event_col='Status')

                        if 'High_Group' in cph.summary.index:
                            hazard_ratio_estimate = cph.summary.loc['High_Group', 'exp(coef)']
                        else:
                            hazard_ratio_estimate = None
                    else:
                        hazard_ratio_estimate = None

                    data_list.append({
                        'Threshold_Percentage': threshold_per,
                        'Threshold_Value': threshold,
                        'P_Value': p_value,
                        'HR': hazard_ratio_estimate,
                        'HR_reciprocal': (1 / hazard_ratio_estimate)
                    })

                    if p_value < min_p_value:
                        min_p_value = p_value
                        optimal_threshold = threshold

                results_df = pd.DataFrame(data_list)

                pvals = results_df['P_Value'].values
                fdr_corrected_pvals = multipletests(pvals, method='fdr_bh')[1]

                results_df['P_Value_FDR'] = fdr_corrected_pvals
                results_df.set_index("Threshold_Percentage", inplace=True)
                st.write(f":blue[Optimal threshold: ] :red[{optimal_threshold}]")

                merged['Group'] = np.where(merged['ssGSEA'] < optimal_threshold, "Low", "High")
                surv_df = merged[['sample', 'Time', 'Status', 'Group']]
                surv_df.rename(columns={'Group': 'Cluster'}, inplace=True)
                surv_df.set_index('sample', inplace=True)
                th_num = optimal_threshold
                auto_correct = True

            except Exception as e:

                st.warning(f":blue[Convergence halted. ] :red[Attempting to set threshold automatically.] "
    
                           f":blue[Error Details: {e}]")
                # Initialize list of alternative percentiles to try
                alternative_percentiles = [0.5, 0.25, 0.75]  # median, 25th, 75th
                threshold_set = False
                for perc in alternative_percentiles:
                    try:
                        threshold = merged['ssGSEA'].quantile(perc)
                        merged['Group'] = np.where(merged['ssGSEA'] <= threshold, "Low", "High")
                        unique_groups = merged['Group'].nunique()
                        if unique_groups >= 2:
                            st.write(f":blue[Threshold set to {int(perc * 100)}th percentile: ] :red[{threshold}]")
                            th_num = threshold
                            threshold_set = True
                            surv_df = merged[['sample', 'Time', 'Status', 'Group']]
                            surv_df.rename(columns={'Group': 'Cluster'}, inplace=True)
                            surv_df.set_index('sample', inplace=True)
                            st.info(
                                f":blue[The threshold has been automatically set to the {int(perc * 100)}th percentile to create two groups.]")
                            break
                        else:
                            st.warning(
                                f":blue[Threshold at {int(perc * 100)}th percentile ({threshold}) did not create two groups.]")
                    except Exception as inner_e:
                        st.warning(
                            f":blue[Failed to set threshold at {int(perc * 100)}th percentile due to error: {inner_e}]")

                if not threshold_set:
                    st.error(":red[Unable to separate samples into two groups. Please adjust your data or settings.]")

                    st.stop()

        plot_ext = ssgsea_form.expander(":blue[Advanced Settings of Plots]", expanded=False)
        color_cols = plot_ext.columns(3)

        if sep_opt == 'Percentage decomposition (lower and upper limit)':
            color1 = color_cols[2].color_picker(':blue[Group ] :red[High]', value="#e41a1d", key='color1_ss')
            color2 = color_cols[0].color_picker(':blue[Group ] :red[Low]', value="#377eb8", key='color2_ss')
            color3 = color_cols[1].color_picker(':blue[Group ] :red[Intermediate]', value='#4daf4a', key='color3_ss')
            color_list_plotly = [color2, color3, color1]
            groups = surv_df["Cluster"].unique()
            surv_colors = pd.DataFrame()
            surv_colors["cluster"] = groups
            surv_colors["color"] = color_list_plotly
        else:
            color1 = color_cols[2].color_picker(':blue[Group ] :red[High]', value="#e41a1d", key='color1_b_ss')
            color2 = color_cols[0].color_picker(':blue[Group ] :red[Low]', value="#377eb8", key='color2_b_ss')
            color3 = None
            color_list_plotly = [color2, color1]
            groups = surv_df["Cluster"].unique()
            surv_colors = pd.DataFrame()
            surv_colors["cluster"] = groups
            surv_colors["color"] = color_list_plotly

        surv_df_clust_unique = surv_df["Cluster"].unique()
        surv_colors = surv_colors[surv_colors['cluster'].isin(groups)]

        st.subheader(" ", divider="blue")
        ssgsea_button = ssgsea_form.form_submit_button("Accept Settings and Start Plotting", help=dsc.act_btn())

        if ssgsea_button:
            ssgsea_form.info("The changes have been :red[SUCCESSFULLY] refreshed")

        mv_group_df = pd.DataFrame()
        mv_group_df.index = surv_df.index
        mv_group_df["Cluster"] = surv_df["Cluster"]

        advancet_cols = plot_ext.columns(2)
        min_num = advancet_cols[0].number_input(":blue[Select] :red[minimum number of groups]", min_value=1, max_value=10,
                                                value=3,
                                                key="min_num_ssgsea")

        mv_db1, cols_to_filter, available_columns = sg_ssgsea_multi_data(dataset_option, mv_group_df, surv_df, clinic_df, anno_file)

        invalid_values = {}
        for col in cols_to_filter:
            if col in mv_db1.columns and mv_db1[col].dtype == 'object':
                value_counts = mv_db1[col].value_counts()
                invalid_values[col] = value_counts[value_counts < min_num].index.tolist()
            if col in mv_db1.columns and mv_db1[col].dtype == 'int':
                value_counts = mv_db1[col].value_counts()
                invalid_values[col] = value_counts[value_counts < min_num].index.tolist()

        for col, values in invalid_values.items():
            mv_db1.loc[mv_db1[col].isin(values), col] = 'nan'

        columns_with_multiple_values = [
            col for col in mv_db1.columns
            if (mv_db1[col].nunique() >= 2) or
               (mv_db1[col].nunique() == 2 and not mv_db1[col].isin(['nan']).any())
        ]

        selected_features_option = list(set(columns_with_multiple_values) & set(available_columns))

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

        selected_features = advancet_cols[1].multiselect(":blue[Select the] :red[Features]",
                                                         options=selected_features_option,
                                                         key="multi_adv_g_sgssea", default=selected_features_option_default,
                                                         max_selections=5)

        if df_melt.nunique()[0] > 200:
            default_dendro_sample = False
        else:
            default_dendro_sample = True
        if gene_clust.shape[0] > 400:
            default_dendro_gene = False
        else:
            default_dendro_gene = True

        dendogram = advancet_cols[0].toggle(":blue[Use ] :red[Dendogram ] :blue[for Heatmap's Samples]",
                                              key="dendogram_sample", value=default_dendro_sample)
        dendogram_gene = advancet_cols[1].toggle(":blue[Use ] :red[Dendogram ] :blue[for Heatmap's Genes]",
                                                   key="dendogram_gene", value=default_dendro_gene)

        selectbox_df = pd.DataFrame(columns=['Selected References'])

        if len(selected_features) == 0:
            selected_features = ['Cluster']

        for feature in selected_features:
            unique_values = [value for value in mv_db1[feature].unique().tolist() if value != 'nan']
            most_common_value = mv_db1[feature].mode().dropna().iloc[0] if not mv_db1[
                feature].mode().dropna().empty else None

            if most_common_value is None:
                continue

            initial_index = unique_values.index(most_common_value) if most_common_value in unique_values else 0

            selected_val = plot_ext.selectbox(f":blue[Select] :red[reference] :blue[value for] :red[{feature}]",
                                              options=unique_values, key=f"select_ssgsea_{feature}",
                                              index=initial_index)

            selectbox_df.loc[feature] = selected_val

        if len(selected_features) == mv_db1.shape[1]:
            plot_ext.info(":red[All] features are selected")
        else:
            plot_ext.info(f":blue[Selected feature(s): ] :red[{str(selected_features)}]")

        if len(selected_features) != 0:
            selected_features.extend(["Status", "Time"])
            mv_db1 = mv_db1[selected_features]
        mv_db1.replace("NaN", "nan", inplace=True)
        mv_db1.to_csv(multi_db_name)

        if sep_opt == 'Percentage decomposition (lower and upper limit)':
            int_exclusion = plot_ext.toggle(":blue[Exclude ] :red[Intermediate ] :blue[Group]", key='int_exclusion')

            if int_exclusion:
                surv_df = surv_df[surv_df['Cluster'] != 'Intermediate']
                surv_colors = surv_colors[surv_colors['cluster'] != 'Intermediate']

        surv_df.to_csv(surv_dataframe_path)
        surv_colors.to_csv(surv_colors_path, index=False)

        user_id = get_session_id()
        file_name_list = naming(user_id)
        surv_plot(surv_dataframe_path, surv_png, surv_pdf, user_id, dataset_option, surv_colors_path)

        surv_col = st.columns(3)
        try:
            surv_col[0].subheader(":blue[ssGSEA Cluster SurvPlot]", divider="blue")
            surv_col[0].image(surv_png, use_container_width=True)
        except Exception as e:
            surv_col[0].image("style_items/error.svg", use_container_width=True)
            surv_col[0].warning(
                f":blue[There is problem with the ] :red[chart! ] :blue[Please try another ] :red[settings ] "
                f":blue[or ] :red[refresh ] :blue[page] Error: {e}")

        hm_ssgsea_groups = pd.DataFrame()
        hm_ssgsea_groups.index = surv_df.index
        hm_ssgsea_groups["Groups"] = surv_df["Cluster"]
        hm_ssgsea_groups_path = file_name_list[32]

        hm_ssgsea_groups.to_csv(hm_ssgsea_groups_path)

        merged.set_index("sample", inplace=True)
        selected_cluster = str(merged["Cluster"].unique()[0])
        ssgsea_dot_path = naming(user_id)[53]
        fig = dot_scatter_ssgsea(merged, color_list_plotly, th_num, th_num1, th_num2, sep_opt, ssgsea_dot_path)

        try:
            surv_col[1].subheader(f":blue[ssGSEA profile of Cluster {selected_cluster}]", divider="blue")
            surv_col[1].plotly_chart(fig, use_container_width=True)
        except Exception as e:
            surv_col[1].image("style_items/error.svg", use_container_width=True)
            surv_col[1].warning(
                f":blue[There is problem with the ] :red[chart! ] :blue[Please try another ] :red[settings ] "
                f":blue[or ] :red[refresh ] :blue[page] Error: {e}")

        anno_colors_path = file_name_list[42]  # 42
        anno_colors_path2 = file_name_list[55]  # 55
        ht_expression_path = file_name_list[54]  # 54
        ht_top_annotation_path = file_name_list[2]  # 2
        ht_top_annotation_path2 = file_name_list[57]  # 57
        ssgsea_row_cluster = file_name_list[58]  # 58
        ht_png = file_name_list[59]  # 59
        ht_pdf = file_name_list[60]  # 60
        ht_row_cluster_path = file_name_list[6]  # 6

        anno_colors = pd.read_csv(anno_colors_path)
        anno_colors = anno_colors[anno_colors['column'] != 'Cluster']

        surv_colors['column'] = 'Cluster'
        surv_colors.rename(columns={'cluster': 'type'}, inplace=True)
        surv_colors = surv_colors[['column', 'type', 'color']]

        result_df = pd.concat([anno_colors, surv_colors], ignore_index=True)

        result_df.to_csv(anno_colors_path2, index=False)

        top_annotation = pd.read_csv(ht_top_annotation_path)

        if dataset_option == "TCGA":
            top_annotation.rename(columns={"bcr_patient_barcode": "Samples"}, inplace=True)

        top_annotation.drop(['Cluster'], axis=1, inplace=True)
        top_annotation.set_index("Samples", inplace=True)
        top_annotation['Cluster'] = surv_df['Cluster']
        top_annotation.reset_index(inplace=True)

        top_annotation.to_csv(ht_top_annotation_path2, index=False)

        if gene_cluster:
            gene_df = pd.read_csv(ht_row_cluster_path)
        else:
            gene_df = pd.DataFrame()
            gene_df["Samples"] = gene_clust.index
            gene_df["Cluster"] = "All Genes"

        if clust_opt == "All Genes":
            gene_ssgsea_cluster = pd.DataFrame()
            gene_ssgsea_cluster["Names"] = gene_df["Samples"]
            gene_ssgsea_cluster["Cluster"] = "Genes"
            gene_ssgsea_cluster.to_csv(ssgsea_row_cluster, index=False)
        else:
            gene_ssgsea_cluster = pd.DataFrame()
            gene_ssgsea_cluster["Names"] = gene_df["Samples"]
            gene_ssgsea_cluster["Cluster"] = gene_df["Cluster"]
            if gene_ssgsea_cluster["Cluster"].dtype != "object":
                gene_ssgsea_cluster = gene_ssgsea_cluster[gene_ssgsea_cluster['Cluster'] == int(clust_opt)]
                gene_ssgsea_cluster.to_csv(ssgsea_row_cluster, index=False)
            else:
                gene_ssgsea_cluster = gene_ssgsea_cluster[gene_ssgsea_cluster['Cluster'] == clust_opt]
                gene_ssgsea_cluster.to_csv(ssgsea_row_cluster, index=False)

        py_hm_ssgsea_surv(ht_name="Expression",
                          ht_top_annotation_color=anno_colors_path2,
                          ht_expression_path=ht_expression_path, column_cluster_path=hm_ssgsea_groups_path,
                          ht_top_annotation_path=ht_top_annotation_path2,
                          ht_png=ht_png, ht_pdf=ht_pdf,
                          user_id=user_id,
                          ht_row_cluster_path=ssgsea_row_cluster, dataset_option=dataset_option,
                          dendrogram_display=dendogram, dendrogram_display_gene=dendogram_gene)

        try:
            surv_col[2].subheader(":blue[ssGSEA Cluster Heatmap]", divider="blue")
            surv_col[2].image(ht_png, use_container_width=True)
        except Exception as e:
            surv_col[0].image("style_items/error.svg", use_container_width=True)
            surv_col[0].warning(
                f":blue[There is problem with the ] :red[chart! ] :blue[Please try another ] :red[settings ] "
                f":blue[or ] :red[refresh ] :blue[page] Error: {e}")

        plotly_col = st.columns(3)
        try:
            multi_var(user_id, selectbox_df, multi_info_name, forest_pdf, forest_png, multi_db_name)
            image = Image.open(forest_png)
            if sep_opt == "Automatic":
                plotly_col[0].subheader(":blue[Multivariate Plot]", divider="blue")
                plotly_col[0].image(image, use_container_width=True, output_format='PNG')
            else:
                multi_plot_col = st.columns(3)
                multi_plot_col[1].subheader(":blue[Multivariate Plot]", divider="blue")
                multi_plot_col[1].image(image, use_container_width=True, output_format='PNG')
        except Exception as e:
            if sep_opt == "Automatic":
                plotly_col[0].subheader(":blue[Multivariate Plot]", divider="blue")
                plotly_col[0].image("style_items/error.svg", use_container_width=True)
                plotly_col[0].warning(
                    f":blue[You have ] :red[infinity ] :blue[value(s). Try to remove] :red[features] "
                    f":blue[ or try to use other] :red[ clustering / dimensionality method]! Error: {e}")
            else:
                multi_plot_col = st.columns(3)
                multi_plot_col[1].subheader(":blue[Multivariate Plot]", divider="blue")
                multi_plot_col[1].image("style_items/error.svg", use_column_width=True)
                multi_plot_col[1].warning(
                    f":blue[You have ] :red[infinity ] :blue[value(s). Try to remove] :red[features] "
                    f":blue[ or try to use other] :red[ clustering / dimensionality method]! Error: {e}")

        if auto_correct:
            if sep_opt == 'Automatic':
                plotly_col[1].subheader(":blue[P-Values Comparison Table]", divider="blue")
                plotly_col[1].dataframe(results_df, use_container_width=True)

                lower_limit = results_df.index.min()
                upper_limit = results_df.index.max()

                results_df['hover_text_pv'] = results_df.apply(
                    lambda x: f"Threshold: {x['Threshold_Value']}<br>P_Value: {x['P_Value']}", axis=1)

                results_df['hover_text_pvfdr'] = results_df.apply(
                    lambda x: f"Threshold: {x['Threshold_Value']}<br>P_Value_FDR: {x['P_Value_FDR']}", axis=1)

                fig_p_value = fdr_p_plot(results_df, lower_limit, upper_limit, pv_dot_name)

                try:
                    plotly_col[2].subheader(":blue[P-Values Comparison]", divider="blue")
                    plotly_col[2].plotly_chart(fig_p_value, use_container_width=True)
                except Exception as e:
                    plotly_col[2].image("style_items/error.svg", use_container_width=True)
                    plotly_col[2].warning(
                        f":blue[There is problem with the ] :red[chart! ] :blue[Please try another ] :red[settings ] "
                        f":blue[or ] :red[refresh ] :blue[page] Error: {e}")

                st.toast(":red[ssGSEA Survival Analysis ] :blue[is Completed!]", icon='🎉')

                survival_info(surv_dataframe_path, user_id, surv_type="ssgsea")

                return gene_ssgsea_cluster
            else:
                st.warning(":blue[Please click on the ] :red['Accept Changes and Start Plotting' button]")

        file_copy(user_id)

    else:
        st.warning(":blue[Insufficient data available for] :red[Survival Analysis].")


@st.fragment
def ssgsea_analysis(seed, exp_data, gene_clust, dataset_option, clinic_df, strip_path,
                    ssgsea_result_path, user_id, anno_file, gene_cluster, subset_gene_exp_df,
                    st_anno, sample_cluster):
    """
    Perform ssGSEA analysis and generate visualizations.

    Args:
        seed (int): Random seed for reproducibility.
        exp_data (pd.DataFrame): Expression data DataFrame.
        gene_clust (pd.DataFrame): DataFrame containing gene clusters.
        dataset_option (str): Dataset option selected by the user.
        clinic_df (pd.DataFrame): Clinical data DataFrame.
        strip_path (str): Path to save the strip plot image.
        ssgsea_result_path (str): Path to save the ssGSEA results.
        user_id (str): User session ID.
        anno_file (pd.DataFrame): Annotation file DataFrame.
        gene_cluster (bool): Whether to use gene clustering.
        subset_gene_exp_df (pd.DataFrame): Subset of gene expression data.
        st_anno (pd.DataFrame): Streamlit annotation DataFrame.

    Returns:
        None
    """
    ssgsea_col = st.columns(2)
    cluster_select_container = st.container()
    st.subheader(" ", divider="blue")

    exp_data_hm = exp_data[exp_data.index.isin(subset_gene_exp_df.index)]
    z_score_toggle = ssgsea_col[0].selectbox(":blue[Heatmap ] :red[Visualisation:]", key="ztog_ssgsea",
                                             options=["Log2 Transformed", "Z-Scored"], index=0)
    if z_score_toggle == "Z-Scored":
        df_z_sg = zscore(exp_data_hm, axis=1, nan_policy='omit')
        df_z_sg = df_z_sg.dropna(axis=0)
    else:
        df_z_sg = exp_data_hm.copy()
    exp_data_hm = df_z_sg

    ht_expression_path = naming(user_id)[54]
    exp_data_hm.to_csv(ht_expression_path)

    unique_clusters = sorted(gene_clust['Cluster'].unique())
    unique_clusters.append('All Genes')

    if gene_cluster:
        ssgsea_clust = cluster_select_container.multiselect(":blue[Select ] :red[Cluster]", options=unique_clusters,
                                                            default=unique_clusters, key="ssgsea_ms")
        if ssgsea_clust:
            ssgsea_opt = ssgsea_clust
        else:
            ssgsea_opt = ["All Genes"]

        def Merge(dict1, dict2):
            for i in dict2.keys():
                dict1[i] = dict2[i]
            return dict1

        # Build the gene_sets based on the selected clusters
        gene_sets_cluster = {str(cluster): gene_clust[gene_clust['Cluster'] == cluster].index.astype(str).tolist() for
                             cluster in ssgsea_opt}

        if 'All Genes' in unique_clusters:
            all_gene_set = {'All Genes': gene_clust.index.tolist()}
            gene_sets_clusters = Merge(gene_sets_cluster, all_gene_set)

    else:
        ssgsea_opt = unique_clusters
        gene_sets_clusters = {str(cluster): gene_clust[gene_clust['Cluster'] == cluster].index.astype(str).tolist() for
                              cluster in ssgsea_opt}

    gsea_res_clusters = start_ssgsea(exp_data, gene_sets_clusters, seed)

    # Pivot the dataframe to make clusters as columns and samples as rows
    df_pivot = gsea_res_clusters.res2d.pivot(index='Name', columns='Term', values='NES')

    # Melt the dataframe for visualization
    df_melt = df_pivot.reset_index().melt(id_vars=["Name"], var_name="Cluster", value_name="ssGSEA")
    df_melt.rename(columns={'Name': 'sample'}, inplace=True)
    df_melt.sort_values(by=['ssGSEA'], inplace=True)

    df_melt['Cluster'] = df_melt['Cluster'].astype(str)
    df_melt_display = pd.DataFrame(df_melt)
    df_melt_display.rename(columns={"sample": "Sample ID"}, inplace=True)
    df_melt_display.set_index("Sample ID", inplace=True)

    selected_clusters = [str(c) for c in ssgsea_opt]
    df_melt = df_melt[df_melt['Cluster'].isin(selected_clusters)]

    fig = plot_custom_strip(df_melt, strip_path, anno_file, st_anno, ssgsea_col[1], key="ssgsea_analisys",
                            sample_cluster=sample_cluster)

    try:
        st.plotly_chart(fig, use_container_width=True, key="ssgsea_analisys")
    except:
        st.image("style_items/error.svg", use_container_width=True)
        st.warning(
            ":blue[There is problem with the ] :red[chart! ] :blue[Please try another ] :red[settings ] "
            ":blue[or ] :red[refresh ] :blue[page]")

    ssgsea_surv(df_melt, clinic_df, dataset_option, anno_file, gene_cluster, gene_clust)

    file_copy(user_id)

    st.header(":blue[ssGSEA Results]", divider='blue')
    user_dir = f"/result_data_{user_id}"
    file_path_ssgsea_p = f"result_data{user_dir}/ssgsea_melt_df_surv{user_id}.csv"
    df_melt_display.to_csv(file_path_ssgsea_p)
    st.dataframe(df_melt_display, use_container_width=True)


    file_copy(user_id)


@st.cache_resource(show_spinner="ssGSEA Calculating...", max_entries=10, ttl=1800)
def start_ssgsea(exp_data, gene_sets_clusters, seed):
    """
    Run ssGSEA for given clusters.

    Args:
        exp_data (pd.DataFrame): Expression data DataFrame.
        gene_sets_clusters (dict): Dictionary of gene sets for ssGSEA.
        seed (int): Random seed for reproducibility.

    Returns:
        pd.DataFrame: ssGSEA results DataFrame.
    """
    gsea_res_clusters = ssgsea(data=exp_data, gene_sets=gene_sets_clusters, min_size=1,
                               sample_norm_method=None, max_size=len(exp_data.index) - 1,
                               format='dataframe', outdir=None, seed=seed)
    return gsea_res_clusters

@st.fragment
def ssgsea_compare(df, strip_path, seed, exp_data, gene_clust, dataset_option, clinic_df, ssgsea_result_path, user_id,
                   anno_file, gene_cluster, st_anno, sample_cluster):
    """
    Compare ssGSEA results across different clusters and generate visualizations.

    Args:
        df (pd.DataFrame): DataFrame for ssGSEA comparison.
        strip_path (str): Path to save the strip plot image.
        seed (int): Random seed for reproducibility.
        exp_data (pd.DataFrame): Expression data DataFrame.
        gene_clust (pd.DataFrame): DataFrame containing gene clusters.
        dataset_option (str): Dataset option selected by the user.
        clinic_df (pd.DataFrame): Clinical data DataFrame.
        ssgsea_result_path (str): Path to save the ssGSEA results.
        user_id (str): User session ID.
        anno_file (pd.DataFrame): Annotation file DataFrame.
        gene_cluster (bool): Whether to use gene clustering.
        st_anno (pd.DataFrame): Streamlit annotation DataFrame.

    Returns:
        None
    """
    cols_to_drop_anno = anno_file.columns[anno_file.nunique() < 2]
    cols_to_drop_st = st_anno.columns[st_anno.nunique() < 2]
    cols_to_drop_clin = clinic_df.columns[clinic_df.nunique() < 2]

    anno_file = anno_file.drop(columns=cols_to_drop_anno)
    st_anno = st_anno.drop(columns=cols_to_drop_st)
    clinic_df = clinic_df.drop(columns=cols_to_drop_clin)

    df = df.transpose()
    exp_data = exp_data.transpose()

    df = df[df.index.isin(exp_data.index)]

    df = df.transpose()
    exp_data = exp_data.transpose()

    ssgsea_col = st.columns(2)
    st.subheader(" ", divider="blue")

    unique_clusters = sorted(gene_clust['Cluster'].unique())
    unique_clusters.append('All Genes')

    if gene_cluster:
        ssgsea_clust = ssgsea_col[1].multiselect(":blue[Select ] :red[Cluster]", options=unique_clusters,
                                                            default=unique_clusters, key="ssgsea_ms_comp")
        if ssgsea_clust:
            ssgsea_opt = ssgsea_clust
        else:
            ssgsea_opt = ["All Genes"]

        def Merge(dict1, dict2):
            for i in dict2.keys():
                dict1[i] = dict2[i]
            return dict1

        # Build the gene_sets based on the selected clusters
        gene_sets_cluster = {str(cluster): gene_clust[gene_clust['Cluster'] == cluster].index.astype(str).tolist() for
                             cluster in ssgsea_opt}

        if 'All Genes' in unique_clusters:
            all_gene_set = {'All Genes': gene_clust.index.tolist()}
            gene_sets_clusters = Merge(gene_sets_cluster, all_gene_set)

    else:
        ssgsea_opt = unique_clusters
        gene_sets_clusters = {str(cluster): gene_clust[gene_clust['Cluster'] == cluster].index.astype(str).tolist() for
                              cluster in ssgsea_opt}

    gsea_res_clusters = start_ssgsea(exp_data, gene_sets_clusters, seed)

    # Pivot the dataframe to make clusters as columns and samples as rows
    df_pivot = gsea_res_clusters.res2d.pivot(index='Name', columns='Term', values='NES')

    # Melt the dataframe for visualization
    df_melt = df_pivot.reset_index().melt(id_vars=["Name"], var_name="Cluster", value_name="ssGSEA")
    df_melt.rename(columns={'Name': 'sample'}, inplace=True)
    df_melt.sort_values(by=['ssGSEA'], inplace=True)

    df_melt['Cluster'] = df_melt['Cluster'].astype(str)
    df_melt_display = pd.DataFrame(df_melt)
    df_melt_display.rename(columns={"sample": "Sample ID"}, inplace=True)
    df_melt_display.set_index("Sample ID", inplace=True)

    selected_clusters = [str(c) for c in ssgsea_opt]
    df_melt = df_melt[df_melt['Cluster'].isin(selected_clusters)]
    ssgsea_cluster_options = df_melt['Cluster'].unique()

    fig = plot_custom_strip(df_melt, strip_path, anno_file, st_anno, ssgsea_col[0], key="ssgsea_compare",
                            sample_cluster=sample_cluster)

    try:
        st.plotly_chart(fig, use_container_width=True, key="ssgsea_compare")
    except:
        st.image("style_items/error.svg", use_container_width=True)
        st.warning(
            ":blue[There is problem with the ] :red[chart! ] :blue[Please try another ] :red[settings ] "
            ":blue[or ] :red[refresh ] :blue[page]")

    st.header(":blue[ssGSEA Results]", divider='blue')
    user_dir = f"/result_data_{user_id}"
    file_path_ssgsea_p = f"result_data{user_dir}/ssgsea_melt_df{user_id}.csv"
    df_melt_display.to_csv(file_path_ssgsea_p)
    st.dataframe(df_melt_display, use_container_width=True)

    try:
        settings_col = st.columns(4)
        axis_gene_col = st.columns(2)
        # Identify overlapping columns
        corr_comp_df = None
        default_index = 0
        overlap_columns = anno_file.columns.intersection(clinic_df.columns)

        # Rename overlapping columns in clini_df before merging
        clinic_df_renamed = clinic_df.rename(columns={col: f"{col}_from_clinic_df" for col in overlap_columns})

        # Merge anno_file and clini_df_renamed
        merged_df = pd.merge(anno_file, clinic_df_renamed, left_index=True, right_index=True, how='outer')
        # Assuming df1 has the column 'sample_cluster' and df2 is 'merged_df'
        merged_df['Sample_Cluster'] = sample_cluster['Cluster'].reindex(merged_df.index)

        for col in overlap_columns:
            merged_df.drop(col, axis=1, inplace=True)  # Drop anno_file's original overlapping columns
            merged_df.rename(columns={f"{col}_from_clinic_df": col}, inplace=True)  # Rename clini_df's columns back

        # For selectbox
        # Ensure your DataFrame index is a list or convert it to a list
        df_index = df.index.tolist()

        # Find the index of 'TP53' in the list
        try:
            selected_index = df_index.index('TP53')
        except ValueError:
            selected_index = 110  # Default to the first item if TP53 is not found

        selected_gene_y = axis_gene_col[0].selectbox(":blue[Select a Y Axis ] :red[ssGSEA Cluster]",
                                                     options=ssgsea_cluster_options,
                                                     key="gene_selx_ssgsea")

        if len(ssgsea_cluster_options) >= 2:
            compare_type = settings_col[2].selectbox(":blue[Basis of ] :red[Comparison]",
                                                     options=["Gene", "ssGSEA Cluster"],
                                                     key="gene_vs_ssgsea")
        else:
            compare_type = "Gene"

        df_y_gene_tmp = df_melt_display[df_melt_display["Cluster"] == selected_gene_y]
        df_y_gene_tmp.drop("Cluster", axis=1, inplace=True)
        df_y_gene_tmp = df_y_gene_tmp.transpose()
        df_y_gene_tmp = df_y_gene_tmp.rename(index={'ssGSEA': selected_gene_y})
        df_y_gene_tmp = df_y_gene_tmp[df.columns]
        df_y_gene_tmp = df_y_gene_tmp.astype(float)
        df = pd.concat([df, df_y_gene_tmp])
        df = df.astype(float)
        df_melt_display.index.rename("index", inplace=True)

        if compare_type == "Gene":
            corr_comp_text = 'Genes'
        else:
            corr_comp_text = 'ssGSEA Cluster(s)'

        corr_comp_chb = settings_col[0].toggle(
            f":blue[Compare ] :red[All ] :blue[{corr_comp_text}: ] :red[{selected_gene_y} Cluster]",
            key='corr_comp_chb_ssgea')

        if compare_type != "Gene":
            df = df_melt_display
            df.reset_index(inplace=True)
            df = df.pivot_table(index='Cluster', columns='index', values='ssGSEA', aggfunc='first')
            df = df.astype(float)

        if corr_comp_chb:
            corr_meth = st.selectbox(":blue[Select the ] :red[Correlation Method]", options=["Spearman", "Pearson"],
                                     index=0, key='corr_meth_ssgsea')
            corr_comp_df = all_vs_one_correlation(df, selected_gene_y, corr_meth)
            corr_comp_df.set_index("Gene", inplace=True)

        if compare_type == "Gene":
            if corr_comp_chb:
                selected_gene_x = axis_gene_col[1].selectbox(":blue[Select a ] :red[X Axis Gene]",
                                                             options=corr_comp_df.index, index=0,
                                                             key="ssgsea_sely")
                axis_gene_col[1].caption(":blue[Check your ] :red[X Axis Gene]",
                                         help="https://www.genecards.org/Search/Keyword?queryString=" + selected_gene_x)
            else:
                selected_gene_x = axis_gene_col[1].selectbox(":blue[Select a ] :red[X Axis Gene]", options=df_index,
                                                             index=selected_index)
                axis_gene_col[1].caption(":blue[Check your ] :red[X Axis Gene]",
                                         help="https://www.genecards.org/Search/Keyword?queryString=" + selected_gene_x)
        else:
            if corr_comp_chb:
                selected_gene_x = axis_gene_col[1].selectbox(":blue[Select a X Axis ] :red[ssGSEA Cluster]",
                                                             options=corr_comp_df.index, index=0,
                                                             key="ssgsea_sely")
            else:
                x_ssgsea_cluster_options = ssgsea_cluster_options[ssgsea_cluster_options != selected_gene_y]
                selected_gene_x = axis_gene_col[1].selectbox(":blue[Select a X Axis ] :red[ssGSEA Cluster]",
                                                             options=x_ssgsea_cluster_options,
                                                             index=0)

        st.subheader(" ", divider="blue")

        df_x_gene = df.loc[[selected_gene_x]]
        df_y_gene = df.loc[[selected_gene_y]]

        advanced_set_col = st.columns(2)

        gc_corr_meth = advanced_set_col[0].selectbox(":blue[Select the ] :red[Correlation Method]",
                                                     options=["Spearman", "Pearson"],
                                                     index=0, key='ssgsea_corr_method')

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
                                                     options=color_options, index=default_index, key='ssgsea_dot_color')

        gc_mc_cp = None

        if gc_dot_color == 'Monochrome':
            gc_mc_cp = advanced_set_col[1].color_picker(":blue[Pick the ] :red[Color of Dots ] :blue[in the Plot]",
                                                        value='#00AAF0',
                                                        key='ssgsea_mc_cp')

        # Before calling gc_corr, ensure df_x_gene and df_y_gene are Series or 1-dimensional arrays
        df_x_gene = df.loc[selected_gene_x]  # Use single brackets to get a Series
        df_y_gene = df.loc[selected_gene_y]  # Use single brackets to get a Series

        # Ensure that the indices align and only the common samples are used for correlation
        common_samples = df_y_gene.index.intersection(df_x_gene.index)

        df_x_gene = df_x_gene[common_samples]
        df_y_gene = df_y_gene[common_samples]

        corr, p_value = gc_corr(df_y_gene, df_x_gene, gc_corr_meth, selected_gene_x)

        st.write(f':blue[{gc_corr_meth} Correlation Coefficient: ] :red[{corr}]')
        st.write(f':blue[{gc_corr_meth} P Value: ] :red[{p_value}]')

        if compare_type == 'Gene':
            plot_compare_type = 'gc_ssgsea'
        else:
            plot_compare_type = 'ssgsea'

        fig_scatter = plot_gene_correlation(df_x_gene, df_y_gene, selected_gene_x, selected_gene_y, corr, p_value,
                                            gc_dot_color, gc_mc_cp, merged_df, comp_type=plot_compare_type,
                                            user_id=user_id)

        order_col = st.columns(2)

        fig_violin_y, p_value_violin_y = plot_custom_fusion(df_y_gene, selected_gene_y, gc_dot_color, gc_mc_cp,
                                                            merged_df, comp_type='ssgsea', user_id=user_id,
                                                            axis_id="y_axis", key="y_axix_ssgsea", multi_col=order_col[1])
        fig_violin_x, p_value_violin_x = plot_custom_fusion(df_x_gene, selected_gene_x, gc_dot_color, gc_mc_cp,
                                                            merged_df, comp_type='ssgsea', user_id=user_id,
                                                            axis_id="x_axis", key="x_axix_ssgsea", multi_col=order_col[0])

        # Show plot
        if corr_comp_chb:
            plot_col = st.columns(2)
            if compare_type == 'Gene':
                plot_col[1].subheader(":blue[Top 1000 Lowest P Value]", divider="blue")
            else:
                plot_col[1].subheader(":blue[P Value Table]", divider="blue")

            file_path_ssgsea_p = f"result_data{user_dir}/1000_lowest_ssgsea{user_id}.csv"
            corr_comp_df.to_csv(file_path_ssgsea_p)
            plot_col[1].dataframe(corr_comp_df.head(999), use_container_width=True, hide_index=False)
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

        st.toast(":red[ssGSEA Compare ] :blue[is Completed!]", icon='🎉')

    except Exception as e:
        st.warning(f":blue[SVD did ] :red[NOT ] :blue[converge in Linear Least Squares. ] "
                   f":blue[Please choose ] :red[another gene!] Error: {e}")
        st.stop()

    file_copy(user_id)