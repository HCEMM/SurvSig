import streamlit as st
import pandas as pd
import numpy as np
from Scripts.user_interaction import naming
from Scripts.multivariate_analysis import sg_ssgsea_multi_data
from lifelines.statistics import logrank_test
from lifelines import CoxPHFitter
from statsmodels.stats.multitest import multipletests
from Scripts.rpy2_heatmap_plots import surv_plot
from Scripts.rpy2_heatmap_plots import multi_var
from PIL import Image
from pandas.api.types import CategoricalDtype
from Scripts.plotly_plots import dot_scatter_sg, fdr_p_plot
import Scripts.descriptions as dsc
from Scripts.survival_info import survival_info
import concurrent.futures
from pathlib import Path


@st.fragment
def sg_analysis_surv(user_id, df, clinic_df, dataset_option, anno_file):
    """
    This function performs survival analysis on single gene expression data.

    Parameters:
    - user_id: Unique identifier for the user session.
    - df: DataFrame containing gene expression data.
    - clinic_df: DataFrame containing clinical data.
    - dataset_option: Selected dataset option for analysis.
    - anno_file: Annotation file containing additional data. - TCGA

    Returns:
    - None. The results are displayed using Streamlit components.
    """
    if dataset_option == "TCGA":
        sample_limit = anno_file[anno_file.index.isin(df.columns)].shape[0]
    else:
        sample_limit = clinic_df[clinic_df.index.isin(df.columns)].shape[0]

    if sample_limit >= 8:

        th_num = None
        th_num1 = None
        th_num2 = None
        results_df = None
        step_size = None
        hazard_ratio_estimate = None
        auto_correct = False

        sg_cont = st.container()
        sg_col = sg_cont.columns(2)
        sg_form = st.form("sg_form")

        # File paths using f-strings
        multi_db_name = naming(user_id)[43]
        multi_info_name = naming(user_id)[44]
        forest_pdf = naming(user_id)[45]
        forest_png = naming(user_id)[46]
        surv_colors_path = naming(user_id)[47]
        surv_dataframe_path = naming(user_id)[48]
        surv_pdf = naming(user_id)[49]
        surv_png = naming(user_id)[50]
        sg_dot_name = naming(user_id)[51]
        pv_dot_name = naming(user_id)[52]

        # Adjusting clinical data based on the selected dataset
        if dataset_option == 'George-SCLC':
            clinic_df.index.rename("sample", inplace=True)
            clinic_df = pd.DataFrame(clinic_df)
            clinic_df.reset_index(inplace=True)
            clinic_df.rename(columns={"Overall Survival (Months)": "Time"}, inplace=True)
            clinic_df.rename(columns={"Overall Survival Status": "Status"}, inplace=True)
        elif dataset_option == 'Anish-SCLC':
            clinic_df.index.rename("sample", inplace=True)
            clinic_df = pd.DataFrame(clinic_df)
            clinic_df.reset_index(inplace=True)
            clinic_df.rename(columns={"Overall survival since diagnosis": "Time"}, inplace=True)
            clinic_df.rename(columns={"Mortality status": "Status"}, inplace=True)
        elif dataset_option == 'Jiang-SCLC':
            clinic_df.index.rename("sample", inplace=True)
            clinic_df = pd.DataFrame(clinic_df)
            clinic_df.reset_index(inplace=True)
            clinic_df.rename(columns={"Months": "Time"}, inplace=True)
            clinic_df.rename(columns={"Survival Status": "Status"}, inplace=True)
        elif dataset_option == 'George-LCNEC':
            clinic_df.index.rename("sample", inplace=True)
            clinic_df = pd.DataFrame(clinic_df)
            clinic_df.reset_index(inplace=True)
            clinic_df.rename(columns={"Months": "Time"}, inplace=True)
            clinic_df.rename(columns={"Survival Status": "Status"}, inplace=True)
        elif dataset_option == 'Fernandez-Carcinoid':
            clinic_df.index.rename("sample", inplace=True)
            clinic_df = pd.DataFrame(clinic_df)
            clinic_df.reset_index(inplace=True)
            clinic_df.rename(columns={"Months": "Time"}, inplace=True)
            clinic_df.rename(columns={"Survival Status": "Status"}, inplace=True)
        elif dataset_option == 'Liu-SCLC':
            clinic_df.index.rename("sample", inplace=True)
            clinic_df = pd.DataFrame(clinic_df)
            clinic_df.reset_index(inplace=True)
            clinic_df.rename(columns={"Months": "Time"}, inplace=True)
            clinic_df.rename(columns={"Survival Status": "Status"}, inplace=True)
        elif dataset_option == 'Alcala-Carcinoid':
            clinic_df.index.rename("sample", inplace=True)
            clinic_df = pd.DataFrame(clinic_df)
            clinic_df.reset_index(inplace=True)
            clinic_df.rename(columns={"Months": "Time"}, inplace=True)
            clinic_df.rename(columns={"Survival Status": "Status"}, inplace=True)
        elif dataset_option == 'Rousseaux-Mixed(non-NE)':
            clinic_df.index.rename("sample", inplace=True)
            clinic_df = pd.DataFrame(clinic_df)
            clinic_df.reset_index(inplace=True)
            clinic_df.rename(columns={"Months": "Time"}, inplace=True)
            clinic_df.rename(columns={"Survival Status": "Status"}, inplace=True)
        elif dataset_option == 'TCGA':
            clinic_df.rename(columns={"bcr_patient_barcode": "sample"}, inplace=True)
            clinic_df.dropna(subset=['Status', 'Time'])
            clinic_df["Time"] = clinic_df["Time"].astype('float')
            clinic_df["Status"] = clinic_df["Status"].astype('float')
            if 'sample' in clinic_df.columns:
                clinic_df.set_index('sample', inplace=True)
            clinic_df = pd.DataFrame(clinic_df)
            clinic_df.reset_index(inplace=True)
        else:
            pass


        # For selectbox

        # Ensure your DataFrame index is a list or convert it to a list
        df_index = df.index.tolist()

        # Find the index of 'TP53' in the list
        try:
            selected_index = df_index.index('TP53')
        except ValueError:
            selected_index = 0  # Default to the first item if TP53 is not found

        selected_gene = sg_col[0].selectbox(":blue[Select a ] :red[Gene]", options=df_index, index=selected_index,
                                            key="gene_sel")
        df.loc[selected_gene, ].values.min()

        all_same = df.loc[selected_gene, ].values.min() == df.loc[selected_gene, ].values.max()

        if all_same:
            st.warning("blue[This gene has ] :red[NO varience. ] :blue[Please select another gene]")
            st.stop()

        median_error = df.loc[selected_gene, ].median() == df.loc[selected_gene, ].values.max()
        if median_error:
            st.warning("blue[This gene has ] :red[LOW varience. ] :blue[Please select another gene]")
            st.stop()

        sg_col[0].caption(":blue[Check your ] :red[Gene]",
                          help="https://www.genecards.org/Search/Keyword?queryString=" + selected_gene)
        df_single_gene = df.loc[[selected_gene]]

        separation_opt = ['Median', 'Percentage decomposition', 'Percentage decomposition (lower and upper limit)',
                          'Automatic', 'Expression values cutoff']
        sep_opt = sg_col[1].selectbox(":blue[Please Select the ] :red[Method of Sample Separation]", options=separation_opt,
                                      key='sep_opt_sg', index=3, help=dsc.sep_method_sg())

        # Handling different separation options
        if sep_opt == "Median":
            surv_df = clinic_df[['sample', 'Time', 'Status']]
            surv_df.set_index('sample', inplace=True)
            df_single_gene = df_single_gene.transpose()
            merged_df = df_single_gene.join(surv_df)
            median_value = merged_df[selected_gene].median()
            merged_df['Group'] = np.where(merged_df[selected_gene] <= median_value, 1, 2)
            st.write(f":blue[Median value: ] :red[{median_value}]")
            merged_df = merged_df.replace({'Group': {1: 'Low', 2: 'High'}})
            merged_df.sort_values(by=['Group'], inplace=True)
            surv_df = merged_df[['Time', 'Status', 'Group']]
            surv_df.rename(columns={'Group': 'Cluster'}, inplace=True)
            th_num = median_value

        elif sep_opt == 'Percentage decomposition':
            surv_df = clinic_df[['sample', 'Time', 'Status']]
            surv_df.set_index('sample', inplace=True)
            threshold_per = sg_form.slider(":blue[Set the ] :red[Threshold]", max_value=90, min_value=10, value=50,
                                           key='ssgsea_per_sg')
            df_single_gene = df_single_gene.transpose()
            merged_df = df_single_gene.join(surv_df)
            threshold = merged_df[selected_gene].quantile(threshold_per / 100)
            merged_df['Group'] = np.where(merged_df[selected_gene] < threshold, 1, 2)
            st.write(f":blue[Threshold value: ] :red[{threshold}]")
            if int(merged_df['Group'].nunique()) < 2:
                st.error(':blue[You have only ] :red[ONE ] :blue[group. Please adjust the threshold!]')
            merged_df = merged_df.replace({'Group': {1: 'Low', 2: 'High'}})
            merged_df.sort_values(by=['Group'], inplace=True)
            surv_df = merged_df[['Time', 'Status', 'Group']]
            surv_df.rename(columns={'Group': 'Cluster'}, inplace=True)
            th_num = threshold

        elif sep_opt == 'Percentage decomposition (lower and upper limit)':
            surv_df = clinic_df[['sample', 'Time', 'Status']]
            surv_df.set_index('sample', inplace=True)
            threshold_per = sg_form.slider(":blue[Set the ] :red[Threshold]", 10, 90, (33, 66),
                                           key="perc_ssgsea_up_down_sg")
            df_single_gene = df_single_gene.transpose()
            merged_df = df_single_gene.join(surv_df)
            lower_threshold = merged_df[selected_gene].quantile(threshold_per[0] / 100)
            upper_threshold = merged_df[selected_gene].quantile(threshold_per[1] / 100)
            threshold_up_down_col = st.columns(5)
            threshold_up_down_col[0].write(f":blue[Lower threshold: ] :red[{lower_threshold}]")
            threshold_up_down_col[4].write(f":blue[Upper threshold: ] :red[{upper_threshold}]")
            th_num = threshold_per
            conditions = [
                (merged_df[selected_gene] < lower_threshold),
                (merged_df[selected_gene] >= lower_threshold) & (merged_df[selected_gene] <= upper_threshold),
                (merged_df[selected_gene] > upper_threshold)
            ]
            values = [1, 2, 3]
            merged_df['Group'] = np.select(conditions, values)

            try:
                merged_df = merged_df.replace({'Group': {1: 'Low', 2: 'Intermediate', 3: 'High'}})
                merged_df.sort_values(by=['Group'], inplace=True)
                surv_df = merged_df[['Time', 'Status', 'Group']].copy()
                surv_df.rename(columns={'Group': 'Cluster'}, inplace=True)
                th_num1 = lower_threshold
                th_num2 = upper_threshold

                if surv_df["Cluster"].nunique() < 3:
                    st.warning(
                        ":blue[You have fewer than ] :red[three groups. ] "
                        ":blue[Please adjust the thresholds or choose other genes because the variance of expression is low or 0!]"
                    )
            except:
                st.error(
                    ":blue[You have fewer than ] :red[two groups. ] "
                    ":blue[Please adjust the thresholds or choose other genes because the variance of expression is low or 0!]"
                )

        elif sep_opt == 'Expression values cutoff':
            surv_df = clinic_df[['sample', 'Time', 'Status']]
            surv_df.set_index('sample', inplace=True)
            df_single_gene = df_single_gene.transpose()
            merged_df = df_single_gene.join(surv_df)
            threshold = sg_form.slider(":blue[Set the ] :red[Threshold]",
                                       max_value=df_single_gene[selected_gene].max() - 0.01,
                                       min_value=df_single_gene[selected_gene].min() + 0.01,
                                       value=df_single_gene[selected_gene].median(), key='expression_th_sg')
            merged_df['Group'] = np.where(merged_df[selected_gene] < threshold, 1, 2)
            st.write(f":blue[Threshold value: ] :red[{threshold}]")
            th_num = threshold
            if int(merged_df['Group'].nunique()) < 2:
                st.error(':blue[You have only ] :red[ONE ] :blue[group. Please adjust the threshold!]')
            merged_df = merged_df.replace({'Group': {1: 'Low', 2: 'High'}})
            merged_df.sort_values(by=['Group'], inplace=True)
            surv_df = merged_df[['Time', 'Status', 'Group']]
            surv_df.rename(columns={'Group': 'Cluster'}, inplace=True)

        elif sep_opt == 'Automatic':
            try:
                surv_df = clinic_df[['sample', 'Time', 'Status']]
                surv_df.set_index('sample', inplace=True)
                auto_sg_col = sg_form.columns(2)
                threshold_interval = auto_sg_col[0].slider(":blue[Set the ] :red[Interval]", 5, 95, (10, 90),
                                                           key="perc_automatic_sg", help=dsc.perc_int())

                step_size = int(
                    auto_sg_col[1].select_slider(":blue[Set the ] :red[Step]", options=['1', '5', '10'], value='1',
                                                 key="step_sg", help=dsc.sep_step()))
                df_single_gene = df_single_gene.transpose()
                merged_df = df_single_gene.join(surv_df)
                data_list = []
                min_p_value = float('inf')
                optimal_threshold = None

                for threshold_per in range(threshold_interval[0], threshold_interval[1] + 1, step_size):
                    threshold = merged_df[selected_gene].quantile(threshold_per / 100)
                    merged_df['Group'] = np.where(merged_df[selected_gene] < threshold, "Low", "High")
                    merged_df['High_Group'] = np.where(merged_df['Group'] == "High", 1, 0)
                    group1_time = merged_df[merged_df['Group'] == "Low"]['Time']
                    group2_time = merged_df[merged_df['Group'] == "High"]['Time']
                    event_observed_A = merged_df[merged_df['Group'] == "Low"]['Status']
                    event_observed_B = merged_df[merged_df['Group'] == "High"]['Status']
                    result = logrank_test(group1_time, group2_time, event_observed_A=event_observed_A,
                                          event_observed_B=event_observed_B)
                    p_value = result.p_value

                    cox_df = merged_df[['Time', 'Status', 'High_Group']].dropna(
                        subset=['Time', 'Status', 'High_Group'])

                    if not cox_df.empty:
                        cph = CoxPHFitter()
                        cph.fit(cox_df, 'Time', event_col='Status', fit_options=dict(step_size=0.1))

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
                        'HR_reciprocal': (1 / hazard_ratio_estimate) if hazard_ratio_estimate else None
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
                merged_df['Group'] = np.where(merged_df[selected_gene] < optimal_threshold, "Low", "High")
                merged_df.sort_values(by=['Group'], inplace=True)
                surv_df = merged_df[['Time', 'Status', 'Group']]
                surv_df.rename(columns={'Group': 'Cluster'}, inplace=True)
                th_num = optimal_threshold
                auto_correct = True

            except Exception as e:
                st.warning(f":blue[Convergence halted. ] :red[Attempting to set threshold automatically.] "
                           f":blue[Error Details: {e}]. :red[Gene has low varience, choose another one.]")

                # Initialize list of alternative percentiles to try
                alternative_percentiles = [0.5, 0.25, 0.75]  # median, 25th, 75th
                threshold_set = False
                for perc in alternative_percentiles:
                    try:
                        threshold = merged_df[selected_gene].quantile(perc)
                        merged_df['Group'] = np.where(merged_df[selected_gene] <= threshold, "Low", "High")
                        unique_groups = merged_df['Group'].nunique()

                        if unique_groups >= 2 and merged_df['Group'].isin(['Low', 'High']).all():
                            st.write(f":blue[Threshold set to {int(perc * 100)}th percentile: ] :red[{threshold}]")
                            th_num = threshold
                            threshold_set = True
                            surv_df = merged_df[['Time', 'Status', 'Group']]
                            surv_df.rename(columns={'Group': 'Cluster'}, inplace=True)
                            st.info(
                                f":blue[The threshold has been automatically set to the {int(perc * 100)}th percentile to create two groups.]")
                            break
                        else:
                            st.warning(
                                f":blue[Threshold at {int(perc * 100)}th percentile ({threshold}) did not create two valid groups.]")

                    except Exception as inner_e:
                        st.warning(
                            f":blue[Failed to set threshold at {int(perc * 100)}th percentile due to error: {inner_e}]")

                if not threshold_set:
                    # Additional handling if all percentiles fail
                    st.error(":red[Unable to separate samples into two groups. Please adjust your data or settings.]")
                    st.stop()

                # Additional check: Ensure that 'Low' and 'High' groups have sufficient samples
                group_counts = merged_df['Group'].value_counts()
                if group_counts.min() < 2:
                    st.warning(":blue[One of the groups has fewer than 2 samples. Consider adjusting the threshold.]")

        # Advanced settings for plots
        plot_ext = sg_form.expander(":blue[Advanced Settings of Plots]", expanded=False)
        color_cols = plot_ext.columns(3)

        # Define default colors
        default_colors = {"High": "#e41a1d", "Low": "#377eb8", "Intermediate": "#4daf4a"}

        # Use color picker for selecting colors, falling back to default if not changed by the user
        color1 = color_cols[1].color_picker(":blue[Select Color for ] :red[High Group:]", value=default_colors["High"],
                                            key='color_high')
        color2 = color_cols[0].color_picker(":blue[Select Color for ] :red[Low Group:]", value=default_colors["Low"],
                                            key='color_low')
        color3 = color_cols[2].color_picker(":blue[Select Color for ] :red[Intermediate Group:]",
                                            value=default_colors["Intermediate"],
                                            key='color_intermediate') if 'Intermediate' in surv_df[
            "Cluster"].unique() else None

        # Map selected or default colors to cluster names
        color_mapping = {"High": color1, "Low": color2}
        if color3:
            color_mapping["Intermediate"] = color3

        unique_clusters = surv_df["Cluster"].unique()
        cluster_colors = [color_mapping.get(cluster, "#000000") for cluster in unique_clusters]

        surv_colors = pd.DataFrame({"cluster": unique_clusters, "color": cluster_colors})
        surv_colors = surv_colors[surv_colors['cluster'].isin(unique_clusters)]
        surv_colors.to_csv(surv_colors_path, index=False)

        mv_group_df = pd.DataFrame()
        mv_group_df.index = surv_df.index
        mv_group_df["Cluster"] = surv_df["Cluster"]

        adv_set_multi_sg_col = plot_ext.columns(2)
        min_num = adv_set_multi_sg_col[0].number_input(":blue[Select] :red[minimum number of groups]", min_value=1,
                                                       max_value=10, value=3, key="min_num_sg")

        mv_db1, cols_to_filter, available_columns = sg_ssgsea_multi_data(dataset_option, mv_group_df, surv_df, clinic_df,
                                                                         anno_file)

        # Identify invalid values and set them to NaN
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
            if (mv_db1[col].nunique() >= 2) or (mv_db1[col].nunique() == 2 and not mv_db1[col].isin(['nan']).any())
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

        selected_features = adv_set_multi_sg_col[1].multiselect(":blue[Select the] :red[Features]",
                                                         options=selected_features_option,
                                                         key="multi_adv_g_sg", default=selected_features_option_default,
                                                         max_selections=5)

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
                                              options=unique_values, key=f"select_sg_{feature}", index=initial_index)
            selectbox_df.loc[feature] = selected_val

        if len(selected_features) == mv_db1.shape[1]:
            plot_ext.info(":red[All] :blue[features are selected]")
        else:
            plot_ext.info(f":blue[Selected feature(s): ] :red[{str(selected_features)}]")
        if len(selected_features) != 0:
            selected_features.extend(["Status", "Time"])
            mv_db1 = mv_db1[selected_features]

        mv_db1.replace("NaN", "nan", inplace=True)
        mv_db1.to_csv(multi_db_name)

        if sep_opt == 'Percentage decomposition (lower and upper limit)':
            int_exclusion = plot_ext.toggle(":blue[Exclude ] :red[Intermediate ] :blue[Group]", key='int_exclusion_sg')
            if int_exclusion:
                surv_df = surv_df[surv_df['Cluster'] != 'Intermediate']
                surv_colors = surv_colors[surv_colors['cluster'] != 'Intermediate']

        surv_df.to_csv(surv_dataframe_path)
        surv_colors.to_csv(surv_colors_path, index=False)

        sg_submit_btn = sg_form.form_submit_button("Accept Settings and Start Plotting", help=dsc.act_btn())
        if sg_submit_btn:
            sg_form.info("The changes have been :red[SUCCESSFULLY] refreshed")

        if sep_opt == 'Percentage decomposition (lower and upper limit)':
            cat_type = CategoricalDtype(categories=["High", "Intermediate", "Low"], ordered=True)
        else:
            cat_type = CategoricalDtype(categories=["High", "Low"], ordered=True)

        merged_df['Group'] = merged_df['Group'].astype(cat_type)

        fig = dot_scatter_sg(merged_df, selected_gene, cluster_colors, th_num, th_num1, th_num2, sep_opt)

        if sg_submit_btn:
            try:
                plotly_col = st.columns(3)
                plotly_col[0].subheader(":blue[Expression Profile]", divider='blue')
                plotly_col[0].plotly_chart(fig, use_container_width=True)
                fig.write_image(sg_dot_name, width=600, height=600, scale=1)
            except Exception as e:
                plotly_col = st.columns(3)
                plotly_col[0].image("style_items/error.svg", use_container_width=True)
                plotly_col[0].warning(
                    f":blue[There is problem with the ] :red[chart! ] :blue[Please try another ] :red[settings ] :blue[or ] :red[refresh ] :blue[page] Error: {e}")

            try:
                surv_plot(surv_dataframe_path, surv_png, surv_pdf, user_id, dataset_option, surv_colors_path)
                plotly_col[1].subheader(":blue[Survival Plot]", divider='blue')
                plotly_col[1].image(surv_png, use_container_width=True)
            except Exception as e:
                plotly_col[1].image("style_items/error.svg", use_container_width=True)
                plotly_col[1].warning(
                    f":blue[There is problem with the ] :red[chart! ] :blue[Please try another ] :red[settings ] :blue[or ] :red[refresh ] :blue[page] Error: {e}")

            try:
                multi_var(user_id, selectbox_df, multi_info_name, forest_pdf, forest_png, multi_db_name)
                image = Image.open(forest_png)
                plotly_col[2].subheader(":blue[Multivariate Plot]", divider="blue")
                plotly_col[2].image(image, use_container_width=True, output_format='PNG')
            except Exception as e:
                plotly_col[2].subheader(":blue[Multivariate Plot]", divider="blue")
                plotly_col[2].image("style_items/error.svg", use_container_width=True)
                plotly_col[2].warning(
                    f":blue[You have ] :red[infinity ] :blue[value(s). Try to remove] :red[features] :blue[ or try to use other] :red[ clustering / dimensionality method] Error: {e}")

            if auto_correct:
                if sep_opt == 'Automatic':
                    p_value_col = st.columns(2)
                    p_value_col[0].subheader(":blue[P-Values Comparison Table]", divider="blue")
                    p_value_col[0].dataframe(results_df)
                    lower_limit = results_df.index.min()
                    upper_limit = results_df.index.max()
                    results_df['hover_text_pv'] = results_df.apply(
                        lambda x: f"Threshold: {x['Threshold_Value']}<br>P_Value: {x['P_Value']}", axis=1)
                    results_df['hover_text_pvfdr'] = results_df.apply(
                        lambda x: f"Threshold: {x['Threshold_Value']}<br>P_Value_FDR: {x['P_Value_FDR']}", axis=1)
                    fig_p_value = fdr_p_plot(results_df, lower_limit, upper_limit, pv_dot_name)

                    try:
                        p_value_col[1].subheader(":blue[P-Values Comparison]", divider="blue")
                        p_value_col[1].plotly_chart(fig_p_value, use_container_width=True)
                        st.toast(":red[Single Gene Survival Analysis] :blue[is Completed!]", icon='🎉')
                    except Exception as e:
                        p_value_col[1].image("style_items/error.svg", use_container_width=True)
                        p_value_col[1].warning(
                            f":blue[There is problem with the ] :red[chart! ] :blue[Please try another ] :red[settings ] :blue[or ] :red[refresh ] :blue[page] Error: {e}")

            survival_info(surv_dataframe_path, user_id, surv_type="sg")

        else:
            st.info('If you have entered the plot and analysing settings, please click on the '
                    ':red["Accept Changes and Start Plotting"] button to accept the changes. ')

    else:
        st.warning(":blue[Insufficient data available for] :red[Survival Analysis].")