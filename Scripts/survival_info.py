from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test, pairwise_logrank_test
import numpy as np
import pandas as pd
import streamlit as st
from Scripts.user_interaction import naming


def survival_info(surv_df_path, user_id, surv_type):
    # Load the dataset
    data = pd.read_csv(surv_df_path, index_col=None)


    if surv_type == 'cluster':
        survival_info_pw_path = naming(user_id)[65]
        surv_info_path = naming(user_id)[62]
    elif surv_type == 'ssgsea':
        survival_info_pw_path = naming(user_id)[66]
        surv_info_path = naming(user_id)[63]
    else:
        survival_info_pw_path = naming(user_id)[67]
        surv_info_path = naming(user_id)[64]

    # Handle sample IDs: use index or first unnamed column
    if data.index.name or data.index.dtype == 'object':
        data = data.reset_index().rename(columns={"index": "Sample_ID"})
    else:
        data = data.rename(columns={data.columns[0]: "Sample_ID"})

    # Make column names lowercase for case-insensitive matching
    data.columns = [col.lower() for col in data.columns]

    # Define possible keywords for other columns
    time_keywords = [ "time", "months", "duration", "overall"]
    event_keywords = ["os", "event", "status", "censor"]
    cluster_keywords = ["cluster", "group", "category"]

    # Dynamically detect columns
    def get_column(data, keywords, exclude=[]):
        for keyword in keywords:
            matching_columns = [
                col for col in data.columns if keyword.lower() in col and col not in exclude
            ]
            if matching_columns:
                return matching_columns[0]
        raise ValueError(f"None of the keywords {keywords} matched any column in the dataset.")

    try:
        # Detect event column first to avoid conflicts with time column
        event_col = get_column(data, event_keywords)
        time_col = get_column(data, time_keywords, exclude=[event_col])
        cluster_col = get_column(data, cluster_keywords)

    except ValueError as e:
        st.error(f"Error: {e}")
        return

    # Rename columns for consistent processing
    data = data.rename(columns={
        time_col: "survival_time",
        event_col: "event",
        cluster_col: "cluster_id"
    })

    # Drop rows with missing values in required columns
    data = data.dropna(subset=["survival_time", "event", "cluster_id"])

    kmf = KaplanMeierFitter()
    clusters = data['cluster_id'].unique()

    # Store survival information
    median_survival = {}
    cluster_sizes = {}

    for cluster in clusters:
        cluster_data = data[data['cluster_id'] == cluster]
        kmf.fit(cluster_data['survival_time'], cluster_data['event'])
        median_survival[cluster] = kmf.median_survival_time_
        cluster_sizes[cluster] = len(cluster_data)

    # Pairwise log-rank tests and hazard ratios
    pairwise_results = []
    for i, cluster1 in enumerate(clusters):
        for j, cluster2 in enumerate(clusters):
            if i < j:
                cluster1_data = data[data['cluster_id'] == cluster1]
                cluster2_data = data[data['cluster_id'] == cluster2]

                # Ensure valid shapes for log-rank test
                if len(cluster1_data) == 0 or len(cluster2_data) == 0:
                    continue

                test_results = logrank_test(
                    cluster1_data['survival_time'],
                    cluster2_data['survival_time'],
                    cluster1_data['event'],
                    cluster2_data['event']
                )

                # Calculate HR using Cox Proportional Hazard model
                cox_data = data[data['cluster_id'].isin([cluster1, cluster2])].copy()
                cox_data['cluster'] = np.where(cox_data['cluster_id'] == cluster1, 1, 0)
                cph = CoxPHFitter()
                cph.fit(cox_data[['survival_time', 'event', 'cluster']], duration_col='survival_time',
                        event_col='event')

                # Calculate HR and confidence intervals explicitly
                hr = cph.hazard_ratios_['cluster']
                log_hr = np.log(hr)
                ci_lower = np.exp(log_hr - 1.96 * cph.standard_errors_['cluster'])
                ci_upper = np.exp(log_hr + 1.96 * cph.standard_errors_['cluster'])

                pairwise_results.append({
                    "Comparison": f"{cluster1} vs {cluster2}",
                    "p-value": test_results.p_value,
                    "HR": hr,
                    "HR_CI_Lower": ci_lower,
                    "HR_CI_Upper": ci_upper
                })

    # Create a summary dataframe
    summary_df = pd.DataFrame({
        "Cluster_ID": clusters,
        "Cluster_Size": [cluster_sizes[c] for c in clusters],
        "Median_Survival_Time": [median_survival[c] for c in clusters]
    })

    pairwise_df = pd.DataFrame(pairwise_results)

    # Display results
    surv_info_expander = st.expander(":red[Survival ] :blue[Info]", expanded=True)
    surv_info_cols = surv_info_expander.columns(2)

    surv_info_cols[0].subheader(":blue[Cluster Survival Summary]")
    surv_info_cols[0].dataframe(summary_df, use_container_width=True, height=200, hide_index=True)
    surv_info_cols[1].subheader(":blue[Pairwise Comparisons]")
    surv_info_cols[1].dataframe(pairwise_df, use_container_width=True, height=200, hide_index=True)

    summary_df.to_csv(surv_info_path)
    pairwise_df.to_csv(survival_info_pw_path)
