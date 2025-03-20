import pandas as pd
import numpy as np
import streamlit as st
from scipy.stats import chi2_contingency
from Scripts.r_code_settings import file_copy
from Scripts.plotly_plots import plot_heatmap, plot_component_counts
from Scripts.user_interaction import get_session_id


@st.fragment
def chi_sqrt_test(clust_db, sub_clin_df, cont_name, comp_name, chi_hm_name, chi_bp_name):
    """
    Perform Chi-Squared tests on the given cluster and clinical data.

    Parameters:
    - clust_db (pd.DataFrame): DataFrame containing cluster data.
    - sub_clin_df (pd.DataFrame): DataFrame containing clinical data.
    - cont_name (str): Name for the contingency matrix CSV file.
    - comp_name (str): Name for the component counts CSV file.
    - chi_hm_name (str): Name for the Chi-Squared heatmap image file.
    - chi_bp_name (str): Name for the component counts bar plot image file.

    Returns:
    - None
    """
    nan_values = ['[Discrepancy]', 'nan', '[Not Applicable]', '[Not Evaluated]', '[Not Available]', '[Unknown]']
    user_id = get_session_id()
    try:
        # Joining the cluster data and clinical data
        t_test_df = clust_db[clust_db.index.isin(sub_clin_df.index)]
        t_test_df = t_test_df.join(sub_clin_df)

        # Replacing the specific NaN values
        for nan_value in nan_values:
            t_test_df.replace(nan_value, float('nan'), inplace=True)

        # List of categorical columns to test
        categorical_columns = [
            col
            for col in t_test_df.columns[1:]
            if t_test_df[col].nunique() > 1
        ]

        chi_col = st.columns(2)
        selected_category = chi_col[0].selectbox(":blue[Select a ] :red[Category]", options=categorical_columns, key='cat_sel', index=1)
        selected_category = [selected_category]

        if int(t_test_df[selected_category].nunique().unique()) <= 1:
            st.warning(":blue[Not correct ] :red[data] :blue[for this calculation! Select another ] :red[category]!")
            return None, None
        else:
            corr_boolean = chi_col[1].selectbox(":blue[Use ] :red[P-value Correction]", key='corr_boolean', options=["Apply", "Not Apply"], index=0)
            st.subheader(" ", divider="blue")

            corr_bool = True if corr_boolean == "Apply" else False

            # Split the DataFrame by the 'Cluster' column into a dictionary
            groups = {k: v for k, v in t_test_df.groupby('Cluster')}
            clusters = sorted(groups.keys())

            # Initialize a dictionary to store counts for each cluster
            component_counts = {cluster: {} for cluster in clusters}

            # Initialize a dictionary to hold the matrices
            p_value_matrices = {col: np.zeros((len(clusters), len(clusters))) for col in selected_category}

            # Iterate through all unique combinations of 2 clusters
            for i, group1 in enumerate(clusters):
                for j, group2 in enumerate(clusters[i + 1:], start=i + 1):
                    group1_df = groups[group1]
                    group2_df = groups[group2]

                    # Get the selected category column
                    category_column1 = group1_df[selected_category[0]]
                    category_column2 = group2_df[selected_category[0]]

                    # Count the occurrences of each unique component in the category for each cluster
                    for component, count in category_column1.value_counts().items():
                        component_counts[group1][component] = count
                    for component, count in category_column2.value_counts().items():
                        component_counts[group2][component] = count

                    # Iterate through the categorical columns and perform Chi-Squared tests
                    for col in selected_category:
                        # Create a contingency table
                        contingency_table = pd.crosstab(pd.concat([group1_df[col], group2_df[col]]),
                                                        [len(group1_df) * [group1] + len(group2_df) * [group2]])
                        if contingency_table.size == 0:
                            st.warning(f":red[Warning: ] :blue[Contingency table is empty. Skipping chi-square test between ]"
                                       f":red[{group1_df['Cluster'].unique()[0]}] "
                                       f":blue[and ] :red[{group2_df['Cluster'].unique()[0]}] :blue[clusters. Try another group, or ] "
                                       ":blue[filter the data in the ] :red[Subtype Selector Tab]")
                            # Handle the situation as you see fit
                        else:
                            chi2, p_value, _, _ = chi2_contingency(contingency_table, correction=corr_bool)
                            p_value_matrices[col][j, i] = p_value

            # Convert the component_counts dictionary to a Pandas DataFrame
            component_counts_df = pd.DataFrame.from_dict(component_counts, orient='index')

            # Fill any NaN values with 0 (if some clusters don't contain certain components)
            component_counts_df.fillna(0, inplace=True)
            component_counts_df.index.rename("Cluster", inplace=True)

            chi_table_col = st.columns(2)

            # Display the DataFrame in Streamlit
            chi_table_col[0].subheader(":blue[Component Counts]")
            chi_table_col[0].dataframe(component_counts_df, use_container_width=True)
            component_counts_df.to_csv(f'{comp_name}', index=True)

            # Convert the matrix of the selected category to a DataFrame
            selected_col = selected_category[0]
            p_value_matrix = p_value_matrices[selected_col]
            p_value_df = pd.DataFrame(p_value_matrix, index=clusters, columns=clusters)
            p_value_df.index.rename("Clusters", inplace=True)

            # Display the DataFrame in Streamlit
            chi_table_col[1].subheader(":blue[Contingency Matrix]")
            chi_table_col[1].dataframe(p_value_df, use_container_width=True)
            p_value_df.to_csv(f'{cont_name}', index=True)

            if p_value_df is not None:
                st.subheader(" ", divider="blue")
                chi_plot_col = st.columns(2)
                plot_heatmap(p_value_df, f'{chi_hm_name}', chi_plot_col)
                plot_component_counts(component_counts_df, f'{chi_bp_name}', chi_plot_col)
                st.subheader(" ", divider="blue")
            else:
                pass
        file_copy(user_id)
    except Exception as e:
        st.error(f"An error occurred during the Chi-Squared test: {e}")
