import streamlit as st
import pandas as pd
from scipy.stats import zscore
from Scripts.user_interaction import naming


def zscore_calc(df):
    """
    Calculate z-scores for the dataframe.

    Parameters:
    df (DataFrame): Input dataframe to calculate z-scores.

    Returns:
    df_z (DataFrame): Z-scored dataframe with positive adjustment.
    df_z_neg (DataFrame): Z-scored dataframe without adjustment.
    """
    try:
        # Compute z-scores once and wrap the result in DataFrames to preserve index/columns
        z = zscore(df, axis=1, nan_policy='omit')
        df_z = pd.DataFrame(z, index=df.index, columns=df.columns)
        df_z_neg = df_z.copy()

        # Find the minimum value in the z-scored dataframe
        df_min = df_z.min(numeric_only=True).min()

        # Adjust z-scores to make all values positive if the minimum value is negative
        if df_min < 0:
            df_z += abs(df_min)

        # Drop rows with missing values and record the dropped indices
        dropped_indices = df_z.index[df_z.isnull().any(axis=1)].tolist()
        df_z.dropna(axis=0, inplace=True)
        df_z_neg.dropna(axis=0, inplace=True)

        if dropped_indices:
            st.sidebar.write(
                f":blue[After Z-Scoring we ] :red[dropped ] :blue[following gene(s) (no standard deviation): ] "
                f":red[{dropped_indices}]"
            )

        return df_z_neg, df_z
    except Exception as e:
        st.error(f"Error calculating z-scores: {e}")
        return None, None


def _process_zscore_data(df, csv_name, csv_name_nmf, user_id):
    """
    Internal helper function to process data with z-score toggles and perform final transformations.

    Parameters:
    df (DataFrame): Input gene expression data.
    csv_name (str): File path to save selected genes.
    csv_name_nmf (str): File path to save selected genes for NMF.
    user_id (str): User identifier.

    Returns:
    tuple: (processed_df, processed_nmf, df_names, random_state_samples, df_names_gene)
    """
    # Create columns for layout
    z_col = st.columns(6)
    z_col[3].write(" ")
    z_col[3].write(":blue[Use ] :red[**Z-Score:**]")

    # Checkboxes for using z-score in clustering and heatmap
    z_score_chb_clust = z_col[4].toggle(":blue[Clustering]", value=True, key="z_score_toggle_clustering")
    z_score_chb_hm = z_col[5].toggle(":blue[Heatmap]", value=True, key="z_score_toggle_heatmap_heatmap")

    # Process based on user selection
    if z_score_chb_clust or z_score_chb_hm:
        df_zscore, df_zscore_neg = zscore_calc(df)
        if z_score_chb_clust and z_score_chb_hm:
            out_df = df_zscore
            out_nmf = df_zscore_neg
            out_df.to_csv(csv_name)
            out_nmf.to_csv(csv_name_nmf)
        elif z_score_chb_clust and not z_score_chb_hm:
            out_df = df_zscore
            out_nmf = df_zscore_neg
            # Save original data instead of z-scored version for CSV export
            df.to_csv(csv_name)
            df.to_csv(csv_name_nmf)
            out_df = out_df.loc[out_df.index.isin(df.index)]
        elif not z_score_chb_clust and z_score_chb_hm:
            out_df = df_zscore
            out_nmf = df_zscore_neg
            out_df.to_csv(csv_name)
            out_nmf.to_csv(csv_name_nmf)
            df = df.loc[df.index.isin(out_df.index)]
            out_df = df
        else:
            out_df = df
            out_nmf = df
    else:
        df.to_csv(csv_name)
        df.to_csv(csv_name_nmf)
        out_df = df
        out_nmf = df

    # Finalize: transpose the data once and extract gene names
    proc_df = out_df.transpose()
    proc_nmf = out_nmf.transpose()
    df_names = proc_df.index
    # df_names_gene is taken as the index of the original orientation
    df_names_gene = proc_df.transpose().index

    # Number input for random state (for reproducibility)
    random_state_samples = z_col[0].number_input(
        ":blue[Set the ] :red[Random State]",
        min_value=1,
        max_value=2**32 - 1,
        value=42,
        help=("In machine learning algorithms, the random state is a parameter that controls the randomness "
              "to ensure reproducibility. By setting this value, users can achieve consistent results across "
              "multiple runs of the same algorithm, aiding in debugging and model comparison.")
    )
    proc_df = proc_df.dropna(axis=1)

    # Save the actual gene list to a CSV file
    final_gene_list_path = naming(user_id)[36]
    actual_gene_list = pd.DataFrame(proc_df.transpose().index)
    actual_gene_list.to_csv(final_gene_list_path, index=False, header=None)

    return proc_df, proc_nmf, df_names, random_state_samples, df_names_gene


def data_preprocess(selected_genes_name, selected_genes_name_nmf, selected_genes_filtered, name_list, user_id):
    """
    Preprocess the data by applying z-score normalization and saving the results for Neuroendocrine Lung Cohorts.

    Parameters:
    selected_genes_name (str): File path to save selected genes.
    selected_genes_name_nmf (str): File path to save selected genes for NMF.
    selected_genes_filtered (DataFrame): Filtered selected genes.
    name_list (list): List of sample names.
    user_id (str): User identifier.

    Returns:
    tuple: (selected_genes1, selected_genes_nmf, df_names, random_state_samples, df_names_gene)
    """
    try:
        if len(name_list) >= 8:
            return _process_zscore_data(selected_genes_filtered, selected_genes_name, selected_genes_name_nmf, user_id)
        else:
            st.warning(":blue[The number of samples should be at ] :red[least 8.]")
            return None, None, None, None, None
    except Exception as e:
        st.error(f"Error in data preprocessing: {e}")
        return None, None, None, None, None


def post_process(target_cancer, selected_genes_name, selected_genes_name_nmf, user_id):
    """
    Post-process the target cancer (TCGA) data by applying z-score normalization and saving the results for TCGA cohorts.

    Parameters:
    target_cancer (DataFrame): Target cancer data.
    selected_genes_name (str): File path to save selected genes.
    selected_genes_name_nmf (str): File path to save selected genes for NMF.
    user_id (str): User identifier.

    Returns:
    tuple: (selected_genes1, selected_genes_nmf, df_names, random_state_samples, df_names_gene)
    """
    try:
        name_list = target_cancer.columns
        if len(name_list) >= 8:
            return _process_zscore_data(target_cancer, selected_genes_name, selected_genes_name_nmf, user_id)
        else:
            st.warning(":blue[The number of samples should be at ] :red[least 8.]")
            return None, None, None, None, None
    except Exception as e:
        st.error(f"Error in post-processing: {e}")
        return None, None, None, None, None
