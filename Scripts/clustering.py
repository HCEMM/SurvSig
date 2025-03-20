import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from Scripts.clust_dim_def import (hierarchical_clustering, gmm, kmeans, som, dtc, hdb, optics,
                                   hierarchical_clustering_gene, gmm_gene, kmeans_gene, som_gene, dtc_gene, hdb_gene,
                                   optics_gene)
from Scripts.plotly_plots import clust_scatter, clust_scatter3d
from Scripts.r_code_settings import ht_adv_set
from Scripts.rpy2_heatmap_plots import corr_heatmap
from Scripts.r_code_settings import rgb_to_hex
import Scripts.descriptions as dscript
from Scripts.user_interaction import naming


def save_disp(df_names, cluster_name, lat_vec, selected_genes_dim_red,
              labels, option_dataset, clinic_df, surv_dataframe_path, user_session_id,
              ht_expression_path, column_cluster_path, ht_top_annotation_path,
              ht_png, ht_pdf, surv_png, surv_pdf, ht_row_cluster_path,
              selected_subtypes, scatter_2d_name, scatter_3d_name, adv_col2, selected_genes_dm_path,
              opt_dim_red, cluster_toggle, e, cancer_opt=None):
    """
    Save and display clustering results.

    Parameters:
    - df_names (list): List of sample names.
    - cluster_name (str): Filename to save clustering results.
    - lat_vec (int): Number of latent vectors.
    - selected_genes_dim_red (DataFrame): Dimensionally reduced data.
    - labels (array): Cluster labels.
    - option_dataset (str): Selected dataset option.
    - clinic_df (DataFrame): Clinical data.
    - surv_dataframe_path (str): Path to save survival analysis data.
    - user_session_id (str): User identifier.
    - ht_expression_path (str): Path to save heatmap expression data.
    - column_cluster_path (str): Path to save column cluster data.
    - ht_top_annotation_path (str): Path to save top annotation data.
    - ht_png (str): Path to save heatmap PNG.
    - ht_pdf (str): Path to save heatmap PDF.
    - surv_png (str): Path to save survival analysis PNG.
    - surv_pdf (str): Path to save survival analysis PDF.
    - ht_row_cluster_path (str): Path to save heatmap row cluster data.
    - selected_subtypes (str): Selected subtypes.
    - scatter_2d_name (str): Filename for 2D scatter plot.
    - scatter_3d_name (str): Filename for 3D scatter plot.
    - adv_col2 (streamlit.Column): Streamlit column for advanced settings.
    - selected_genes_dm_path (str): Path to save dimensionally reduced data.
    - opt_dim_red (str): Selected dimensionality reduction method.
    - cancer_opt (str, optional): Selected cancer option.

    Returns:
    - display (DataFrame): DataFrame with clustering results.
    """
    try:
        sample_type = "sample"
        display = pd.DataFrame()
        display["Samples"] = df_names
        display["Cluster"] = pd.DataFrame(labels)
        display.sort_values(by=['Cluster'], inplace=True)
        display.to_csv(f'{cluster_name}', index=False)
        exp_data = pd.read_csv(f'{ht_expression_path}')

        if len(selected_genes_dim_red) != 0 and selected_genes_dim_red.shape[1] <= 3:
            if int(exp_data.shape[1]) != 1:
                if lat_vec == 3:
                    clust_scatter3d(selected_genes_dim_red, labels, df_names, sample_type=sample_type,
                                    cluster_toggle=cluster_toggle, e=e)
                    dl_info_col = st.columns(3)
                    dl_info_col[1].info(
                        ":blue[You can ] :red[download] :blue[the 3D plot in ] :red[SVG ] :blue[, if you click on the ]"
                        ':red["Camera "] :blue[icon]')
                elif lat_vec == 2:
                    clust_scatter(selected_genes_dim_red, labels, df_names, scatter_2d_name, sample_type=sample_type,
                                  cluster_toggle=cluster_toggle, e=e)
        else:
            if len(display["Cluster"].unique()) <= 9:
                palette = px.colors.qualitative.Set1
            else:
                palette = px.colors.qualitative.Light24
            color_file = pd.DataFrame()
            palette = palette[:len(display["Cluster"].unique())]
            color_file["type"] = display["Cluster"].unique()
            color_file["column"] = "Cluster"
            color_file.sort_values(ascending=True, by="type", inplace=True)
            color_file["color"] = palette
            color_file['color'] = color_file['color'].apply(rgb_to_hex)
            user_dir = f"/result_data_{user_session_id}"
            anno_color_corr = naming(user_session_id)[29]
            color_file.to_csv(anno_color_corr, index=0)
            title = "Samples"

            if opt_dim_red == "Correlation":
                corr_heatmap(user_session_id, selected_genes_dm_path, column_cluster_path, title, anno_color_corr, e,
                             cluster_toggle)

        if option_dataset in ["Anish-SCLC", "Jiang-SCLC", "George-LCNEC", "Fernandez-Carcinoid",
                              "Alcala-Carcinoid", "Liu-SCLC", "Rousseaux-Mixed(non-NE)", "George-SCLC"]:
            display.reset_index()
            display.set_index("Samples", inplace=True)
            clinic_df_to_r = clinic_df[clinic_df.index.isin(display.index)]
            clinic_df_to_r["Cluster"] = display["Cluster"]
            clinic_df_to_r.to_csv(f'{surv_dataframe_path}')

            ht_adv_set(column_cluster_path, ht_row_cluster_path, user_session_id, ht_top_annotation_path,
                       option_dataset,
                       f'{surv_dataframe_path}', ht_png, ht_pdf, surv_png, surv_pdf, ht_expression_path,
                       cancer_opt)
        elif option_dataset == "TCGA":
            display.reset_index()
            display.set_index("Samples", inplace=True)
            clinic_df_to_r = clinic_df[clinic_df.index.isin(display.index)]
            clinic_df_to_r["Cluster"] = display["Cluster"]
            clinic_df_to_r.to_csv(f'{surv_dataframe_path}')

            ht_adv_set(column_cluster_path, ht_row_cluster_path, user_session_id, ht_top_annotation_path,
                       option_dataset,
                       f'{surv_dataframe_path}', ht_png, ht_pdf, surv_png, surv_pdf, ht_expression_path,
                       cancer_opt)

        return display
    except Exception as e:
        st.error(f"Error in save_disp: {e}")
        return None


def save_disp_genes(df_names, cluster_name, lat_vec, selected_genes_dim_red, labels, ht_expression_path,
                    gene_exp, scatter_2d_name, scatter_3d_name, adv_col2_gene, user_session_id,
                    column_cluster_path, selected_genes_dm_path, cluster_toggle, e):
    """
    Save and display clustering results for genes.

    Parameters:
    - df_names (list): List of gene names.
    - cluster_name (str): Filename to save clustering results.
    - lat_vec (int): Number of latent vectors.
    - selected_genes_dim_red (DataFrame): Dimensionally reduced data.
    - labels (array): Cluster labels.
    - ht_expression_path (str): Path to save heatmap expression data.
    - gene_exp (streamlit.Column): Streamlit column for gene expression settings.
    - scatter_2d_name (str): Filename for 2D scatter plot.
    - scatter_3d_name (str): Filename for 3D scatter plot.
    - adv_col2_gene (streamlit.Column): Streamlit column for advanced settings.
    - user_session_id (str): User identifier.
    - column_cluster_path (str): Path to save column cluster data.
    - selected_genes_dm_path (str): Path to save dimensionally reduced data.

    Returns:
    - display (DataFrame): DataFrame with clustering results.
    """
    try:
        sample_type = "gene"
        display = pd.DataFrame()
        display["Samples"] = pd.DataFrame(df_names)
        display["Cluster"] = pd.DataFrame(labels)
        display.sort_values(by=['Cluster'], inplace=True)
        display.to_csv(f'{cluster_name}', index=False)
        exp_data = pd.read_csv(f'{ht_expression_path}')

        if exp_data.shape[0] > 3 >= selected_genes_dim_red.shape[1]:
            gene_scatter_chb = gene_exp.toggle(":blue[Show ] :red[Plot]", value=True, key="plot_chb")
            if gene_scatter_chb:
                if len(selected_genes_dim_red) != 0:
                    if int(exp_data.shape[0]) >= 1:
                        if lat_vec == 3:
                            clust_scatter3d(selected_genes_dim_red, labels, df_names, sample_type=sample_type,
                                            cluster_toggle=cluster_toggle, e=e)
                        elif lat_vec == 2:
                            clust_scatter(selected_genes_dim_red, labels, df_names, scatter_2d_name, sample_type,
                                          cluster_toggle=cluster_toggle, e=e)
        else:
            gene_scatter_chb = gene_exp.toggle(":blue[Show ] :red[Plot]", value=True, key="plot_chb2")
            if gene_scatter_chb:
                if len(display["Cluster"].unique()) <= 9:
                    palette = px.colors.qualitative.Set1
                else:
                    palette = px.colors.qualitative.Light24
                color_file = pd.DataFrame()
                palette = palette[:len(display["Cluster"].unique())]
                color_file["type"] = display["Cluster"].unique()
                color_file["column"] = "Cluster"
                color_file.sort_values(ascending=True, by="type", inplace=True)
                color_file["color"] = palette
                color_file['color'] = color_file['color'].apply(rgb_to_hex)
                user_dir = f"/result_data_{user_session_id}"
                anno_color_corr = naming(user_session_id)[35]
                color_file.to_csv(anno_color_corr, index=0)
                title = "Genes"
                corr_heatmap(user_session_id, selected_genes_dm_path, column_cluster_path, title, anno_color_corr, e,
                             cluster_toggle)

        return display
    except Exception as e:
        st.error(f"Error in save_disp_genes: {e}")
        return None


def clustering(selected_genes_dim_red, option_dm, filtered_selected_genes, option_dataset, df_names, cluster_name,
               lat_vec, user_session_id, random_state_samples, adv_col1, adv_col2, ht_expression_path,
               column_cluster_path, ht_top_annotation_path, ht_png, ht_pdf, surv_png, surv_pdf, ht_row_cluster_path,
               selected_subtypes, cancer_opt, scatter_2d_name, scatter_3d_name, selected_genes_dm_path,  cluster_toggle,
               sidebar_expander,
               e, clinic_df=None,
               surv_dataframe_path=None):
    """
    Perform clustering on selected genes.

    Parameters:
    - selected_genes_dim_red (DataFrame): Dimensionally reduced data.
    - option_dm (str): Selected dimensionality reduction method.
    - filtered_selected_genes (DataFrame): Filtered gene expression data.
    - option_dataset (str): Selected dataset option.
    - df_names (list): List of sample names.
    - cluster_name (str): Filename to save clustering results.
    - lat_vec (int): Number of latent vectors.
    - user_session_id (str): User identifier.
    - random_state_samples (int): Random state for reproducibility.
    - adv_col1 (streamlit.Column): Streamlit column for advanced settings.
    - adv_col2 (streamlit.Column): Streamlit column for advanced settings.
    - ht_expression_path (str): Path to save heatmap expression data.
    - column_cluster_path (str): Path to save column cluster data.
    - ht_top_annotation_path (str): Path to save top annotation data.
    - ht_png (str): Path to save heatmap PNG.
    - ht_pdf (str): Path to save heatmap PDF.
    - surv_png (str): Path to save survival analysis PNG.
    - surv_pdf (str): Path to save survival analysis PDF.
    - ht_row_cluster_path (str): Path to save heatmap row cluster data.
    - selected_subtypes (str): Selected subtypes.
    - cancer_opt (str): Selected cancer option.
    - scatter_2d_name (str): Filename for 2D scatter plot.
    - scatter_3d_name (str): Filename for 3D scatter plot.
    - selected_genes_dm_path (str): Path to save dimensionally reduced data.
    - clinic_df (DataFrame, optional): Clinical data.
    - surv_dataframe_path (str, optional): Path to save survival analysis data.

    Returns:
    - cluster_db (DataFrame): Clustering database.
    """
    try:
        cluster_db = None
        n_clust = None
        cluster_db = pd.DataFrame()
        labels = None

        cluster_options, clust_index = get_cluster_options(filtered_selected_genes, option_dm)

        num_of_clust = get_number_of_clusters(filtered_selected_genes)

        option_cluster = sidebar_expander.selectbox(':blue[Clustering ] :red[Method]',
                                              options=cluster_options,
                                              help="Clustering methods are unsupervised machine learning techniques used to"
                                                   "group similar data points together based on their attributes or "
                                                   "characteristics. There are several types of clustering algorithms, "
                                                   "including partition-based methods (such as k-means), hierarchical "
                                                   "methods (such as agglomerative clustering, dynamic cut tree), "
                                                   "model-based methods (such as Gaussian mixture models) and"
                                                   "artificial intelligence "
                                                   "model-based methods (such as Self-Organizing Map)",
                                              key="pot_clust", index=clust_index)

        if option_cluster in ["K-means", "Gaussian Mixture Model", "Agglomerative Clustering", 'Self-Organizing Map']:
            n_clust = st.sidebar.select_slider(
                ':blue[Select the ] :red[Number of Clusters]', options=num_of_clust, key="lust_num")
            st.sidebar.write(f':blue[Number of cluster] :red[{str(n_clust)}]')
        st.sidebar.header(" ", divider="blue")

        if option_cluster == "K-means":
            labels = perform_kmeans_clustering(selected_genes_dim_red, n_clust, random_state_samples, adv_col1)
        elif option_cluster == 'Gaussian Mixture Model':
            labels = perform_gmm_clustering(selected_genes_dim_red, n_clust, random_state_samples, adv_col1)
        elif option_cluster == 'Agglomerative Clustering':
            labels = perform_hierarchical_clustering(selected_genes_dim_red, n_clust, adv_col1)
        elif option_cluster == 'Self-Organizing Map':
            labels = perform_som_clustering(selected_genes_dim_red, n_clust, random_state_samples, adv_col1)
        elif option_cluster == 'DynamicTreeCut':
            labels = perform_dtc_clustering(selected_genes_dim_red, filtered_selected_genes, adv_col1)
        elif option_cluster == 'HDBSCAN':
            labels = perform_hdb_clustering(selected_genes_dim_red, filtered_selected_genes, lat_vec, adv_col1)
        elif option_cluster == "OPTICS":
            labels = perform_optics_clustering(selected_genes_dim_red, filtered_selected_genes, adv_col1)
        else:
            st.error(
                ":blue[Select a clustering method and the number of clusters ] "
                ":red[(Except DynamicCutTree, HDBSCAN and OPTICS)! ] "
                ":blue[You can apply dimensionality reduction method if you want!]")
            return cluster_db

        labels = labels + 1
        labels = np.array(["S" + str(number) for number in labels])

        cluster_db = save_disp(df_names, cluster_name, lat_vec, selected_genes_dim_red,
                               labels, option_dataset, clinic_df, surv_dataframe_path, user_session_id,
                               ht_expression_path, column_cluster_path, ht_top_annotation_path,
                               ht_png, ht_pdf, surv_png, surv_pdf, ht_row_cluster_path,
                               selected_subtypes, scatter_2d_name, scatter_3d_name, adv_col2, selected_genes_dm_path,
                               option_dm, cluster_toggle, e, cancer_opt)

        return cluster_db

    except Exception as e:
        st.error(f"Error during clustering: {e}")
        return None


def perform_kmeans_clustering(selected_genes_dim_red, n_clust, random_state_samples, adv_col1):
    """
    Perform K-means clustering.

    Parameters:
    - selected_genes_dim_red (DataFrame): Dimensionally reduced data.
    - n_clust (int): Number of clusters.
    - random_state_samples (int): Random state for reproducibility.
    - adv_col1 (streamlit.Column): Streamlit column for advanced settings.

    Returns:
    - labels (array): Cluster labels.
    """
    try:
        adv_col1.subheader(":blue[Advanced Settings of K-means]", help=dscript.kmeans_adv())
        km_init = adv_col1.selectbox(":blue[K-means ] :red[Initialisation ] :blue[Method]",
                                     options=("random", "k-means++"), index=0, key="km_init")
        km_tol = adv_col1.number_input(":blue[Relative ] :red[Tolerance]", min_value=0.0001, max_value=1.0000,
                                       step=0.0001, value=0.0001, key="km_tol", format="%f")
        mi_km = adv_col1.number_input(":red[Max ] :blue[Iteration]", min_value=200, max_value=2000, step=1, value=1000,
                                      key="km_max_iter")
        km_alg = adv_col1.selectbox(":red[K-means Algorithm] :blue[to Use]", options=("lloyd", "elkan"))

        return kmeans(selected_genes_dim_red, n_clust, km_init, km_tol, mi_km, km_alg, random_state_samples)
    except Exception as e:
        st.error(f"Error in perform_kmeans_clustering: {e}")
        return None


def perform_gmm_clustering(selected_genes_dim_red, n_clust, random_state_samples, adv_col1):
    """
    Perform Gaussian Mixture Model clustering.

    Parameters:
    - selected_genes_dim_red (DataFrame): Dimensionally reduced data.
    - n_clust (int): Number of clusters.
    - random_state_samples (int): Random state for reproducibility.
    - adv_col1 (streamlit.Column): Streamlit column for advanced settings.

    Returns:
    - labels (array): Cluster labels.
    """
    try:
        adv_col1.subheader(":blue[Advanced Settings of Gaussian Mixture Model]", help=dscript.gmm_adv())
        init_method = adv_col1.selectbox(":red[Initialisation ] :blue[Methods]",
                                         options=("kmeans", "k-means++", "random", "random_from_data"),
                                         index=1, key="in_met_gmm")
        covariance_type = adv_col1.selectbox(":red[Convariance] :blue[Type]", options=("full", "tied"), index=0,
                                             key="cov_typ")
        tol_gmm = adv_col1.select_slider(":blue[Convergence ] :red[Threshold]", options=(
            0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1),
                                         value=0.0001, key="tol_gmm_ss")
        reg_cov_gmm = adv_col1.select_slider(':blue[Non-negative ] :red[Regularization]',
                                             options=(1e-3, 1e-4, 1e-5, 1e-6, 1e-7), value=1e-5,
                                             key='reg_covar_gmm')

        return gmm(selected_genes_dim_red, n_clust, init_method, covariance_type, tol_gmm, random_state_samples,
                   reg_cov_gmm)
    except Exception as e:
        st.error(f"Error in perform_gmm_clustering: {e}")
        return None


def perform_hierarchical_clustering(selected_genes_dim_red, n_clust, adv_col1):
    """
    Perform hierarchical clustering.

    Parameters:
    - selected_genes_dim_red (DataFrame): Dimensionally reduced data.
    - n_clust (int): Number of clusters.
    - adv_col1 (streamlit.Column): Streamlit column for advanced settings.

    Returns:
    - labels (array): Cluster labels.
    """
    try:
        adv_col1.subheader(":blue[Advanced Settings of Agglomerative Clustering]", help=dscript.agglo_adv())
        link = adv_col1.selectbox(":red[Linkage Method ] :blue[of Agglomerative Clustering]",
                                  options=("ward", "complete", "average", 'single'), index=0, key="link_hc")
        if link == "ward":
            hc_met = "euclidean"
        else:
            hc_met = adv_col1.selectbox(":red[Metric ] :blue[of Agglomerative Clustering]",
                                        options=('euclidean', 'l1', 'l2', 'manhattan', 'cosine'),
                                        index=0, key="metric_hc")

        return hierarchical_clustering(selected_genes_dim_red, n_clust, link, hc_met)
    except Exception as e:
        st.error(f"Error in perform_hierarchical_clustering: {e}")
        return None


def perform_som_clustering(selected_genes_dim_red, n_clust, random_state_samples, adv_col1):
    """
    Perform Self-Organizing Map clustering.

    Parameters:
    - selected_genes_dim_red (DataFrame): Dimensionally reduced data.
    - n_clust (int): Number of clusters.
    - random_state_samples (int): Random state for reproducibility.
    - adv_col1 (streamlit.Column): Streamlit column for advanced settings.

    Returns:
    - labels (array): Cluster labels.
    """
    try:
        adv_col1.subheader(":blue[Advanced Settings of Self-Organizing Map]", help=dscript.som_adv())
        neg_som = adv_col1.selectbox(":red[Neighbours ] :blue[Function Type]",
                                     options=('gaussian', 'mexican_hat', 'bubble', 'triangle'), index=0,
                                     key="neg_som")
        sig_som = adv_col1.slider(":red[Sigma ] :blue[Value]", key="sig_som", min_value=0.10,
                                  max_value=50.00, value=1.50, step=0.10)
        lr_som = adv_col1.slider(":red[Learning] :blue[ Rate]", key="lr_som", min_value=0.10,
                                 max_value=5.00, value=0.50)
        epoch = adv_col1.number_input(":blue[Number of ] :red[Epochs]", min_value=1, max_value=2000, value=500,
                                      step=1, key="som_epo")
        act_dist = adv_col1.selectbox(":blue[Choose the ] :red[Activation Distance]", key='act_dist_s',
                                      options=("euclidean", "cosine", "manhattan", "chebyshev"), index=0)

        topology = adv_col1.selectbox(":blue[Select the ] :red[Topology of the Map]", key="tplgy_s",
                                      options=("rectangular", "hexagonal"), index=0)

        return som(selected_genes_dim_red, n_clust, neg_som, sig_som, lr_som, epoch, random_state_samples, act_dist,
                   topology)
    except Exception as e:
        st.error(f"Error in perform_som_clustering: {e}")
        return None


def perform_dtc_clustering(selected_genes_dim_red, filtered_selected_genes, adv_col1):
    """
    Perform Dynamic Tree Cut clustering.

    Parameters:
    - selected_genes_dim_red (DataFrame): Dimensionally reduced data.
    - filtered_selected_genes (DataFrame): Filtered gene expression data.
    - adv_col1 (streamlit.Column): Streamlit column for advanced settings.

    Returns:
    - labels (array): Cluster labels.
    """
    try:
        adv_col1.subheader(":blue[Advanced Settings of DynamicTreeCut]", help=dscript.dtc_adv())
        link = adv_col1.selectbox(":red[Linkage ] :blue[Method]", options=(
            "single", "complete", "average", "weighted", "centroid", "median", "ward"), index=6, key="link_dct")
        distance = adv_col1.selectbox(":red[Distance ] :blue[Metric]",
                                      options=("braycurtis", "canberra", "chebyshev", "cityblock",
                                               "correlation", "cosine", "euclidean",
                                               "matching", "minkowski", "seuclidean",
                                               "sqeuclidean"), key="dct_dist_1" if selected_genes_dim_red.shape[1] >
                                      selected_genes_dim_red.shape[0]
                                      else "dct_dist_2", index=6)

        min_clust_size_dtc = adv_col1.number_input(":red[Minimum ] :blue[Cluster Size]", min_value=2, step=1, value=2,
                                                   max_value=int(filtered_selected_genes.shape[0] / 2), key="min_clust")

        labels = dtc(selected_genes_dim_red, link, distance, min_clust_size_dtc)
        return labels["labels"] - 1
    except Exception as e:
        st.error(f"Error in perform_dtc_clustering: {e}")
        return None


def perform_hdb_clustering(selected_genes_dim_red, filtered_selected_genes, lat_vec, adv_col1):
    """
    Perform HDBSCAN clustering.

    Parameters:
    - selected_genes_dim_red (DataFrame): Dimensionally reduced data.
    - filtered_selected_genes (DataFrame): Filtered gene expression data.
    - lat_vec (int): Number of latent vectors.
    - adv_col1 (streamlit.Column): Streamlit column for advanced settings.

    Returns:
    - labels (array): Cluster labels.
    """
    try:
        adv_col1.subheader(":blue[Advanced Settings of HDBSCAN]", help=dscript.hdb_adv())
        alfa = adv_col1.number_input(":red[Alpha ] :blue[Value]", min_value=0.1, max_value=5.0, value=1.0, step=0.1,
                                     key="alfa_dtc", format='%f')
        hdb_met = adv_col1.selectbox(":red[Metric ] :blue[Method]", options=(
            "braycurtis", "canberra", "chebyshev", "cityblock", "euclidean", "hamming",
            "infinity", "manhattan") if lat_vec == 3 else (
            "braycurtis", "canberra", "chebyshev", "cityblock", "euclidean", "haversine", "hamming",
            "infinity", "l1", "l2", "manhattan", "matching"), key="hdb_metric_2d" if lat_vec == 3 else "hdb_metric_3d",
                                     index=4)
        mcs_hdb = adv_col1.number_input(":red[Minimum Number of Samples ] :blue[in Each Cluster]", min_value=2, step=1,
                                        max_value=int(filtered_selected_genes.shape[0] / 2), value=2, key="mcs_key")
        min_sam_hdb = adv_col1.number_input(":red[Minimum Number of Samples ] :blue[in Each Cluster]", min_value=2,
                                            step=1, max_value=int(filtered_selected_genes.shape[0] / 2), key="min_sam_")
        eps = adv_col1.number_input(":blue[Cluster Selection ] :red[Epsilon Value]", min_value=0.0,
                                    max_value=float(filtered_selected_genes.shape[0] / 2), value=1.0, step=0.1,
                                    key="eps_slider", format="%f")
        clust_meth = adv_col1.selectbox(":blue[Select the ] :red[Cluster Selection Method]", key='clust_sel',
                                        options=("eom", "leaf"), index=0)

        return hdb(selected_genes_dim_red, alfa, hdb_met, mcs_hdb, min_sam_hdb, eps, clust_meth)
    except Exception as e:
        st.error(f"Error in perform_hdb_clustering: {e}")
        return None


def perform_optics_clustering(selected_genes_dim_red, filtered_selected_genes, adv_col1):
    """
    Perform OPTICS clustering.

    Parameters:
    - selected_genes_dim_red (DataFrame): Dimensionally reduced data.
    - filtered_selected_genes (DataFrame): Filtered gene expression data.
    - adv_col1 (streamlit.Column): Streamlit column for advanced settings.

    Returns:
    - labels (array): Cluster labels.
    """
    try:
        adv_col1.subheader(":blue[Advanced Settings of OPTICS]", help=dscript.opt_adv())
        ms_optics = adv_col1.number_input(":red[Minimum ] :blue[Samples]", min_value=2, step=1,
                                          max_value=int(filtered_selected_genes.shape[0] / 2), value=5, key="ms_optics")
        metric_opt = adv_col1.selectbox(":red[Distance ] :blue[Metric]", options=(
            "minkowski", "cityblock", "cosine", "euclidean", "l1", "l2", "manhattan", "braycurtis",
            "canberra", "chebyshev", "correlation", "sqeuclidean"), key="metric_opt", index=3)
        cluster_met = adv_col1.selectbox(":red[Initialisation ] :blue[Method]", options=("xi", "dbscan"),
                                         key="clust_met", index=1)
        eps_opt = adv_col1.number_input(":blue[Cluster Selection ] :red[Epsilon Value]", min_value=0.0, step=0.1,
                                        key="eps_opt", max_value=float(filtered_selected_genes.shape[0]),
                                        value=1.0) if cluster_met == "dbscan" else None
        alg_opt = "auto"

        return optics(selected_genes_dim_red, ms_optics, metric_opt, eps_opt, alg_opt, cluster_met)
    except Exception as e:
        st.error(f"Error in perform_optics_clustering: {e}")
        return None


def clustering_genes(selected_genes_dim_red, filtered_selected_genes, df_names, cluster_name, lat_vec,
                     random_state_samples, adv_col1_gene, ht_expression_path, gene_clust_col3, scatter_2d_name,
                     scatter_3d_name, gene_exp, adv_col2_gene, option_dm, user_session_id, column_cluster_path,
                     selected_genes_dm_path, user_gene_cluster_path, cluster_toggle, chb_gene_clust,e):
    """
    Perform clustering on selected genes (gene-specific).

    Parameters:
    - selected_genes_dim_red (DataFrame): Dimensionally reduced data.
    - filtered_selected_genes (DataFrame): Filtered gene expression data.
    - df_names (list): List of gene names.
    - cluster_name (str): Filename to save clustering results.
    - lat_vec (int): Number of latent vectors.
    - random_state_samples (int): Random state for reproducibility.
    - adv_col1_gene (streamlit.Column): Streamlit column for advanced settings.
    - ht_expression_path (str): Path to save heatmap expression data.
    - gene_clust_col3 (streamlit.Column): Streamlit column for gene clustering settings.
    - scatter_2d_name (str): Filename for 2D scatter plot.
    - scatter_3d_name (str): Filename for 3D scatter plot.
    - gene_exp (streamlit.Column): Streamlit column for gene expression settings.
    - adv_col2_gene (streamlit.Column): Streamlit column for advanced settings.
    - option_dm (str): Selected dimensionality reduction method.
    - user_session_id (str): User identifier.
    - column_cluster_path (str): Path to save column cluster data.
    - selected_genes_dm_path (str): Path to save dimensionally reduced data.
    - user_gene_cluster_path (str): Path to save user gene cluster data.
    - cluster_toggle (bool): Cluster toggle flag.

    Returns:
    - cluster_db (DataFrame): Clustering database.
    """
    try:
        if not cluster_toggle:
            n_clust = None
            labels = None

            cluster_options, clust_index = get_cluster_options(filtered_selected_genes, option_dm, gene=True)

            num_of_clust = get_number_of_clusters(filtered_selected_genes, gene=True)

            option_cluster = gene_exp.selectbox(':blue[Clustering ] :red[Method]',
                                                options=cluster_options, key="pot_clust_genes", index=clust_index)

            if option_cluster in ["K-means", "Gaussian Mixture Model", "Agglomerative Clustering",
                                  'Self-Organizing Map']:
                n_clust = gene_exp.select_slider('Select the number of clusters', options=num_of_clust,
                                                 key="genes_clust_num")
                gene_exp.write(f':blue[Number of cluster:] :red[{str(n_clust)}]')

            if option_cluster == "K-means":
                labels = perform_kmeans_clustering_genes(selected_genes_dim_red, n_clust, random_state_samples,
                                                         adv_col1_gene)
            elif option_cluster == 'Gaussian Mixture Model':
                labels = perform_gmm_clustering_genes(selected_genes_dim_red, n_clust, random_state_samples,
                                                      adv_col1_gene)
            elif option_cluster == 'Agglomerative Clustering':
                labels = perform_hierarchical_clustering_genes(selected_genes_dim_red, n_clust, adv_col1_gene)
            elif option_cluster == 'Self-Organizing Map':
                labels = perform_som_clustering_genes(selected_genes_dim_red, n_clust, random_state_samples,
                                                      adv_col1_gene)
            elif option_cluster == 'DynamicTreeCut':
                labels = perform_dtc_clustering_genes(selected_genes_dim_red, filtered_selected_genes, adv_col1_gene)
            elif option_cluster == 'HDBSCAN':
                labels = perform_hdb_clustering_genes(selected_genes_dim_red, filtered_selected_genes, lat_vec,
                                                      adv_col1_gene)
            elif option_cluster == "OPTICS":
                labels = perform_optics_clustering_genes(selected_genes_dim_red, filtered_selected_genes, adv_col1_gene)
            else:
                st.error(":blue[Select a clustering method and the number of clusters ] "
                         ":red[(Except DynamicCutTree, HDBSCAN and OPTICS)! ] "
                         ":blue[You can apply dimensionality reduction method if you want!]")

            labels = labels + 1
            labels = np.array(["G" + str(number) for number in labels])

            cluster_db = save_disp_genes(df_names, cluster_name, lat_vec, selected_genes_dim_red, labels,
                                         ht_expression_path,
                                         gene_exp, scatter_2d_name, scatter_3d_name, adv_col2_gene, user_session_id,
                                         column_cluster_path, selected_genes_dm_path, chb_gene_clust, e)

            return cluster_db

        else:
            user_cluster_df = pd.read_csv(f'{user_gene_cluster_path}', index_col='Samples')
            user_cluster_df = user_cluster_df[user_cluster_df.index.isin(filtered_selected_genes.columns)]
            user_cluster_df = user_cluster_df.reindex(df_names)
            user_cluster_df.index.rename("Samples", inplace=True)
            user_cluster_df.reset_index(inplace=True)

            labels = user_cluster_df["Cluster"].values
            df_names = user_cluster_df["Samples"].values
            cluster_db = save_disp_genes(df_names, cluster_name, lat_vec, selected_genes_dim_red, labels,
                                         ht_expression_path, gene_exp, scatter_2d_name, scatter_3d_name, adv_col2_gene,
                                         user_session_id, column_cluster_path, selected_genes_dm_path, chb_gene_clust, e)

            return cluster_db

    except Exception as e:
        st.error(f"Error in clustering_genes: {e}")
        return None


def get_cluster_options(filtered_selected_genes, option_dm, gene=False):
    """
    Get cluster options and index based on the filtered selected genes and option_dm.

    Parameters:
    - filtered_selected_genes (DataFrame): Filtered gene expression data.
    - option_dm (str): Selected dimensionality reduction method.
    - gene (bool): Flag for gene-specific clustering.

    Returns:
    - cluster_options (DataFrame): DataFrame with cluster options.
    - clust_index (int): Index of the selected cluster option.
    """
    try:
        if filtered_selected_genes.shape[1] < 20 if gene else filtered_selected_genes.shape[0] < 20:
            cluster_options = pd.DataFrame(
                {'clust_opt': ["K-means", "Gaussian Mixture Model", "Agglomerative Clustering", 'Self-Organizing Map',
                               "HDBSCAN", "OPTICS"]})
            clust_index = 1
        else:
            cluster_options = pd.DataFrame(
                {'clust_opt': ["K-means", "Gaussian Mixture Model", "Agglomerative Clustering", 'Self-Organizing Map',
                               "DynamicTreeCut", "HDBSCAN", "OPTICS"]})
            clust_index = 4
        if option_dm == "Correlation":
            cluster_options = pd.DataFrame(
                {'clust_opt': ["K-means", "Gaussian Mixture Model", "Agglomerative Clustering", "DynamicTreeCut",
                               "HDBSCAN", "OPTICS"]})
            clust_index = 3

        return cluster_options, clust_index
    except Exception as e:
        st.error(f"Error in get_cluster_options (genes): {e}")
        return None, None


def get_number_of_clusters(filtered_selected_genes, gene=False):
    """
    Get the range of number of clusters based on the filtered selected genes.

    Parameters:
    - filtered_selected_genes (DataFrame): Filtered gene expression data.
    - gene (bool): Flag for gene-specific clustering.

    Returns:
    - num_of_clust (list): List of number of clusters.
    """
    try:
        if filtered_selected_genes.shape[1] < 30 if gene else filtered_selected_genes.shape[0] < 30:
            min_num = 2
            max_num = int(filtered_selected_genes.shape[1 if gene else 0])
            num_of_clust = list(range(min_num, max_num + 1))
        else:
            num_of_clust = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
        return num_of_clust
    except Exception as e:
        st.error(f"Error in get_number_of_clusters (genes): {e}")
        return None


def perform_kmeans_clustering_genes(selected_genes_dim_red, n_clust, random_state_samples, adv_col1_gene):
    """
    Perform K-means clustering for genes.

    Parameters:
    - selected_genes_dim_red (DataFrame): Dimensionally reduced data.
    - n_clust (int): Number of clusters.
    - random_state_samples (int): Random state for reproducibility.
    - adv_col1_gene (streamlit.Column): Streamlit column for advanced settings.

    Returns:
    - labels (array): Cluster labels.
    """
    try:
        adv_col1_gene.subheader(":blue[Advanced Settings of K-means (Genes)]", help=dscript.kmeans_adv())
        km_init = adv_col1_gene.selectbox(":blue[K-means ] :red[Initialisation ] :blue[Method]",
                                          options=("random", "k-means++"), index=0, key="km_init_genes")
        km_tol = adv_col1_gene.number_input(":blue[Relative ] :red[Tolerance]", min_value=0.0001, max_value=1.0000,
                                            step=0.0001, value=0.0001, key="km_tol_genes", format="%f")
        mi_km = adv_col1_gene.number_input(":red[Max ] :blue[Iteration]", min_value=200, max_value=2000, step=1,
                                           value=1000, key="km_max_iter_genes")
        km_alg = adv_col1_gene.selectbox(":red[K-means Algorithm] :blue[to Use]", options=("lloyd", "elkan"),
                                         key="gene_kmeans")

        return kmeans_gene(selected_genes_dim_red, n_clust, km_init, km_tol, mi_km, km_alg, random_state_samples)
    except Exception as e:
        st.error(f"Error in perform_kmeans_clustering_genes: {e}")
        return None


def perform_gmm_clustering_genes(selected_genes_dim_red, n_clust, random_state_samples, adv_col1_gene):
    """
    Perform Gaussian Mixture Model clustering for genes.

    Parameters:
    - selected_genes_dim_red (DataFrame): Dimensionally reduced data.
    - n_clust (int): Number of clusters.
    - random_state_samples (int): Random state for reproducibility.
    - adv_col1_gene (streamlit.Column): Streamlit column for advanced settings.

    Returns:
    - labels (array): Cluster labels.
    """
    try:
        adv_col1_gene.subheader(":blue[Advanced Settings of Gaussian Mixture Model (Genes)]",
                                help=dscript.gmm_adv())
        init_method = adv_col1_gene.selectbox(":red[Initialisation ] :blue[Methods]",
                                              options=("kmeans", "k-means++", "random", "random_from_data"),
                                              index=1, key="in_met_gmm_genes")
        covariance_type = adv_col1_gene.selectbox(":red[Convariance] :blue[Type]", options=("full", "tied"),
                                                  index=0, key="cov_typ_gene")
        tol_gmm = adv_col1_gene.select_slider(":blue[Convergence ] :red[Threshold]", options=(
            0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1),
                                              value=0.0001, key="tol_gmm_ss_genes")
        reg_cov_gmm = adv_col1_gene.select_slider(':red[Non-negative ] :blue[Regularization]',
                                                  options=(1e-3, 1e-4, 1e-5, 1e-6, 1e-7), value=1e-5,
                                                  key='reg_covar_gmm_genes')

        return gmm_gene(selected_genes_dim_red, n_clust, init_method, covariance_type, tol_gmm, random_state_samples,
                        reg_cov_gmm)
    except Exception as e:
        st.error(f"Error in perform_gmm_clustering_genes: {e}")
        return None


def perform_hierarchical_clustering_genes(selected_genes_dim_red, n_clust, adv_col1_gene):
    """
    Perform hierarchical clustering for genes.

    Parameters:
    - selected_genes_dim_red (DataFrame): Dimensionally reduced data.
    - n_clust (int): Number of clusters.
    - adv_col1_gene (streamlit.Column): Streamlit column for advanced settings.

    Returns:
    - labels (array): Cluster labels.
    """
    try:
        adv_col1_gene.subheader(":blue[Advanced Settings of Agglomerative Clustering (Genes)]",
                                help=dscript.agglo_adv())
        link = adv_col1_gene.selectbox(":red[Linkage Method ] :blue[of Agglomerative Clustering]",
                                       options=("ward", "complete", "average", 'single'), index=0,
                                       key="link_hc_genes")
        if link == "ward":
            hc_met = "euclidean"
        else:
            hc_met = adv_col1_gene.selectbox("Metric of Agglomerative Clustering",
                                             options=('euclidean', 'l1', 'l2', 'manhattan', 'cosine'),
                                             index=0, key="metric_hc_genes")

        return hierarchical_clustering_gene(selected_genes_dim_red, n_clust, link, hc_met)
    except Exception as e:
        st.error(f"Error in perform_hierarchical_clustering_genes: {e}")
        return None


def perform_som_clustering_genes(selected_genes_dim_red, n_clust, random_state_samples, adv_col1_gene):
    """
    Perform Self-Organizing Map clustering for genes.

    Parameters:
    - selected_genes_dim_red (DataFrame): Dimensionally reduced data.
    - n_clust (int): Number of clusters.
    - random_state_samples (int): Random state for reproducibility.
    - adv_col1_gene (streamlit.Column): Streamlit column for advanced settings.

    Returns:
    - labels (array): Cluster labels.
    """
    try:
        adv_col1_gene.subheader(":blue[Advanced Settings of Self-Organizing Map (Genes)]", help=dscript.som_adv())
        neg_som = adv_col1_gene.selectbox(":red[Neighbours ] :blue[Function Type]",
                                          options=('gaussian', 'mexican_hat', 'bubble', 'triangle'), index=0,
                                          key="neg_som_genes")
        sig_som = adv_col1_gene.slider(":red[Sigma ] :blue[Value]", key="sig_som_genes", min_value=0.10,
                                       max_value=50.00, value=1.50, step=0.10)
        lr_som = adv_col1_gene.slider(":red[Learning] :blue[ Rate]", key="lr_som_gene", min_value=0.10,
                                      max_value=5.00, value=0.50)
        epoch = adv_col1_gene.number_input(":blue[Number of ] :red[Epochs]", min_value=1, max_value=2000,
                                           value=500, step=1,
                                           key="som_epo_genes", format="%d")
        act_dist = adv_col1_gene.selectbox(":blue[Choose the ] :red[Activation Distance]", key='act_dist_gene',
                                           options=("euclidean", "cosine", "manhattan", "chebyshev"), index=0)

        topology = adv_col1_gene.selectbox(":blue[Select the ] :red[Topology of the Map]", key="tplgy_gene",
                                           options=("rectangular", "hexagonal"), index=0)

        return som_gene(selected_genes_dim_red, n_clust, neg_som, sig_som, lr_som, epoch, random_state_samples,
                        act_dist, topology)
    except Exception as e:
        st.error(f"Error in perform_som_clustering_genes: {e}")
        return None


def perform_dtc_clustering_genes(selected_genes_dim_red, filtered_selected_genes, adv_col1_gene):
    """
    Perform Dynamic Tree Cut clustering for genes.

    Parameters:
    - selected_genes_dim_red (DataFrame): Dimensionally reduced data.
    - filtered_selected_genes (DataFrame): Filtered gene expression data.
    - adv_col1_gene (streamlit.Column): Streamlit column for advanced settings.

    Returns:
    - labels (array): Cluster labels.
    """
    try:
        adv_col1_gene.subheader(":blue[Advanced Settings of DynamicTreeCut (Genes)]", help=dscript.dtc_adv())
        link = adv_col1_gene.selectbox(":red[Linkage ] :blue[Method]", options=(
            "single", "complete", "average", "weighted", "centroid", "median", "ward"), index=6,
                                       key="link_dct_genes")
        distance = adv_col1_gene.selectbox(":red[Distance ] :blue[Metric]",
                                           options=("braycurtis", "canberra", "chebyshev", "cityblock",
                                                    "correlation", "cosine", "euclidean",
                                                    "matching", "minkowski", "seuclidean",
                                                    "sqeuclidean"),
                                           key="dct_dist_1_genes" if selected_genes_dim_red.shape[1] >
                                           selected_genes_dim_red.shape[0] else
                                           "dct_dist_2_genes", index=6)

        min_clust_size_dtc = adv_col1_gene.number_input(":red[Minimum ] :blue[Cluster Size]", min_value=2, step=1,
                                                        value=2,
                                                        max_value=int(filtered_selected_genes.shape[1] / 2),
                                                        key="min_clust_gene")

        labels = dtc_gene(selected_genes_dim_red, link, distance, min_clust_size_dtc)
        return labels["labels"] - 1
    except Exception as e:
        st.error(f"Error in perform_dtc_clustering_genes: {e}")
        return None


def perform_hdb_clustering_genes(selected_genes_dim_red, filtered_selected_genes, lat_vec, adv_col1_gene):
    """
    Perform HDBSCAN clustering for genes.

    Parameters:
    - selected_genes_dim_red (DataFrame): Dimensionally reduced data.
    - filtered_selected_genes (DataFrame): Filtered gene expression data.
    - lat_vec (int): Number of latent vectors.
    - adv_col1_gene (streamlit.Column): Streamlit column for advanced settings.

    Returns:
    - labels (array): Cluster labels.
    """
    try:
        adv_col1_gene.subheader(":blue[Advanced Settings of HDBSCAN (Genes)]", help=dscript.hdb_adv())
        alfa = adv_col1_gene.number_input(":red[Alpha ] :blue[Value]", min_value=0.1, max_value=5.0, value=1.0,
                                          step=0.1, format='%f',
                                          key="alfa_dtc_genes")
        hdb_met = adv_col1_gene.selectbox(":red[Metric ] :blue[Method]", options=(
            "braycurtis", "canberra", "chebyshev", "cityblock", "euclidean", "hamming",
            "infinity", "manhattan") if lat_vec == 3 else (
            "braycurtis", "canberra", "chebyshev", "cityblock", "euclidean", "haversine", "hamming",
            "infinity", "l1", "l2", "manhattan", "matching"),
                                          key="hdb_metric_2d_genes" if lat_vec == 3 else "hdb_metric_3d_genes", index=4)
        mcs_hdb = adv_col1_gene.number_input(":red[Minimum Number of Samples ] :blue[in Each Cluster]",
                                             min_value=1,
                                             max_value=int(filtered_selected_genes.shape[1] / 2), value=2,
                                             key="mcs_key_genes")
        min_sam_hdb = adv_col1_gene.number_input(":red[Minimum Number of Samples ] :blue[in Each Cluster]",
                                                 min_value=2,
                                                 max_value=int(filtered_selected_genes.shape[1] / 2), value=2,
                                                 key="min_sam_genes")
        eps = adv_col1_gene.number_input(":blue[Cluster Selection ] :red[Epsilon Value]", min_value=0.0,
                                         max_value=float(filtered_selected_genes.shape[1] / 2),
                                         value=1.0, step=0.1,
                                         key="eps_slider_gene", format="%f")
        clust_meth = adv_col1_gene.selectbox(":blue[Select the ] :red[Cluster Selection Method]",
                                             key='clust_sel_gene', options=("eom", "leaf"), index=0)

        return hdb_gene(selected_genes_dim_red, alfa, hdb_met, mcs_hdb, min_sam_hdb, eps, clust_meth)
    except Exception as e:
        st.error(f"Error in perform_hdb_clustering_genes: {e}")
        return None


def perform_optics_clustering_genes(selected_genes_dim_red, filtered_selected_genes, adv_col1_gene):
    """
    Perform OPTICS clustering for genes.

    Parameters:
    - selected_genes_dim_red (DataFrame): Dimensionally reduced data.
    - filtered_selected_genes (DataFrame): Filtered gene expression data.
    - adv_col1_gene (streamlit.Column): Streamlit column for advanced settings.

    Returns:
    - labels (array): Cluster labels.
    """
    try:
        adv_col1_gene.subheader(":blue[Advanced Settings of OPTICS (Genes)]", help=dscript.opt_adv())
        ms_optics = adv_col1_gene.number_input(":red[Minimum ] :blue[Samples]", min_value=2, step=1, value=5,
                                               key="ms_optics_genes",
                                               max_value=int(filtered_selected_genes.shape[1] / 2))
        metric_opt = adv_col1_gene.selectbox(":red[Distance ] :blue[Metric]", options=(
            "minkowski", "cityblock", "cosine", "euclidean", "l1", "l2", "manhattan", "braycurtis",
            "canberra", "chebyshev", "correlation", "sqeuclidean"), key="metric_opt_genes", index=3)
        cluster_met = adv_col1_gene.selectbox(":red[Initialisation ] :blue[Method]", options=("xi", "dbscan"),
                                              key="clust_met_gene", index=1)
        eps_opt = adv_col1_gene.number_input(":blue[Cluster Selection ] :red[Epsilon Value]", min_value=0.0,
                                             step=0.1, key="eps_opt_genes",
                                             max_value=float(filtered_selected_genes.shape[1]),
                                             value=1.0) if cluster_met == "dbscan" else None
        alg_opt = "auto"

        return optics_gene(selected_genes_dim_red, ms_optics, metric_opt, eps_opt, alg_opt, cluster_met)
    except Exception as e:
        st.error(f"Error in perform_optics_clustering_genes: {e}")
        return None
