import streamlit as st
from umap.umap_ import UMAP
from sklearn.manifold import TSNE, MDS
from sklearn.decomposition import PCA, NMF
from sklearn.cluster import KMeans, OPTICS, AgglomerativeClustering
from sklearn.mixture import GaussianMixture
from minisom import MiniSom
from scipy.cluster.hierarchy import linkage
from dynamicTreeCut import cutreeHybrid
from scipy.spatial.distance import pdist
from sklearn.cluster import HDBSCAN
from bignmf.models.jnmf.integrative import IntegrativeJnmf
import pandas as pd
import numpy as np
import functools
from Scripts.user_interaction import get_session_id
import inspect
import phate


# Reusable function for logging arguments
def create_logger(log_file_template, log_msg_template):
    def log_args(func):
        sig = inspect.signature(func)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            bound = sig.bind(*args, **kwargs)
            bound.apply_defaults()
            arg_list = list(bound.arguments.items())
            arg_strs = [f"{arg}={value}" for arg, value in arg_list[1:]]  # ignore the first argument
            user_id = get_session_id()
            log_file = log_file_template.format(user_id=user_id)
            with open(log_file, "w") as f:
                print(log_msg_template.format(func_name=func.__name__, args=", ".join(arg_strs)), file=f)
            return func(*args, **kwargs)

        return wrapper

    return log_args


log_args_clust = create_logger(
    "result_data/result_data_{user_id}/clust_info_{user_id}.txt",
    'Calling "{func_name}" clustering(samples) with settings: {args}'
)

log_args_dimred = create_logger(
    "result_data/result_data_{user_id}/dim_red_info_{user_id}.txt",
    'Calling "{func_name}" (samples) with settings: {args}'
)

log_args_clust_gene = create_logger(
    "result_data/result_data_{user_id}/clust_info_gene_{user_id}.txt",
    'Calling "{func_name}" clustering(genes) with settings: {args}'
)

log_args_dimred_gene = create_logger(
    "result_data/result_data_{user_id}/dim_red_info_gene_{user_id}.txt",
    'Calling "{func_name}" (genes) with settings: {args}'
)


def nmf_data(nmf_name):
    """
    Load and preprocess NMF data.

    Parameters:
    nmf_name (str): Filename of the NMF data.

    Returns:
    tuple: Processed data dictionary, sample names, and name order.
    """
    selected_genes_nmf = pd.read_csv(nmf_name, index_col=0, header=0)
    data_dict = {"sim1": selected_genes_nmf}
    sample_names1 = selected_genes_nmf.transpose().index
    sample_names2 = selected_genes_nmf.index
    name_order = list(sample_names1)
    return data_dict, sample_names1, sample_names2, name_order


@log_args_clust
def nmf_clust(_data_dict, k, iter_nmf, tri_nmf, lamb_nmf):
    """
    Perform NMF clustering on the provided data dictionary.

    Parameters:
    _data_dict (dict): Data dictionary for NMF.
    k (int): Number of clusters.
    iter_nmf (int): Number of iterations for NMF.
    tri_nmf (int): Number of trials for NMF.
    lamb_nmf (float): Regularization parameter for NMF.

    Returns:
    model: Trained NMF model.
    clust_nmf: Clustering result.
    """
    try:
        model = IntegrativeJnmf(_data_dict, k, lamb_nmf)
        model.run(tri_nmf, iter_nmf)
        model.cluster_data()
        clust_nmf = model.h_cluster
        model.calc_consensus_matrices()
        return model, clust_nmf
    except Exception as e:
        st.error(f"Error during NMF clustering: {e}")
        return None, None


@log_args_dimred
@st.cache_data(show_spinner="MDS Calculating...", max_entries=10, ttl=1800)
def MDS_dimension_reduction(data, lat_vec, mi, eps_mds, rs_mds, metric_mds):
    """
    Perform MDS dimension reduction on the provided data.

    Parameters:
    data (DataFrame): Input data for MDS.
    lat_vec (int): Number of latent vectors (components) for MDS.
    mi (int): Maximum number of iterations for MDS.
    eps_mds (float): Stopping criterion for MDS.
    rs_mds (int): Random state for reproducibility.
    metric_mds (bool): Whether to use metric MDS.

    Returns:
    DataFrame: Reduced data after MDS.
    """
    try:
        mds = MDS(n_components=lat_vec, max_iter=mi, eps=eps_mds, random_state=rs_mds, metric=metric_mds)
        reduced_data = mds.fit_transform(data)
        return reduced_data
    except Exception as e:
        st.error(f"Error during MDS dimension reduction: {e}")
        return None


@log_args_dimred
@st.cache_data(show_spinner="NMF Calculating...", max_entries=10, ttl=1800)
def NMF_dimension_reduction(data, lat_vec, max_it, tolerance, init_alg, b_loss, rs_nmf, solver, shuffle):
    """
    Perform NMF dimension reduction on the provided data.

    Parameters:
    data (DataFrame): Input data for NMF.
    lat_vec (int): Number of latent vectors (components) for NMF.
    max_it (int): Maximum number of iterations for NMF.
    tolerance (float): Tolerance for NMF.
    init_alg (str): Initialization algorithm for NMF.
    b_loss (str): Beta loss function for NMF.
    rs_nmf (int): Random state for reproducibility.
    solver (str): Solver for NMF.
    shuffle (bool): Whether to shuffle the data.

    Returns:
    DataFrame: Reduced data after NMF.
    """
    try:
        nmf = NMF(n_components=lat_vec, max_iter=max_it, tol=tolerance, random_state=rs_nmf, shuffle=shuffle,
                  init=init_alg,
                  solver=solver, beta_loss=b_loss, alpha_W=0, alpha_H="same")
        reduced_data = nmf.fit_transform(data)
        return reduced_data
    except Exception as e:
        st.error(f"Error during NMF dimension reduction: {e}")
        return None


@log_args_dimred
@st.cache_data(show_spinner="PCA Calculating...", max_entries=10, ttl=1800)
def PCA_dimension_reduction(data, lat_vec, tolerance, rs_pca, whiten, svd_solver_app):
    """
    Perform PCA dimension reduction on the provided data.

    Parameters:
    data (DataFrame): Input data for PCA.
    lat_vec (int): Number of latent vectors (components) for PCA.
    tolerance (float): Tolerance for PCA.
    rs_pca (int): Random state for reproducibility.
    whiten (bool): Whether to whiten the components.

    Returns:
    DataFrame: Reduced data after PCA.
    """
    try:
        pca = PCA(n_components=lat_vec, random_state=rs_pca, tol=tolerance, power_iteration_normalizer="auto",
                  whiten=whiten, svd_solver=svd_solver_app)
        reduced_data = pca.fit_transform(data)
        return reduced_data
    except Exception as e:
        st.error(f"Error during PCA dimension reduction: {e}")
        return None


@log_args_dimred
@st.cache_data(show_spinner="UMAP Calculating...", max_entries=10, ttl=1800)
def UMAP_dimension_reduction(data, lat_vec, min_dis, n_neigh, rs_umap, metric):
    """
    Perform UMAP dimension reduction on the provided data.

    Parameters:
    data (DataFrame): Input data for UMAP.
    lat_vec (int): Number of latent vectors (components) for UMAP.
    min_dis (float): Minimum distance for UMAP.
    n_neigh (int): Number of neighbors for UMAP.
    rs_umap (int): Random state for reproducibility.
    metric (str): Metric for UMAP.

    Returns:
    DataFrame: Reduced data after UMAP.
    """
    try:
        reduced_data = UMAP(n_components=lat_vec, min_dist=min_dis, n_neighbors=n_neigh,
                            random_state=rs_umap, metric=metric).fit_transform(data)
        return reduced_data
    except Exception as e:
        st.error(f"Error during UMAP dimension reduction: {e}")
        return None


@log_args_dimred
@st.cache_data(show_spinner="PHATE Calculating...", max_entries=10, ttl=1800)
def phate_dimension_reduction(data, lat_vec, n_neight, random_sample, n_pca):
    """
    Perform PHATE dimension reduction on the provided data.

    Parameters:
    data (DataFrame): Input data for PHATE.
    lat_vec (int): Number of latent vectors (components) for PHATE.
    n_neight (int): Number of neighbors for PHATE.
    random_sample (int): Random state for reproducibility.
    n_pca (int): Number of principal components for PHATE.

    Returns:
    DataFrame: Reduced data after PHATE.
    """
    try:
        phate_run = phate.PHATE(n_components=lat_vec, knn=n_neight, n_pca=n_pca, random_state=random_sample)
        reduced_data = phate_run.fit_transform(data)
        return reduced_data
    except Exception as e:
        st.error(f"Error during PHATE dimension reduction: {e}")
        return None


@log_args_dimred
@st.cache_data(show_spinner="t-SNE Calculating...", max_entries=10, ttl=1800)
def t_SNE_dimension_reduction(data, lat_vec, perp, lr, niter, met, rs_tsne, init_met):
    """
    Perform t-SNE dimension reduction on the provided data.

    Parameters:
    data (DataFrame): Input data for t-SNE.
    lat_vec (int): Number of latent vectors (components) for t-SNE.
    perp (float): Perplexity for t-SNE.
    lr (float): Learning rate for t-SNE.
    niter (int): Number of iterations for t-SNE.
    met (str): Metric for t-SNE.
    rs_tsne (int): Random state for reproducibility.
    init_met (str): Initialization method for t-SNE.

    Returns:
    DataFrame: Reduced data after t-SNE.
    """
    try:
        tsne = TSNE(n_components=lat_vec, perplexity=perp, learning_rate=lr, random_state=rs_tsne, n_iter=niter,
                    metric=met,
                    init=init_met)
        reduced_data = tsne.fit_transform(data)
        return reduced_data
    except Exception as e:
        st.error(f"Error during t-SNE dimension reduction: {e}")
        return None


@log_args_clust
@st.cache_data(show_spinner="Agglomerative Clustering is Calculating", max_entries=10, ttl=1800)
def hierarchical_clustering(data, k, link, hc_met):
    """
    Perform hierarchical clustering on the provided data.

    Parameters:
    data (DataFrame): Input data for hierarchical clustering.
    k (int): Number of clusters.
    link (str): Linkage criterion for hierarchical clustering.
    hc_met (str): Metric for hierarchical clustering.

    Returns:
    array: Cluster labels.
    """
    try:
        agg_clustering = AgglomerativeClustering(n_clusters=k, linkage=link, compute_distances=True, metric=hc_met)
        labels_hc = agg_clustering.fit_predict(data)
        return labels_hc
    except Exception as e:
        st.error(f"Error during hierarchical clustering: {e}")
        return None


@log_args_clust
@st.cache_data(show_spinner="Gaussian Mixture Model Clustering is Calculating", max_entries=10, ttl=1800)
def gmm(data, k, init_method, covariance_type, tol_gmm, rs_gmm, reg_cov_gmm):
    """
    Perform Gaussian Mixture Model clustering on the provided data.

    Parameters:
    data (DataFrame): Input data for GMM.
    k (int): Number of clusters.
    init_method (str): Initialization method for GMM.
    covariance_type (str): Covariance type for GMM.
    tol_gmm (float): Tolerance for GMM.
    rs_gmm (int): Random state for reproducibility.
    reg_cov_gmm (float): Regularization for covariance matrix in GMM.

    Returns:
    array: Cluster labels.
    """
    try:
        gmm_ = GaussianMixture(n_components=k, init_params=init_method, tol=tol_gmm, covariance_type=covariance_type,
                               n_init=k, random_state=rs_gmm, reg_covar=reg_cov_gmm, max_iter=100)
        gmm_.fit_predict(data)
        return gmm_.predict(data)
    except Exception as e:
        st.error(f"Error during GMM clustering: {e}")
        return None


@log_args_clust
@st.cache_data(show_spinner="K-means Clustering is Calculating", max_entries=10, ttl=1800)
def kmeans(data, k, km_init, km_tol, mi_km, km_alg, rs_km):
    """
    Perform K-means clustering on the provided data.

    Parameters:
    data (DataFrame): Input data for K-means clustering.
    k (int): Number of clusters.
    km_init (str): Initialization method for K-means.
    km_tol (float): Tolerance for K-means.
    mi_km (int): Maximum number of iterations for K-means.
    km_alg (str): Algorithm for K-means.
    rs_km (int): Random state for reproducibility.

    Returns:
    array: Cluster labels.
    """
    try:
        kmeans_ = KMeans(n_clusters=k, init=km_init, tol=km_tol, n_init="auto", max_iter=mi_km, random_state=rs_km,
                         algorithm=km_alg)
        kmeans_.fit_predict(data)
        return kmeans_.labels_
    except Exception as e:
        st.error(f"Error during K-means clustering: {e}")
        return None


@log_args_clust
@st.cache_data(show_spinner="Self-Organizing Map is Calculating", max_entries=10, ttl=1800)
def som(data, k, neg_som, sig_som, lr_som, epoch, rs_som, act_dist, tplgy):
    """
    Perform Self-Organizing Map clustering on the provided data.

    Parameters:
    data (DataFrame): Input data for SOM.
    k (int): Number of clusters.
    neg_som (str): Neighborhood function for SOM.
    sig_som (float): Sigma for SOM.
    lr_som (float): Learning rate for SOM.
    epoch (int): Number of epochs for training SOM.
    rs_som (int): Random state for reproducibility.
    act_dist (str): Activation distance for SOM.
    tplgy (str): Topology for SOM.

    Returns:
    array: Cluster labels.
    """
    try:
        som_shape = (1, k)
        som_cluster = MiniSom(som_shape[0], som_shape[1], data.shape[1], random_seed=rs_som, sigma=sig_som,
                              learning_rate=lr_som, neighborhood_function=neg_som, activation_distance=act_dist,
                              topology=tplgy)
        som_cluster.random_weights_init(data)
        som_cluster.train_batch(data, epoch, verbose=True)
        winner_coordinates = np.array([som_cluster.winner(x) for x in data]).T
        cluster_index = np.ravel_multi_index(winner_coordinates, som_shape)
        return cluster_index
    except Exception as e:
        st.error(f"Error during SOM clustering: {e}")
        return None


@log_args_clust
@st.cache_data(show_spinner="Dynamic Cut Tree Clustering is Calculating", max_entries=10, ttl=1800)
def dtc(data, link, distance, min_clust_size_dtc):
    """
    Perform Dynamic Tree Cut clustering on the provided data.

    Parameters:
    data (DataFrame): Input data for Dynamic Tree Cut.
    link (str): Linkage criterion for Dynamic Tree Cut.
    distance (str): Distance metric for Dynamic Tree Cut.
    min_clust_size_dtc (int): Minimum cluster size for Dynamic Tree Cut.

    Returns:
    array: Cluster labels.
    """
    try:
        dist_dtc = pdist(data, distance)
        link_dtc = linkage(dist_dtc, link)
        clusters = cutreeHybrid(link_dtc, dist_dtc, minClusterSize=min_clust_size_dtc)
        return clusters
    except Exception as e:
        st.error(f"Error during Dynamic Tree Cut clustering: {e}")
        return None


@log_args_clust
@st.cache_data(show_spinner="HDBSCAN Clustering is Calculating", max_entries=10, ttl=1800)
def hdb(data, alfa, hdb_met, mcs_hdb, min_sam_hdb, eps, clust_meth):
    """
    Perform HDBSCAN clustering on the provided data.

    Parameters:
    data (DataFrame): Input data for HDBSCAN.
    alfa (float): Alpha parameter for HDBSCAN.
    hdb_met (str): Metric for HDBSCAN.
    mcs_hdb (int): Minimum cluster size for HDBSCAN.
    min_sam_hdb (int): Minimum samples for HDBSCAN.
    eps (float): Epsilon parameter for HDBSCAN.
    clust_meth (str): Cluster selection method for HDBSCAN.

    Returns:
    array: Cluster labels.
    """
    try:
        clusters = HDBSCAN(algorithm='auto', alpha=alfa, leaf_size=40, metric=hdb_met, min_cluster_size=mcs_hdb,
                           min_samples=min_sam_hdb, cluster_selection_method=clust_meth,
                           allow_single_cluster=False, cluster_selection_epsilon=eps).fit_predict(data)
        return clusters
    except Exception as e:
        st.error(f"Error during HDBSCAN clustering: {e}")
        return None


@log_args_clust
@st.cache_data(show_spinner="OPTICS Clustering is Calculating", max_entries=10, ttl=1800)
def optics(data, ms_optics, metric_opt, eps_opt, alg_opt, cluster_met):
    """
    Perform OPTICS clustering on the provided data.

    Parameters:
    data (DataFrame): Input data for OPTICS.
    ms_optics (int): Minimum samples for OPTICS.
    metric_opt (str): Metric for OPTICS.
    eps_opt (float): Epsilon parameter for OPTICS.
    alg_opt (str): Algorithm for OPTICS.
    cluster_met (str): Cluster selection method for OPTICS.

    Returns:
    array: Cluster labels.
    """
    try:
        clustering = OPTICS(min_samples=ms_optics, cluster_method=cluster_met, metric=metric_opt,
                            eps=eps_opt, algorithm=alg_opt).fit_predict(data)
        return clustering
    except Exception as e:
        st.error(f"Error during OPTICS clustering: {e}")
        return None


@log_args_dimred_gene
@st.cache_data(show_spinner="MDS Calculating...", max_entries=10, ttl=1800)
def MDS_dimension_reduction_gene(data, lat_vec, mi, eps_mds, rs_mds, metric_mds):
    """
    Perform MDS dimension reduction on gene data.

    Parameters:
    data (DataFrame): Input data for MDS.
    lat_vec (int): Number of latent vectors (components) for MDS.
    mi (int): Maximum number of iterations for MDS.
    eps_mds (float): Stopping criterion for MDS.
    rs_mds (int): Random state for reproducibility.
    metric_mds (bool): Whether to use metric MDS.

    Returns:
    DataFrame: Reduced data after MDS.
    """
    try:
        mds = MDS(n_components=lat_vec, max_iter=mi, eps=eps_mds, random_state=rs_mds, metric=metric_mds)
        reduced_data = mds.fit_transform(data)
        return reduced_data
    except Exception as e:
        st.error(f"Error during MDS dimension reduction: {e}")
        return None


@log_args_dimred_gene
@st.cache_data(show_spinner="NMF Calculating...", max_entries=10, ttl=1800)
def NMF_dimension_reduction_gene(data, lat_vec, max_it, tolerance, init_alg, b_loss, rs_nmf, solver, shuffle):
    """
    Perform NMF dimension reduction on gene data.

    Parameters:
    data (DataFrame): Input data for NMF.
    lat_vec (int): Number of latent vectors (components) for NMF.
    max_it (int): Maximum number of iterations for NMF.
    tolerance (float): Tolerance for NMF.
    init_alg (str): Initialization algorithm for NMF.
    b_loss (str): Beta loss function for NMF.
    rs_nmf (int): Random state for reproducibility.
    solver (str): Solver for NMF.
    shuffle (bool): Whether to shuffle the data.

    Returns:
    DataFrame: Reduced data after NMF.
    """
    try:
        nmf = NMF(n_components=lat_vec, max_iter=max_it, tol=tolerance, random_state=rs_nmf, shuffle=shuffle,
                  init=init_alg, solver=solver, beta_loss=b_loss, alpha_W=0, alpha_H="same")
        reduced_data = nmf.fit_transform(data)
        return reduced_data
    except Exception as e:
        st.error(f"Error during NMF dimension reduction: {e}")
        return None


@log_args_dimred_gene
@st.cache_data(show_spinner="PCA Calculating...", max_entries=10, ttl=1800)
def PCA_dimension_reduction_gene(data, lat_vec, tolerance, rs_pca, whiten):
    """
    Perform PCA dimension reduction on gene data.

    Parameters:
    data (DataFrame): Input data for PCA.
    lat_vec (int): Number of latent vectors (components) for PCA.
    tolerance (float): Tolerance for PCA.
    rs_pca (int): Random state for reproducibility.
    whiten (bool): Whether to whiten the components.

    Returns:
    DataFrame: Reduced data after PCA.
    """
    try:
        pca = PCA(n_components=lat_vec, random_state=rs_pca, tol=tolerance, power_iteration_normalizer="auto",
                  whiten=whiten)
        reduced_data = pca.fit_transform(data)
        return reduced_data
    except Exception as e:
        st.error(f"Error during PCA dimension reduction: {e}")
        return None


@log_args_dimred_gene
@st.cache_data(show_spinner="UMAP Calculating...", max_entries=10, ttl=1800)
def UMAP_dimension_reduction_gene(data, lat_vec, min_dis, n_neigh, rs_umap, metric):
    """
    Perform UMAP dimension reduction on gene data.

    Parameters:
    data (DataFrame): Input data for UMAP.
    lat_vec (int): Number of latent vectors (components) for UMAP.
    min_dis (float): Minimum distance for UMAP.
    n_neigh (int): Number of neighbors for UMAP.
    rs_umap (int): Random state for reproducibility.
    metric (str): Metric for UMAP.

    Returns:
    DataFrame: Reduced data after UMAP.
    """
    try:
        reduced_data = UMAP(n_components=lat_vec, min_dist=min_dis, n_neighbors=n_neigh,
                            random_state=rs_umap, metric=metric).fit_transform(data)
        return reduced_data
    except Exception as e:
        st.error(f"Error during UMAP dimension reduction: {e}")
        return None


@log_args_dimred_gene
@st.cache_data(show_spinner="t-SNE Calculating...", max_entries=10, ttl=1800)
def t_SNE_dimension_reduction_gene(data, lat_vec, perp, lr, niter, met, rs_tsne, init_met):
    """
    Perform t-SNE dimension reduction on gene data.

    Parameters:
    data (DataFrame): Input data for t-SNE.
    lat_vec (int): Number of latent vectors (components) for t-SNE.
    perp (float): Perplexity for t-SNE.
    lr (float): Learning rate for t-SNE.
    niter (int): Number of iterations for t-SNE.
    met (str): Metric for t-SNE.
    rs_tsne (int): Random state for reproducibility.
    init_met (str): Initialization method for t-SNE.

    Returns:
    DataFrame: Reduced data after t-SNE.
    """
    try:
        tsne = TSNE(n_components=lat_vec, perplexity=perp, learning_rate=lr, random_state=rs_tsne, n_iter=niter,
                    metric=met,
                    init=init_met)
        reduced_data = tsne.fit_transform(data)
        return reduced_data
    except Exception as e:
        st.error(f"Error during t-SNE dimension reduction: {e}")
        return None


@log_args_clust_gene
@st.cache_data(show_spinner="Agglomerative Clustering is Calculating", max_entries=10, ttl=1800)
def hierarchical_clustering_gene(data, k, link, hc_met):
    """
    Perform hierarchical clustering on gene data.

    Parameters:
    data (DataFrame): Input data for hierarchical clustering.
    k (int): Number of clusters.
    link (str): Linkage criterion for hierarchical clustering.
    hc_met (str): Metric for hierarchical clustering.

    Returns:
    array: Cluster labels.
    """
    try:
        agg_clustering = AgglomerativeClustering(n_clusters=k, linkage=link, compute_distances=True, metric=hc_met)
        labels_hc = agg_clustering.fit_predict(data)
        return labels_hc
    except Exception as e:
        st.error(f"Error during hierarchical clustering: {e}")
        return None


@log_args_clust_gene
@st.cache_data(show_spinner="Gaussian Mixture Model Clustering is Calculating", max_entries=10, ttl=1800)
def gmm_gene(data, k, init_method, covariance_type, tol_gmm, rs_gmm, reg_cov_gmm):
    """
    Perform Gaussian Mixture Model clustering on gene data.

    Parameters:
    data (DataFrame): Input data for GMM.
    k (int): Number of clusters.
    init_method (str): Initialization method for GMM.
    covariance_type (str): Covariance type for GMM.
    tol_gmm (float): Tolerance for GMM.
    rs_gmm (int): Random state for reproducibility.
    reg_cov_gmm (float): Regularization for covariance matrix in GMM.

    Returns:
    array: Cluster labels.
    """
    try:
        gmm_ = GaussianMixture(n_components=k, init_params=init_method, tol=tol_gmm, covariance_type=covariance_type,
                               n_init=k, random_state=rs_gmm, reg_covar=reg_cov_gmm, max_iter=100)
        gmm_.fit_predict(data)
        return gmm_.predict(data)
    except Exception as e:
        st.error(f"Error during GMM clustering: {e}")
        return None


@log_args_clust_gene
@st.cache_data(show_spinner="K-means Clustering is Calculating", max_entries=10, ttl=1800)
def kmeans_gene(data, k, km_init, km_tol, mi_km, km_alg, rs_km):
    """
    Perform K-means clustering on gene data.

    Parameters:
    data (DataFrame): Input data for K-means clustering.
    k (int): Number of clusters.
    km_init (str): Initialization method for K-means.
    km_tol (float): Tolerance for K-means.
    mi_km (int): Maximum number of iterations for K-means.
    km_alg (str): Algorithm for K-means.
    rs_km (int): Random state for reproducibility.

    Returns:
    array: Cluster labels.
    """
    try:
        kmeans_ = KMeans(n_clusters=k, init=km_init, tol=km_tol, n_init="auto", max_iter=mi_km, random_state=rs_km,
                         algorithm=km_alg)
        kmeans_.fit_predict(data)
        return kmeans_.labels_
    except Exception as e:
        st.error(f"Error during K-means clustering: {e}")
        return None


@log_args_clust_gene
@st.cache_data(show_spinner="Self-Organizing Map is Calculating", max_entries=10, ttl=1800)
def som_gene(data, k, neg_som, sig_som, lr_som, epoch, rs_som, act_dist, tplgy):
    """
    Perform Self-Organizing Map clustering on gene data.

    Parameters:
    data (DataFrame): Input data for SOM.
    k (int): Number of clusters.
    neg_som (str): Neighborhood function for SOM.
    sig_som (float): Sigma for SOM.
    lr_som (float): Learning rate for SOM.
    epoch (int): Number of epochs for training SOM.
    rs_som (int): Random state for reproducibility.
    act_dist (str): Activation distance for SOM.
    tplgy (str): Topology for SOM.

    Returns:
    array: Cluster labels.
    """
    try:
        som_shape = (1, k)
        som_cluster = MiniSom(som_shape[0], som_shape[1], data.shape[1], random_seed=rs_som, sigma=sig_som,
                              learning_rate=lr_som, neighborhood_function=neg_som, activation_distance=act_dist,
                              topology=tplgy)
        som_cluster.random_weights_init(data)
        som_cluster.train_batch(data, epoch, verbose=True)
        winner_coordinates = np.array([som_cluster.winner(x) for x in data]).T
        cluster_index = np.ravel_multi_index(winner_coordinates, som_shape)
        return cluster_index
    except Exception as e:
        st.error(f"Error during SOM clustering: {e}")
        return None


@log_args_clust_gene
@st.cache_data(show_spinner="Dynamic Cut Tree Clustering is Calculating", max_entries=10, ttl=1800)
def dtc_gene(data, link, distance, min_clust_size_dtc):
    """
    Perform Dynamic Tree Cut clustering on gene data.

    Parameters:
    data (DataFrame): Input data for Dynamic Tree Cut.
    link (str): Linkage criterion for Dynamic Tree Cut.
    distance (str): Distance metric for Dynamic Tree Cut.
    min_clust_size_dtc (int): Minimum cluster size for Dynamic Tree Cut.

    Returns:
    array: Cluster labels.
    """
    try:
        dist_dtc = pdist(data, distance)
        link_dtc = linkage(dist_dtc, link)
        clusters = cutreeHybrid(link_dtc, dist_dtc, minClusterSize=min_clust_size_dtc)
        return clusters
    except Exception as e:
        st.error(f"Error during Dynamic Tree Cut clustering: {e}")
        return None


@log_args_clust_gene
@st.cache_data(show_spinner="HDBSCAN Clustering is Calculating", max_entries=10, ttl=1800)
def hdb_gene(data, alfa, hdb_met, mcs_hdb, min_sam_hdb, eps, clust_meth):
    """
    Perform HDBSCAN clustering on gene data.

    Parameters:
    data (DataFrame): Input data for HDBSCAN.
    alfa (float): Alpha parameter for HDBSCAN.
    hdb_met (str): Metric for HDBSCAN.
    mcs_hdb (int): Minimum cluster size for HDBSCAN.
    min_sam_hdb (int): Minimum samples for HDBSCAN.
    eps (float): Epsilon parameter for HDBSCAN.
    clust_meth (str): Cluster selection method for HDBSCAN.

    Returns:
    array: Cluster labels.
    """
    try:
        clusters = HDBSCAN(algorithm='auto', alpha=alfa, leaf_size=40, metric=hdb_met, min_cluster_size=mcs_hdb,
                           min_samples=min_sam_hdb, cluster_selection_method=clust_meth,
                           allow_single_cluster=False, cluster_selection_epsilon=eps).fit_predict(data)
        return clusters
    except Exception as e:
        st.error(f"Error during HDBSCAN clustering: {e}")
        return None


@log_args_clust_gene
@st.cache_data(show_spinner="OPTICS Clustering is Calculating", max_entries=10, ttl=1800)
def optics_gene(data, ms_optics, metric_opt, eps_opt, alg_opt, cluster_met):
    """
    Perform OPTICS clustering on gene data.

    Parameters:
    data (DataFrame): Input data for OPTICS.
    ms_optics (int): Minimum samples for OPTICS.
    metric_opt (str): Metric for OPTICS.
    eps_opt (float): Epsilon parameter for OPTICS.
    alg_opt (str): Algorithm for OPTICS.
    cluster_met (str): Cluster selection method for OPTICS.

    Returns:
    array: Cluster labels.
    """
    try:
        clustering = OPTICS(min_samples=ms_optics, cluster_method=cluster_met, metric=metric_opt,
                            eps=eps_opt, algorithm=alg_opt).fit_predict(data)
        return clustering
    except Exception as e:
        st.error(f"Error during OPTICS clustering: {e}")
        return None


@log_args_dimred_gene
@st.cache_data(show_spinner="PHATE Calculating...", max_entries=10, ttl=1800)
def phate_dimension_reduction_gene(data, lat_vec, n_neight, random_sample, n_pca):
    """
    Perform PHATE dimension reduction on gene data.

    Parameters:
    data (DataFrame): Input data for PHATE.
    lat_vec (int): Number of latent vectors (components) for PHATE.
    n_neight (int): Number of neighbors for PHATE.
    random_sample (int): Random state for reproducibility.
    n_pca (int): Number of principal components for PHATE.

    Returns:
    DataFrame: Reduced data after PHATE.
    """
    try:
        phate_run = phate.PHATE(n_components=lat_vec, knn=n_neight, n_pca=n_pca, random_state=random_sample)
        reduced_data = phate_run.fit_transform(data)
        return reduced_data
    except Exception as e:
        st.error(f"Error during PHATE dimension reduction: {e}")
        return None


@log_args_dimred_gene
@st.cache_data(show_spinner="Correlation Calculating...", max_entries=10, ttl=1800)
def corr_gene(data, method):
    """
    Calculate correlation matrix for gene data.

    Parameters:
    data (DataFrame): Input data for correlation calculation.
    method (str): Method for correlation calculation ('Spearman' or 'Pearson').

    Returns:
    DataFrame: Correlation matrix.
    """
    try:
        if method == "Spearman":
            reduced_data = data.corr(method='spearman')
        else:
            reduced_data = data.corr(method='pearson')
        return reduced_data
    except Exception as e:
        st.error(f"Error during correlation calculation: {e}")
        return None


@log_args_dimred
@st.cache_data(show_spinner="Correlation Calculating...", max_entries=10, ttl=1800)
def corr_sample(data, method):
    """
    Calculate correlation matrix for sample data.

    Parameters:
    data (DataFrame): Input data for correlation calculation.
    method (str): Method for correlation calculation ('Spearman' or 'Pearson').

    Returns:
    DataFrame: Correlation matrix.
    """
    try:
        if method == "Spearman":
            reduced_data = data.corr(method='spearman')
        else:
            reduced_data = data.corr(method='pearson')
        return reduced_data
    except Exception as e:
        st.error(f"Error during correlation calculation: {e}")
        return None
