import streamlit as st
import pandas as pd
from Scripts.plotly_plots import con_hm
from Scripts.clust_dim_def import (
    PCA_dimension_reduction, t_SNE_dimension_reduction, UMAP_dimension_reduction,
    MDS_dimension_reduction, NMF_dimension_reduction, nmf_clust, nmf_data,
    PCA_dimension_reduction_gene, t_SNE_dimension_reduction_gene, UMAP_dimension_reduction_gene,
    MDS_dimension_reduction_gene, NMF_dimension_reduction_gene, phate_dimension_reduction_gene,
    phate_dimension_reduction, corr_gene
)
from Scripts.clustering import save_disp
import Scripts.descriptions as dscript


def dim_reductor(selected_genes1, selected_genes_nmf, random_state_samples, df_names, selected_genes_dm_name, user_id,
                 nmf_name, cluster_name, option_dataset, adv_col1, adv_col2, ht_expression_path,
                 ht_top_annotation_path, ht_png, ht_pdf, surv_png, surv_pdf, ht_row_cluster_path,
                 selected_subtypes, cancer_opt, nmf_cons_name, e, sidebar_expander, gene_clust_chb, clinic_df=None,
                 surv_dataframe_path=None):
    """
    Perform dimensionality reduction on selected genes.

    Parameters:
    - selected_genes1 (DataFrame): Input gene expression data.
    - selected_genes_nmf (DataFrame): NMF data.
    - random_state_samples (int): Random state for reproducibility.
    - df_names (list): List of sample names.
    - selected_genes_dm_name (str): Filename to save dimensionally reduced data.
    - user_id (str): User identifier.
    - nmf_name (str): Filename of NMF data.
    - cluster_name (str): Filename to save clustering results.
    - option_dataset (str): Selected dataset option.
    - adv_col1 (streamlit.Column): Streamlit column for advanced settings.
    - adv_col2 (streamlit.Column): Streamlit column for advanced settings.
    - ht_expression_path (str): Path to save heatmap expression data.
    - ht_top_annotation_path (str): Path to save top annotation data.
    - ht_png (str): Path to save heatmap PNG.
    - ht_pdf (str): Path to save heatmap PDF.
    - surv_png (str): Path to save survival analysis PNG.
    - surv_pdf (str): Path to save survival analysis PDF.
    - ht_row_cluster_path (str): Path to save heatmap row cluster data.
    - selected_subtypes (str): Selected subtypes.
    - cancer_opt (str): Selected cancer option.
    - nmf_cons_name (str): Filename to save NMF consensus matrix.
    - clinic_df (DataFrame, optional): Clinical data.
    - surv_dataframe_path (str, optional): Path to save survival analysis data.

    Returns:
    - selected_genes_dim_red (DataFrame): Dimensionally reduced data.
    - option_dm (str): Selected dimensionality reduction method.
    - lat_vec (int): Number of latent vectors.
    - cluster_db (DataFrame): Clustering database.
    """
    selected_genes_dim_red_scatter = None
    selected_genes_dim_red = None
    cluster_db = pd.DataFrame()

    if int(selected_genes1.shape[0]) <= 60:
        multiplier_dimred = 3
    else:
        multiplier_dimred = int(selected_genes1.shape[0] * 0.05)

    if len(selected_genes1.columns) < 4:
        less_number = True
        dm_selectbox_option = ["--no dimensionality reduction--", "UMAP", "t-SNE", "MDS", "NMF", "PHATE"]

    else:
        dm_selectbox_option = ["--no dimensionality reduction--", "UMAP", "t-SNE", "PCA", "MDS", "NMF", 'NMF Clustering', "PHATE", "Correlation"]
        less_number = False

    three_d_empty = sidebar_expander.empty()
    lat_vec_chb = three_d_empty.toggle(':blue[Analyse in ] :red[3D]', key="3d_chb", value=True)
    if lat_vec_chb:
        lat_vec = 3
    else:
        lat_vec = 2

    option_dm = sidebar_expander.selectbox(
        ':blue[Dimensionality ] :red[Reduction]',
        options=dm_selectbox_option,
        help="Dimension reduction methods are techniques used to reduce the number of "
             "variables or features in a dataset while retaining important information. "
             "They can be categorized into two main types: feature selection, "
             "which involves selecting a subset of the original features, and feature "
             "extraction, which involves creating new features that capture the most "
             "important information from the original features.", key="opt_dm", index=1)

    def save_reduced_data(data, filename):
        df_dr = pd.DataFrame(data)
        df_dr["Samples"] = df_names
        df_dr.to_csv(filename, index=False)

    try:
        if option_dm == "UMAP":
            adv_col2.subheader(":blue[Advanced Settings of UMAP]", help=dscript.umap_adv())
            min_dis = adv_col2.number_input(":red[Minimal Distance] :blue[Between Points]", min_value=0.00, max_value=1.00,
                                            value=0.00, step=0.01, key="umap_md")
            n_neigh = adv_col2.number_input(":blue[Number of ] :red[Neighbours]", min_value=2,
                                            max_value=int(selected_genes1.shape[0] / 2),
                                            value=multiplier_dimred, key="umap_neigh")
            umap_metric = adv_col2.selectbox(":blue[Select ] :red[Metric]", key="umap_metr",
                                             options=['euclidean', 'manhattan', 'chebyshev', 'minkowski', 'canberra',
                                                      'braycurtis', 'mahalanobis', 'cosine', 'correlation'], index=0)
            selected_genes_dim_red = UMAP_dimension_reduction(selected_genes1, lat_vec, min_dis, n_neigh, random_state_samples, umap_metric)
            save_reduced_data(selected_genes_dim_red, selected_genes_dm_name)

        elif option_dm == "t-SNE":
            if selected_genes1.shape[0] <= 10:
                perp_default = 5
            else:
                perp_default = int(selected_genes1.shape[0] / 2)
            adv_col2.subheader(":blue[Advanced Settings of t-SNE]", help=dscript.tsne_adv())
            perp = adv_col2.number_input(":red[Perplexity of ] :blue[t-SNE]", min_value=5, max_value=selected_genes1.shape[0], step=1, value=perp_default, key="perp")
            tsne_lr_chb = adv_col2.toggle(":blue[Automatic ] :red[Learning Rate]", key="tsne_lr_auto", value=True)
            lr = 'auto' if tsne_lr_chb else adv_col2.number_input(":red[Learning ] :blue[Rate]", min_value=10, max_value=1000, value=50, key="lr_ts", step=1)
            niter = adv_col2.number_input(":blue[Number of ] :red[Iterations] :blue[(More Iterations Slow Down the Run)]", min_value=250, max_value=2000, value=1000, key="niter_ts", step=1)
            met = adv_col2.selectbox(":red[Distance] :blue[Metric]", options=("euclidean", "cosine", "cityblock", "manhattan"), key="tsne_sb", index=0)
            if less_number:
                init_meth = "random"
            else:
                init_meth = ("random", "pca")
            init = adv_col2.selectbox(":blue[Initialisation of ] :red[Embedding Method]", options=init_meth, key="tsne_init", index=0)
            selected_genes_dim_red = t_SNE_dimension_reduction(selected_genes1, lat_vec, perp, lr, niter, met, random_state_samples, init)
            save_reduced_data(selected_genes_dim_red, selected_genes_dm_name)

        elif option_dm == "PCA":
            adv_col2.subheader(":blue[Advanced Settings of PCA]", help=dscript.pca_adv())
            tolerance = adv_col2.number_input(":blue[Relative ] :red[Tolerance]", min_value=0.0, max_value=1.0, value=0.1, key="pca_tol", step=0.1, format="%f")
            whiten = adv_col2.selectbox(":blue[Use ] :red[Whiten ]", key="white", options=(True, False), index=1)
            svd_solver = adv_col2.selectbox(":blue[Select ] :red[SVD Solver]", key="svd_solver", options=('auto', 'full', 'covariance_eigh', 'arpack', 'randomized'), index=0)
            selected_genes_dim_red = PCA_dimension_reduction(selected_genes1, lat_vec, tolerance, whiten=whiten, rs_pca=random_state_samples, svd_solver_app=svd_solver)
            save_reduced_data(selected_genes_dim_red, selected_genes_dm_name)

        elif option_dm == "NMF":
            if less_number:
                init_alg_option = "random"
            else:
                init_alg_option = ("nndsvda", "nndsvd", "nndsvdar", "random")
            adv_col2.subheader(":blue[Advanced Settings of NMF]", help=dscript.nmf_adv())
            max_it = adv_col2.number_input(":red[Maximum Number] :blue[ of Iterations]", min_value=200, max_value=1000, value=250, step=1, key="mit")
            toler = adv_col2.select_slider(":red[Tolerance ] :blue[of Stopping Condition]", options=(0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5), value=0.00001, key="nmf_tol")
            init_alg = adv_col2.selectbox(":blue[Method of  ] :red[Initialize]", options=init_alg_option, key="i_alg_df", index=0)
            b_loss = adv_col2.selectbox(":blue[Beta ] :red[Loss]", options=("kullback-leibler", "frobenius"), key="b_loss", index=1)
            solver = adv_col2.selectbox(":blue[Select ] :red[Numerical Solver]", key="nmf_solver", options=("cd", "mu"), index=0) if b_loss == "frobenius" else "mu"
            shuffle = adv_col2.selectbox(":blue[Randomize the ] :red[Order of Coordinates ]", key="nmf_shuffle", options=(True, False), index=1) if solver == "cd" else False
            selected_genes_dim_red = NMF_dimension_reduction(selected_genes_nmf, lat_vec, max_it, toler, init_alg, b_loss, random_state_samples, solver, shuffle)
            save_reduced_data(selected_genes_dim_red, selected_genes_dm_name)

        elif option_dm == "NMF Clustering":
            three_d_empty.empty()
            n_clust = st.sidebar.select_slider(
                ':blue[Select the ] :red[Number of Clusters]',
                options=[2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15], key="nmf_num_clust")
            st.sidebar.write(f':blue[Number of cluster] :red[{str(n_clust)}]')

            adv_col2.subheader(":blue[Advanced Settings of NMF Clustering]", help=dscript.nmf_clust_adv())
            iter_conmat = adv_col2.number_input(":blue[Number of ] :red[Iterations] :blue[(More iterations slow down the run)]",
                                                min_value=5, max_value=500, value=100, key="iter_con")
            tri_nmf = adv_col2.number_input(":blue[Number of ] :red[Trials]", key="tri_nmf",
                                            min_value=1, max_value=100, value=50)
            lam_nmf = adv_col2.number_input(":red[Lamb's ] :blue[Value]", key="nmf_lamb",
                                            min_value=0.01, max_value=1.0, step=0.01, value=0.01)

            data_dict, samples_names1, samples_names2, name_order = nmf_data(nmf_name)
            model, clust_nmf = nmf_clust(data_dict, n_clust, iter_conmat, tri_nmf, lam_nmf)
            clust_nmf = pd.DataFrame(clust_nmf['sim1'], dtype="int")
            clust_nmf = clust_nmf.apply(lambda x: x.idxmax(), axis=0)
            clust_nmf = clust_nmf.str.replace('class-', '')
            clust_nmf = pd.DataFrame({'Cluster': clust_nmf.astype(int) + 1}, index=clust_nmf.index)
            clust_nmf['Cluster'] = clust_nmf['Cluster'].astype("str")
            clust_nmf['Cluster'] = clust_nmf['Cluster'].apply(lambda x: "S" + str(x))
            clust_nmf = clust_nmf.reindex(index=name_order)
            num_of_cluster = len(clust_nmf["Cluster"].unique())
            clust_nmf.index.names = ['Samples']
            clust_nmf.to_csv(cluster_name, index=True)

            adv_col1.subheader(":blue[Color Settings of Consensus Matrix Heatmap]")
            top_color = adv_col1.color_picker(":red[Max ] :blue[Value's Color]", value="#FF0000", key="top_color")
            middle_color = adv_col1.color_picker(":red[Mid ] :blue[Value's Color]", value="#FFFFFF", key="middle_color")
            bot_color = adv_col1.color_picker(":red[Min ] :blue[Value's Color]", value="#3F3F87", key="bot_color")
            cmap = [bot_color, middle_color, top_color]

            con_hm(model, cmap, samples_names1, nmf_cons_name, e, gene_clust_chb)

            selected_genes_dim_red = pd.DataFrame()
            labels = pd.DataFrame()
            labels = clust_nmf.copy()
            labels = labels.reset_index()
            labels["0"] = labels["Cluster"]
            labels = labels.drop(["Samples", "Cluster"], axis=1)
            lat_vec = 2
            cluster_db = save_disp(df_names, cluster_name, lat_vec, selected_genes_dim_red,
                                   labels, option_dataset, clinic_df, surv_dataframe_path, user_id,
                                   ht_expression_path, cluster_name, ht_top_annotation_path,
                                   ht_png, ht_pdf, surv_png, surv_pdf, ht_row_cluster_path,
                                   selected_subtypes, None, None, adv_col2,
                                   selected_genes_dm_name, option_dm, e, cancer_opt)
        elif option_dm == "MDS":
            adv_col2.subheader(":blue[Advanced Settings of MDS]", help=dscript.mds_adv())
            mi = adv_col2.number_input(":red[Maximum Number] :blue[ of Iterations]", min_value=200, max_value=1000, value=1000, key="mds_mi")
            eps_mds = adv_col2.select_slider(":red[Epsilon ] :blue[Value]", options=(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05), value=0.0001, key="eps")
            metric_mds = adv_col2.selectbox(":blue[Use ] :red[Metric MDS]", key="mtrc_mds", options=[True, False], index=0)
            selected_genes_dim_red = MDS_dimension_reduction(selected_genes1, lat_vec, mi, eps_mds, random_state_samples, metric_mds)
            save_reduced_data(selected_genes_dim_red, selected_genes_dm_name)

        elif option_dm == "PHATE":
            adv_col2.subheader(":blue[Advanced Settings of PHATE]", help=dscript.phate_adv())
            n_neight = adv_col2.number_input(":blue[Number of ] :red[Neighbours]", min_value=2, max_value=int(selected_genes1.shape[0] / 2), value=multiplier_dimred, key="phate_neigsh")
            n_pca = adv_col2.number_input(":blue[Number of ] :red[Principal Components]", min_value=2, max_value=int(selected_genes1.shape[0] / 2), value=multiplier_dimred, key="phate_npca")
            selected_genes_dim_red = phate_dimension_reduction(selected_genes1, lat_vec, n_neight, random_state_samples, n_pca)
            save_reduced_data(selected_genes_dim_red, selected_genes_dm_name)

        elif option_dm == "Correlation":
            three_d_empty.empty()
            adv_col2.subheader(":blue[Advanced Settings of Correlation]", help=dscript.corr_adv())
            method = adv_col2.selectbox(":blue[Correlation ] :red[Method]", options=["Spearman", "Pearson"], index=0, key="corr_sample")
            selected_genes_dim_red = corr_gene(selected_genes1.T, method)
            save_reduced_data(selected_genes_dim_red, selected_genes_dm_name)

        else:
            if selected_genes1.shape[1] <= 3:
                selected_genes_dim_red = UMAP_dimension_reduction(selected_genes1, lat_vec, min_dis=0.1, n_neigh=3, rs_umap=42, metric="euclidean")
            else:
                selected_genes_dim_red = PCA_dimension_reduction(selected_genes1, lat_vec, tolerance=10e-3, rs_pca=42, whiten=True)
            save_reduced_data(selected_genes1, selected_genes_dm_name)

    except Exception as e:
        st.error(f"Error during dimensionality reduction: {e}")
        return None, option_dm, lat_vec, cluster_db

    return selected_genes_dim_red, option_dm, lat_vec, cluster_db


def dim_reductor_gene(selected_genes1, selected_genes_nmf, random_state_samples, df_names, selected_genes_dm_name,
                      user_id, nmf_name, cluster_name, option_dataset, adv_col1_gene, adv_col2_gene,
                      ht_png, ht_pdf, surv_png, surv_pdf, gene_clust_col2, gene_exp,
                      clinic_df=None, surv_dataframe_path=None):
    """
    Perform dimensionality reduction on selected genes (gene-specific).

    Parameters:
    - selected_genes1 (DataFrame): Input gene expression data.
    - selected_genes_nmf (DataFrame): NMF data.
    - random_state_samples (int): Random state for reproducibility.
    - df_names (list): List of sample names.
    - selected_genes_dm_name (str): Filename to save dimensionally reduced data.
    - user_id (str): User identifier.
    - nmf_name (str): Filename of NMF data.
    - cluster_name (str): Filename to save clustering results.
    - option_dataset (str): Selected dataset option.
    - adv_col1_gene (streamlit.Column): Streamlit column for advanced settings.
    - adv_col2_gene (streamlit.Column): Streamlit column for advanced settings.
    - ht_png (str): Path to save heatmap PNG.
    - ht_pdf (str): Path to save heatmap PDF.
    - surv_png (str): Path to save survival analysis PNG.
    - surv_pdf (str): Path to save survival analysis PDF.
    - gene_clust_col2 (streamlit.Column): Streamlit column for gene clustering settings.
    - gene_exp (streamlit.Column): Streamlit column for gene expression settings.
    - clinic_df (DataFrame, optional): Clinical data.
    - surv_dataframe_path (str, optional): Path to save survival analysis data.

    Returns:
    - selected_genes_dim_red (DataFrame): Dimensionally reduced data.
    - option_dm (str): Selected dimensionality reduction method.
    - lat_vec (int): Number of latent vectors.
    """
    selected_genes_dim_red_scatter = None
    selected_genes_dim_red = None

    if selected_genes1.shape[0] < 60:
        multiplier_dimred = 3
    else:
        multiplier_dimred = int(selected_genes1.shape[0] * 0.05)

    three_d_empty = gene_exp.empty()
    lat_vec = 3 if three_d_empty.toggle(':blue[Analyse in ] :red[3D]', key="3d_chb_gene", value=True) else 2
    dm_selectbox_option = ["--no dimensionality reduction--", "UMAP", "t-SNE", "PCA", "MDS", "NMF", "PHATE", "Correlation"]

    def save_reduced_data(data, filename):
        df_dr = pd.DataFrame(data)
        df_dr["Samples"] = df_names
        df_dr.to_csv(filename, index=False)

    option_dm = gene_exp.selectbox(
        ':blue[Dimensionality ] :red[Reduction]',
        options=dm_selectbox_option, key="opt_dm_genes", index=1)

    try:
        if option_dm == "UMAP":
            adv_col2_gene.subheader(":blue[Advanced Settings of UMAP (Genes)]", help=dscript.umap_adv())
            min_dis = adv_col2_gene.number_input(":red[Minimal Distance] :blue[Between Points]", min_value=0.00, max_value=1.00, value=0.00, step=0.01, key="gene_umap")
            n_neigh = adv_col2_gene.number_input(":blue[Number of ] :red[Neighbours]", min_value=2, max_value=int(selected_genes1.shape[0] / 2), value=multiplier_dimred, key="umap_gene_neigh")
            umap_metric = adv_col2_gene.selectbox(":blue[Select ] :red[Metric]", key="umap_metr_gene", options=['euclidean', 'manhattan', 'chebyshev', 'minkowski', 'canberra', 'braycurtis', 'mahalanobis', 'cosine', 'correlation'], index=0)
            selected_genes_dim_red = UMAP_dimension_reduction_gene(selected_genes1, lat_vec, min_dis, n_neigh, random_state_samples, umap_metric)
            save_reduced_data(selected_genes_dim_red, selected_genes_dm_name)

        elif option_dm == "t-SNE":
            if selected_genes1.shape[0] <= 10:
                perp_default = 5
            else:
                perp_default = int(selected_genes1.shape[0] / 2)
            adv_col2_gene.subheader(":blue[Advanced Settings of t-SNE (Genes)]", help=dscript.tsne_adv())
            perp = adv_col2_gene.number_input(":red[Perplexity of ] :blue[t-SNE]", min_value=5,
                                              max_value=selected_genes1.shape[0], value=perp_default, step=1, key="perp_gene")
            tsne_lr_chb = adv_col2_gene.toggle(":blue[Automatic ] :red[Learning Rate]", key="tsne_lr_auto_gene", value=True)
            lr = 'auto' if tsne_lr_chb else adv_col2_gene.number_input(":red[Learning ] :blue[Rate]", min_value=10, max_value=1000, value=50, key="lr_ts_gene")
            niter = adv_col2_gene.number_input(":blue[Number of ] :red[Iterations] :blue[(More Iterations Slow Down the Run)]", min_value=250, max_value=2000, value=1000, key="niter_gene", step=1)
            met = adv_col2_gene.selectbox(":red[Distance] :blue[Metric]", options=("euclidean", "cosine", "cityblock", "manhattan"), key="tsne_sb_genes", index=0)
            init = adv_col2_gene.selectbox(":blue[Initialisation of ] :red[Embedding Method]", options=("random", "pca"), key="tsne_init_gene", index=1)
            selected_genes_dim_red = t_SNE_dimension_reduction_gene(selected_genes1, lat_vec, perp, lr, niter, met, random_state_samples, init)
            save_reduced_data(selected_genes_dim_red, selected_genes_dm_name)

        elif option_dm == "PCA":
            adv_col2_gene.subheader(":blue[Advanced Settings of PCA]", help=dscript.pca_adv())
            tolerance = adv_col2_gene.number_input(":blue[Relative ] :red[Tolerance]", min_value=0.0, max_value=1.0, value=0.1, key="pca_tol_gene", step=0.1, format="%f")
            whiten = adv_col2_gene.selectbox(":blue[Use ] :red[Whiten ]", key="white_gene", options=(True, False), index=1)
            selected_genes_dim_red = PCA_dimension_reduction_gene(selected_genes1, lat_vec, tolerance, whiten=whiten, rs_pca=random_state_samples)
            save_reduced_data(selected_genes_dim_red, selected_genes_dm_name)

        elif option_dm == "MDS":
            adv_col2_gene.subheader(":blue[Advanced Settings of MDS (Genes)]", help=dscript.mds_adv())
            mi = adv_col2_gene.number_input(":red[Maximum Number] :blue[ of Iterations]", min_value=200, max_value=1000, value=1000, step=1, key="mi_gene")
            eps_mds = adv_col2_gene.select_slider(":red[Epsilon ] :blue[Value]", options=(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05), value=0.0001, key="eps_gene")
            metric_mds = adv_col2_gene.selectbox(":blue[Use ] :red[Metric MDS]", key="mtrc_mds_gene", options=[True, False], index=0)
            selected_genes_dim_red = MDS_dimension_reduction_gene(selected_genes1, lat_vec, mi, eps_mds, random_state_samples, metric_mds)
            save_reduced_data(selected_genes_dim_red, selected_genes_dm_name)

        elif option_dm == "NMF":
            adv_col2_gene.subheader(":blue[Advanced Settings of NMF (Genes)]", help=dscript.nmf_adv())
            max_it = adv_col2_gene.number_input(":red[Maximum Number] :blue[ of Iterations]", min_value=200, max_value=1000, value=250, step=1, key="mit_gene")
            toler = adv_col2_gene.select_slider(":blue[Relative ] :red[Tolerance]", options=(0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5), value=0.00001, key="nmf_tol_gene")
            init_alg = adv_col2_gene.selectbox(":blue[Choose ] :red[Metric]", options=("nndsvda", "nndsvd", "nndsvdar", "random"), key="i_alg_df_gene", index=0)
            b_loss = adv_col2_gene.selectbox(":blue[Choose ] :red[Metric]", options=("kullback-leibler", "frobenius"), key="b_loss_gene", index=0)
            solver = adv_col2_gene.selectbox(":blue[Select ] :red[Numerical Solver]", key="nmf_solver_gene", options=("cd", "mu"), index=0) if b_loss == "frobenius" else "mu"
            shuffle = adv_col2_gene.selectbox(":blue[Randomize the ] :red[Order of Coordinates ]", key="nmf_shuffle_gene", options=(True, False), index=1) if solver == "cd" else False
            selected_genes_dim_red = NMF_dimension_reduction_gene(selected_genes_nmf, lat_vec, max_it, toler, init_alg, b_loss, random_state_samples, solver, shuffle)
            save_reduced_data(selected_genes_dim_red, selected_genes_dm_name)

        elif option_dm == "PHATE":
            adv_col2_gene.subheader(":blue[Advanced Settings of PHATE]", help=dscript.phate_adv())
            n_neight = adv_col2_gene.number_input(":blue[Number of ] :red[Neighbours]", min_value=2, max_value=int(selected_genes1.shape[0] / 2), value=multiplier_dimred, key="phate_neigh_gene")
            n_pca = adv_col2_gene.number_input(":blue[Number of ] :red[Principal Components]", min_value=2, max_value=int(selected_genes1.shape[0] / 2), value=multiplier_dimred, key="phate_npca_gene")
            selected_genes_dim_red = phate_dimension_reduction_gene(selected_genes1, lat_vec, n_neight, random_state_samples, n_pca)
            save_reduced_data(selected_genes_dim_red, selected_genes_dm_name)

        elif option_dm == "Correlation":
            three_d_empty.empty()
            adv_col2_gene.subheader(":blue[Advanced Settings of Correlation]", help=dscript.corr_adv())
            method = adv_col2_gene.selectbox(":blue[Correlation ] :red[Method]", options=["Spearman", "Pearson"], index=0, key="corr_gene")
            selected_genes_dim_red = corr_gene(selected_genes1.T, method)
            save_reduced_data(selected_genes_dim_red, selected_genes_dm_name)

        else:
            if selected_genes1.shape[1] <= 3:
                selected_genes_dim_red = UMAP_dimension_reduction_gene(selected_genes1, lat_vec, min_dis=0.1, n_neigh=3, rs_umap=42, metric='euclidean')
            else:
                selected_genes_dim_red = PCA_dimension_reduction_gene(selected_genes1, lat_vec, tolerance=10e-3, rs_pca=42, whiten=True)
            save_reduced_data(selected_genes1, selected_genes_dm_name)

    except Exception as e:
        st.error(f"Error during dimensionality reduction: {e}")
        return None, option_dm, lat_vec

    return selected_genes_dim_red, option_dm, lat_vec
