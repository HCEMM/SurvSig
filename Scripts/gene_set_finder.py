import streamlit as st
import pandas as pd
import numpy as np
import math
from sklearn.preprocessing import StandardScaler, KBinsDiscretizer
from sklearn.decomposition import PCA, FastICA, FactorAnalysis, TruncatedSVD, NMF
import Scripts.data_loader as dl
from sklearn.naive_bayes import CategoricalNB
import Scripts.gene_list as gl
import Scripts.descriptions as dsc
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests
from sklearn.neighbors import NearestCentroid
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import classification_report
from Scripts.subtype_selector import subtype_select_gf, subtypes_selector_tcga_gf
from scipy.spatial.distance import cdist


def balanc_dataset(
    expression: pd.DataFrame,
    merged_df: pd.DataFrame,
    top_genes: int = 5000,
    use_pca: bool = True,
    n_components: int = 2,
    representative_method: str = "Centroid",  # or "Medoid"
    distance_metric: str = "Euclidean",       # or "Pearson"
    random_state: int = 42,
    varience_reduction: bool = True,
    gene_reduction: int = 5000
) -> (pd.DataFrame, pd.DataFrame):


    """
    Downsamples each group to the size of the smallest group, returning:
      1) final_expr  : Balanced expression DataFrame
      2) final_merged: Balanced merged_df (annotation table)

    Steps:
      A) Derive group labels from merged_df (assume first column).
      B) Determine min_count = smallest group size.
      C) (Optional) select top variable genes in expression.
      D) (Optional) PCA transform the entire dataset.
      E) For each group:
         - If group size <= min_count, keep them all.
         - Else pick min_count nearest to centroid/medoid.
      F) Subset expression & merged_df to those chosen sample IDs.
    """

    # A) Extract group labels from merged_df
    #    We'll assume the group column is the first column in merged_df:
    group_series = merged_df.iloc[:, 0]  # shape (n_samples,)
    group_series.index = merged_df.index  # ensure same index

    # B) Determine the smallest group size
    group_sizes = group_series.value_counts()
    min_count = group_sizes.min()  # e.g., 19 if NE=62 and non-NE=19

    # Copy to avoid mutating original data
    df_expr = expression.copy()
    df_merge = merged_df.copy()

    # C) Select top variable genes
    if top_genes < df_expr.shape[1]:
        variances = df_expr.var(axis=0)
        top_gene_names = variances.nlargest(top_genes).index
        df_expr = df_expr[top_gene_names]

    # D) (Optional) PCA on entire dataset
    data_matrix = df_expr.values  # shape (n_samples, top_genes)
    if use_pca:
        pca_model = PCA(n_components=n_components, random_state=random_state)
        data_matrix = pca_model.fit_transform(data_matrix)  # shape (n_samples, n_components)
    else:
        data_matrix = df_expr.values

    reduced_df = pd.DataFrame(
        data_matrix,
        index=df_expr.index,
        columns=[f"PC{i+1}" for i in range(data_matrix.shape[1])]
    )

    # E) Per-group downsampling
    final_indices = []
    np.random.seed(random_state)  # for reproducibility if tie-breaking

    for group_label, size in group_sizes.items():
        # All sample IDs for this group
        group_samples = group_series[group_series == group_label].index

        if size <= min_count:
            # Already at or below min_count, keep them all
            final_indices.extend(group_samples)
        else:
            # group size > min_count, pick the top min_count by distance
            group_data = reduced_df.loc[group_samples].values  # shape: (#group_size, n_dims)

            # 1) Compute representative
            if representative_method == "Centroid":
                rep_vec = group_data.mean(axis=0, keepdims=True)  # shape (1, n_dims)
            else:
                # Medoid => sample that minimizes sum of Euclidean distances
                pdist = cdist(group_data, group_data, metric='euclidean')
                sum_dists = pdist.sum(axis=1)
                medoid_idx = np.argmin(sum_dists)
                rep_vec = group_data[medoid_idx:medoid_idx+1, :]

            # 2) Distances
            if distance_metric == "Euclidean":
                dists = cdist(group_data, rep_vec, metric='euclidean').ravel()
            else:
                # 'Pearson' => 1 - correlation
                rep_1d = rep_vec.ravel()
                corr_vals = []
                for i in range(group_data.shape[0]):
                    c = np.corrcoef(group_data[i, :], rep_1d)[0, 1]
                    corr_vals.append(c)
                dists = 1 - np.array(corr_vals)

            # 3) Sort ascending, pick top min_count
            sorted_idx = np.argsort(dists)
            chosen_positions = sorted_idx[:min_count]
            chosen_ids = group_samples[chosen_positions]

            final_indices.extend(chosen_ids)

    # F) Build final subsets of expression & merged
    if varience_reduction:
        variances_2 = expression.var(axis=0)
        top_gene_names_2 = variances_2.nlargest(gene_reduction).index
        expression = expression[top_gene_names_2]

    final_expr = expression.loc[final_indices].copy()
    final_merged = df_merge.loc[final_indices].copy()

    return final_expr, final_merged


def catBN(n_genes, df_clinical, df_expression, random_state):
    df_clinical = df_clinical.dropna().copy()
    common_samples = df_clinical.index.intersection(df_expression.index)
    df_expression = df_expression.loc[common_samples].copy()
    labels = df_clinical.loc[common_samples].copy()
    labels = labels.iloc[:, 0].astype('category')

    df_expression = df_expression.T.copy()

    assert df_expression.shape[1] == labels.shape[0], "Mismatch in number of samples between df_expression and labels"

    X_train, X_test, y_train, y_test = train_test_split(df_expression.T, labels, test_size=0.2,
                                                        random_state=random_state, stratify=labels)

    # Scaling and discretization
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    discretizer = KBinsDiscretizer(n_bins=10, encode='ordinal', strategy='uniform')
    X_train_discretized = discretizer.fit_transform(X_train_scaled)
    X_test_discretized = discretizer.transform(X_test_scaled)

    # Feature selection
    feature_selector = SelectKBest(score_func=f_classif, k=n_genes)
    X_train_selected = feature_selector.fit_transform(X_train_discretized, y_train)
    X_test_selected = feature_selector.transform(X_test_discretized)

    selected_genes = X_train.columns[feature_selector.get_support()]

    # CategoricalNB classification
    classifier = CategoricalNB()
    classifier.fit(X_train_selected, y_train)

    y_pred = classifier.predict(X_test_selected)

    return selected_genes


@st.cache_data(show_spinner="MLP Classifier Calculating...", max_entries=10, ttl=1800)
def mlp_scipy(n_genes, df_clinical, df_expression, random_state):

    df_clinical = df_clinical.dropna().copy()
    common_samples = df_clinical.index.intersection(df_expression.index)
    df_expression = df_expression.loc[common_samples].copy()
    labels = df_clinical.loc[common_samples].copy()
    labels = labels.iloc[:, 0].astype('category')

    # Convert labels to numeric encoding
    from sklearn.preprocessing import LabelEncoder
    label_encoder = LabelEncoder()
    labels = label_encoder.fit_transform(labels)

    df_expression = df_expression.T.copy()

    assert df_expression.shape[1] == len(labels), "Mismatch in number of samples between df_expression and labels"

    X_train, X_test, y_train, y_test = train_test_split(df_expression.T, labels, test_size=0.2,
                                                        random_state=random_state, stratify=labels)

    # Feature selection and scaling
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    feature_selector = SelectKBest(score_func=f_classif, k=n_genes)
    X_train_selected = feature_selector.fit_transform(X_train_scaled, y_train)
    X_test_selected = feature_selector.transform(X_test_scaled)

    selected_genes = X_train.columns[feature_selector.get_support()]

    # MLP classification
    classifier = MLPClassifier(random_state=random_state, solver='adam', max_iter=1000, early_stopping=True,
                               hidden_layer_sizes=(100, 100), activation='relu')

    classifier.fit(X_train_selected, y_train)

    y_pred = classifier.predict(X_test_selected)

    # Debug outputs
    print(f"y_pred: {y_pred[:5]}")
    print(f"y_pred dtype: {np.array(y_pred).dtype}")

    print(classification_report(y_test, y_pred))

    return selected_genes



@st.cache_data(show_spinner="Nearest Centroids Calculating...", max_entries=10, ttl=1800)
def clanc(n_genes, df_clinical, df_expression, random_state):
    df_clinical = df_clinical.dropna().copy()
    common_samples = df_clinical.index.intersection(df_expression.index)
    df_expression = df_expression.loc[common_samples].copy()
    labels = df_clinical.loc[common_samples].copy()
    labels = labels.iloc[:, 0].astype('category')

    df_expression = df_expression.T.copy()

    assert df_expression.shape[1] == labels.shape[0], "Mismatch in number of samples between df_expression and labels"

    X_train, X_test, y_train, y_test = train_test_split(df_expression.T, labels, test_size=0.2,
                                                        random_state=random_state, stratify=labels)

    # Feature selection
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # Feature selection
    feature_selector = SelectKBest(score_func=f_classif, k=n_genes)
    X_train_selected = feature_selector.fit_transform(X_train_scaled, y_train)
    X_test_selected = feature_selector.transform(X_test_scaled)

    # Nearest Centroid classification
    classifier = NearestCentroid()
    classifier.fit(X_train_selected, y_train)

    y_pred = classifier.predict(X_test_selected)

    selected_genes = X_train.columns[feature_selector.get_support()]
    return selected_genes


@st.cache_data(show_spinner="Kruskal Wallis Calculating...", max_entries=10, ttl=1800)
def kruskal_wallis_test(expression_df, merged_df, p_correction_kruskal, kruskal_alpha):
    results = []
    unique_groups = merged_df.iloc[:, 0].unique()  # Get unique groups from the first (and only) column of merged_df
    print(f"Unique groups: {unique_groups}")

    for gene in expression_df.columns:
        group_values = []
        for group in unique_groups:
            indices = merged_df[merged_df.iloc[:, 0] == group].index
            group_data = expression_df.loc[indices, gene].dropna()
            group_values.append(group_data)
            print(f"Group {group} (size {len(group_data)}): {group_data.values}")

        if all(len(g) > 1 for g in group_values) and len(group_values) > 1:
            try:
                statistic, p_value = kruskal(*group_values)
                results.append({'Gene': gene, 'Statistic': statistic, 'P-Value': p_value})
            except ValueError as e:
                print(f"Error with gene {gene}: {str(e)}")
                continue

    results_df = pd.DataFrame(results)
    if results_df.empty:
        return pd.DataFrame(), pd.DataFrame()

    p_values = results_df['P-Value'].values
    reject, pvals_corrected, _, _ = multipletests(p_values, alpha=kruskal_alpha, method=p_correction_kruskal)
    results_df['Adjusted P-Value'] = pvals_corrected
    results_df['Significant'] = reject

    significant_results = results_df[results_df['Significant']].sort_values(by='Adjusted P-Value')

    return significant_results, significant_results["Gene"]


@st.cache_data(show_spinner="PCA Calculating...", max_entries=10, ttl=1800)
def PCA_dimension_reduction(data, lat_vec, rs_pca):
    pca = PCA(n_components=lat_vec, random_state=rs_pca, power_iteration_normalizer="auto")
    reduced_data = pca.fit_transform(data)
    loadings = pca.components_.T
    return loadings


@st.cache_data(show_spinner="ICA Calculating...", max_entries=10, ttl=1800)
def ICA_dimension_reduction(data, lat_vec, rs_ica):
    ica = FastICA(n_components=lat_vec, random_state=rs_ica)
    comps = ica.fit_transform(data)
    loadings = ica.components_.T
    return loadings


@st.cache_data(show_spinner="FA Calculating...", max_entries=10, ttl=1800)
def FA_dimension_reduction(data, lat_vec, rs_fa):
    fa = FactorAnalysis(n_components=lat_vec, random_state=rs_fa)
    comps = fa.fit_transform(data)
    loadings = fa.components_.T
    return loadings


@st.cache_data(show_spinner="SVD Calculating...", max_entries=10, ttl=1800)
def tsvd_dimension_reduction(data, lat_vec, rs_tsvd):
    tsvd = TruncatedSVD(n_components=lat_vec, random_state=rs_tsvd)
    comps = tsvd.fit_transform(data)
    loadings = tsvd.components_.T
    return loadings


@st.cache_data(show_spinner="NMF Calculating...", max_entries=10, ttl=1800)
def NMF_dimension_reduction(data, lat_vec, rs_nmf):
    if np.any(data < 0):
        data -= data.min()
    nmf = NMF(n_components=lat_vec, random_state=rs_nmf)
    reduced_data = nmf.fit_transform(data)
    loadings = nmf.components_.T
    return loadings


def gene_finder(data_set_option, df, gene_sig_name, user_id, clinic_st, clinic_df):
    gene_finder_cols = st.columns(3)

    method_sb = gene_finder_cols[0].selectbox(":blue[Select a ] :red[Method of Analysing]", key="method_sb",
                                              options=["Analysing the Dataset", "Clinical Detail Based Analysing"],
                                              index=0, help=dsc.describe_method_sb_ana_opt())

    subtype_exp = st.expander(":blue[Subtype Selection for ] :red[Analysis]")
    if method_sb == "Analysing the Dataset":
        selected_genes = ml_gene_finder(data_set_option, df, gene_sig_name, user_id, subtype_exp, clinic_st, clinic_df,
                                        gene_finder_cols)
    else:
        selected_genes = clinical_gene_finder(data_set_option, df, gene_sig_name, user_id, clinic_st, clinic_df,
                                              subtype_exp, gene_finder_cols)

    return selected_genes


def ml_gene_finder(data_set_option, df, gene_sig_name, user_id, subtype_exp, clinic_st, clinic_df, gene_finder_cols):
    var_chb = False
    selected_genes_1 = None
    var_chb = None
    rank_info = None
    num_comp_selector = None
    loadings = None
    loadings_df = None
    percentile_selection = None
    merged_df = None
    st_annotation = None
    key = "feature_selector2"

    gene_sig_exp = st.expander(":blue[Advanced Settings]")
    selectbox_options = ["George-SCLC", "Anish-SCLC", 'Jiang-SCLC', "George-LCNEC", "Fernandez-Carcinoid",
                         "Alcala-Carcinoid", "Liu-SCLC", 'Rousseaux-Mixed(non-NE)']

    if data_set_option == "TCGA":
        dataset_index = 8
        selectbox_options = ["George-SCLC", "Anish-SCLC", 'Jiang-SCLC', "George-LCNEC",  "Fernandez-Carcinoid",
                             "Alcala-Carcinoid", "Liu-SCLC", 'Rousseaux-Mixed(non-NE)', "TCGA"]

    elif data_set_option == "George-SCLC":
        dataset_index = 0
    elif data_set_option == "Anish-SCLC":
        dataset_index = 1
    elif data_set_option == "Jiang-SCLC":
        dataset_index = 2
    elif data_set_option == "George-LCNEC":
        dataset_index = 3
    elif data_set_option == "Fernandez-Carcinoid":
        dataset_index = 4
    elif data_set_option == "Alcala-Carcinoid":
        dataset_index = 5
    elif data_set_option == "Liu-SCLC":
        dataset_index = 6
    else:
        dataset_index = 7

    dataset_selection = gene_sig_exp.selectbox(
        ":blue[Choose a ] :red[Database ] :blue[to Serve as the Basis for the Analysis]",
        options=selectbox_options, index=dataset_index, key='ds_select2', help=dsc.describe_ds_select2())

    if dataset_selection == "George-SCLC":
        data, _, merged_df = dl.george_sclc_data_load()
        data.index = data.index.str.upper()
        feature_index = 3
        comp_value = 1
        var_value = 95.00
        def_comp_num = 2000
        rank_meth = 1

    elif dataset_selection == "Anish-SCLC":
        data, _, merged_df = dl.anish_sclc_data_load()
        data.index = data.index.str.upper()
        feature_index = 3
        comp_value = 2
        var_value = 97.50
        def_comp_num = 2000
        rank_meth = 1

    elif dataset_selection == "Jiang-SCLC":
        data, _, merged_df = dl.jiang_sclc_data_load()
        data.index = data.index.str.upper()
        feature_index = 1
        comp_value = 3
        var_value = 90.00
        def_comp_num = 1000
        rank_meth = 1

    elif dataset_selection == "George-LCNEC":
        data, _, merged_df = dl.george_lcnec_data_load()
        data.index = data.index.str.upper()
        feature_index = 1
        comp_value = 3
        var_value = 75.00
        def_comp_num = 1000
        rank_meth = 0

    elif dataset_selection == "Fernandez-Carcinoid":
        data, _, merged_df = dl.carcinoid_data_load()
        data.index = data.index.str.upper()
        feature_index = 3
        comp_value = 1
        var_value = 95.00
        def_comp_num = 2000
        rank_meth = 1

    elif dataset_selection == "Liu-SCLC":
        data, _, merged_df = dl.liu_data_load()
        data.index = data.index.str.upper()
        feature_index = 3
        comp_value = 1
        var_value = 95.00
        def_comp_num = 2000
        rank_meth = 1

    elif dataset_selection == "Alcala-Carcinoid":
        data, _, merged_df = dl.alcala_data_load(False)
        data.index.rename("GeneSymbol", inplace=True)
        data.index = data.index.str.upper()
        feature_index = 3
        comp_value = 1
        var_value = 50.00
        def_comp_num = 1500
        rank_meth = 1

    elif dataset_selection == "Rousseaux-Mixed(non-NE)":
        data, _, merged_df = dl.rousseaux_data_load(False, None)
        data.index = data.index.str.upper()
        feature_index = 4
        comp_value = 1
        var_value = 95.00
        def_comp_num = 500
        rank_meth = 1

    elif dataset_selection == "TCGA":
        data = df
        feature_index = 3
        comp_value = 2
        var_value = 95.00
        def_comp_num = 2000
        rank_meth = 1

    else:
        st.info(":blue[Choose ] :red[dataset]")
        data = None
        return None

    if dataset_selection == data_set_option:                                           
        data = df.copy()

    selected_genes_anno = data.copy()
    selected_genes_anno = selected_genes_anno.transpose()
    if dataset_selection == 'Anish-SCLC':
        st_annotation = dl.data_load_anno_anish(selected_genes_anno)
    elif dataset_selection == 'George-SCLC':
        st_annotation = dl.data_load_anno_george(selected_genes_anno)
    elif dataset_selection == 'George-LCNEC':
        st_annotation = dl.data_load_anno_george_lcnec(selected_genes_anno)
    elif dataset_selection == 'Fernandez-Carcinoid':
        st_annotation = dl.data_load_anno_carcinoid(selected_genes_anno)
    elif dataset_selection == 'Liu-SCLC':
        st_annotation = dl.data_load_anno_liu(selected_genes_anno)
    elif dataset_selection == 'Alcala-Carcinoid':
        st_annotation = dl.data_load_anno_alcala(selected_genes_anno)
    elif dataset_selection == 'Rousseaux-Mixed(non-NE)':
        st_annotation = dl.data_load_anno_rousseaux(selected_genes_anno)
    elif dataset_selection == 'Jiang-SCLC':
        st_annotation = dl.data_load_anno_jiang(selected_genes_anno)

    if dataset_selection != "TCGA":
        data, merged_df = subtype_select_gf(data.copy().transpose(), st_annotation.copy(), data.copy().transpose(),
                                            dataset_selection, merged_df.copy(), subtype_exp, key='subtype_ne_ml')
    else:
        data, merged_df = subtypes_selector_tcga_gf(data.copy(), clinic_df.copy(), clinic_st.copy(), subtype_exp,
                                                    key='subtype_tcga')

    gene_list_opt = gene_finder_cols[2].toggle(":blue[Use ] :red[Gene List]", key='gl_toggle2',
                                        help=dsc.describe_gl_toggle2())

    if gene_list_opt:
        gene_finder_cols[2].info(":blue[You can upload your file on the ] :red[sidebar!]")

    if gene_list_opt:
        gene_sig_exp.info(":blue[You can check your gene list in ] :red[Gene Handling Tab]")
        data = gl.app_maker(user_session_id=user_id, option_dataset=dataset_selection, df=data)

    if data is not None:
        # Separating out the features (assuming your data is just the expression values)
        data.dropna(axis=1, how='all', inplace=True)
        data = data.transpose()
        # Store the column names before converting data to numpy array
        column_names = data.columns
        # Converting data to numpy array
        data.dropna(how='any', axis=0, inplace=True)
        data = data.values

        # Standardizing the features
        x = StandardScaler().fit_transform(data)

        feature_selector = gene_sig_exp.selectbox(":blue[Choose the ] :red[Feature Selector Method]",
                                                  options=["PCA", "ICA", "FA", "SVD", "Standard Deviation", "NMF"],
                                                  index=feature_index, key="feat_selector_meth2",
                                                  help=dsc.describe_feat_selector_meth2())

        if feature_selector not in ["Standard Deviation"]:
            num_comp = gene_sig_exp.slider(":blue[Select the ] :red[Number of Components]", min_value=1,
                                           max_value=50, value=comp_value, key="pca_compo_slider2",
                                           help=dsc.describe_pca_compo_slider2())

            random_state = gene_sig_exp.number_input(":blue[Set the ] :red[Random State ] :blue[of Your Work]",
                                                     min_value=1, max_value=2 ** 32 - 1, value=42,
                                                     key="random_state_gene_sign2",
                                                     help=dsc. describe_random_state_gene_sign2())
        else:
            num_comp = None
            random_state = None

        if data.shape[1] < 2500:
            max_genes = int(data.shape[1])
            if def_comp_num > int(data.shape[1]):
                def_comp_num = int(data.shape[1] / 2)
        else:
            max_genes = 2500

        if gene_list_opt:
            max_genes = data.shape[1]

        if feature_selector == "PCA":
            loadings = PCA_dimension_reduction(x, num_comp, random_state)

        elif feature_selector == "ICA":
            loadings = ICA_dimension_reduction(x, num_comp, random_state)

        elif feature_selector == "FA":
            loadings = FA_dimension_reduction(x, num_comp, random_state)

        elif feature_selector == "SVD":
            loadings = tsvd_dimension_reduction(x, num_comp, random_state)

        elif feature_selector == "NMF":
            loadings = NMF_dimension_reduction(x, num_comp, random_state)

        else:
            std_devs = np.std(data, axis=0)
            loadings_df = pd.DataFrame(std_devs, columns=['Component1'], index=column_names)
            loadings_df = loadings_df.sort_values(ascending=False, by="Component1")
            num_comp_selector = gene_sig_exp.number_input(":blue[Number of ] :red[Most Variable Genes]",
                                                          min_value=10, max_value=max_genes,
                                                          value=def_comp_num, key='num_com_ni_std')
            selected_genes_1 = loadings_df.head(num_comp_selector)

        if feature_selector not in ["Standard Deviation"]:
            # Convert loadings to a DataFrame for better visualization
            loadings_df = pd.DataFrame(loadings, columns=[f'Components{i + 1}' for i in range(loadings.shape[1])],
                                       index=column_names)

            # Sum of squared loadings across principal components
            sum_squared_loadings = (loadings_df ** 2).sum(axis=1)

            # Average absolute loadings across principal components
            average_absolute_loadings = loadings_df.abs().mean(axis=1)

            # Max absolute loadings across principal components
            max_absolute_loading = loadings_df.abs().max(axis=1)

            # Variance loadings across principal components
            variance_of_loadings = loadings_df.var(axis=1)

            if num_comp == 1:
                ranking_methods = ["Max Absolute"]
                rank_meth = 0
            else:
                ranking_methods = ["Sum Square", "Average Absolute", "Max Absolute", "Variance"]

            # Copy one of this: variance_of_loadings, max_absolute_loading, average_absolute_loadings,
            # sum_squared_loadings
            ranking_method = gene_sig_exp.selectbox(":blue[Choose a ] :red[Ranking Method]", key="rank_method2",
                                                    help=dsc.describe_rank_method2(),  options=ranking_methods,
                                                    index=rank_meth)

            rank_info = ranking_method
            if ranking_method == "Sum Square":
                ranking_method = "sum_squared_loadings"
            elif ranking_method == "Average Absolute":
                ranking_method = "average_absolute_loadings"
            elif ranking_method == "Max Absolute":
                ranking_method = "max_absolute_loading"
            elif ranking_method == "Variance":
                ranking_method = "variance_of_loadings"

            ranking_method = eval(ranking_method)

            # Ranking genes
            ranking = ranking_method.sort_values(ascending=False)

            num_comp_selector = gene_sig_exp.number_input(":blue[Number of ] :red[Chosen Genes]", min_value=int(1),
                                                          max_value=max_genes, value=def_comp_num, key='num_com_ni2',
                                                          help=dsc.describe_num_com_ni2())

            # Select the top x genes
            selected_genes_1 = ranking.head(num_comp_selector)

            var_chb = gene_sig_exp.toggle(':blue[Use ] :red[Variance Ranking]', key='var_chb2', value=True,
                                            help=dsc.describe_var_chb2())
            gene_sig_exp.warning(
                ":blue[If you choose ] :red[Variance Ranking] :blue[, the variance calculation is based on the ] "
                " :red[original database].")

            if var_chb:
                selected_genes = df[df.index.isin(selected_genes_1.index)]
                selected_genes.dropna(how='all', axis=1, inplace=True)

                # Get the variance for all genes
                var_genes = selected_genes.var(axis=1)

                # Choose a percentile for variance threshold
                percentile_selection = gene_sig_exp.number_input(":blue[Choose ] :red[Variance Percentile]",
                                                                 min_value=0.01, max_value=99.9, value=var_value,
                                                                 key="var_percentile2", step=0.01,
                                                                 help=dsc.describe_var_percentile2())

                # Calculate the variance threshold based on the chosen percentile
                var_threshold = np.percentile(var_genes, percentile_selection)

                # Filter the genes based on the variance threshold
                selected_genes = selected_genes[var_genes > var_threshold]

            else:
                selected_genes = df[df.index.isin(selected_genes_1.index)]
                selected_genes.dropna(how='all', axis=1, inplace=True)
        else:
            selected_genes = df[df.index.isin(selected_genes_1.index)]
            selected_genes.dropna(how='all', axis=1, inplace=True)

        gene_sig_exp.info(f":blue[You selected: ] :red[{str(len(selected_genes))} ] :blue[genes]")
        gene_sig_exp.info(
            ':blue[The expected number may vary as different databases do ] :red[NOT ] '
            ':blue[contain the same genes!]')
        if gene_list_opt:
            gene_sig_exp.info(
                ':blue[Since the genes in the dataset and the gene list are not identical, the system may ] :red[NOT ] '
                ':blue[return the expected number of genes.]')
        # Print out the selected genes
        name_list = pd.DataFrame(selected_genes.index)

        if len(name_list) != 0:
            st.subheader(":blue[The Chosen Genes]")

            # Determine the number of columns
            num_cols = min(max(math.ceil(len(name_list) / 10), 1), 4)
            gene_list_cols = st.columns(num_cols)

            # Determine the length of each part
            n = math.ceil(len(name_list) / num_cols)

            # Add each part to a column
            for i in range(num_cols):
                start = i * n
                end = (i + 1) * n if i < num_cols - 1 else None  # The last part should include all remaining items
                df = pd.DataFrame(name_list[start:end])
                gene_list_cols[i].dataframe(df, use_container_width=True, hide_index=True)

        gene_list_file = pd.DataFrame(selected_genes.index)
        gene_list_file_dl = gene_list_file.to_csv(index=False, header=None)

        down_load_add_pdf = st.download_button(
            label="Download Chosen Gene List as CSV",
            key='download_add_btn2',
            data=gene_list_file_dl,
            file_name='gene_list.csv',
            mime='text/csv',
        )

        # Create a string with your desired format
        if var_chb:
            content = (f"Analysed dataset: {dataset_selection}\n"
                       f"Analysis method: {feature_selector}\nNumber of components: {num_comp}"
                       f"\nRandom State: {random_state}"
                       f"\nRanking Method: {rank_info}\nNumber of Chosen Genes: {num_comp_selector}"
                       f"\nUse Variance Filtering: {var_chb}\nVariance Percentile: {percentile_selection}"
                       f"\nFinal number of Genes: {len(selected_genes)}")
        else:
            content = (f"Analysed dataset: {dataset_selection}\n"
                       f"Analysis method: {feature_selector}\nNumber of components: {num_comp}"
                       f"\nRandom State: {random_state}"
                       f"\nRanking Method: {rank_info}\nNumber of Chosen Genes: {num_comp_selector}"
                       f"\nUse Variance Filtering: {var_chb}"
                       f"\nFinal number of Genes: {len(selected_genes)}")

        # Write the string to a .txt file
        with open(gene_sig_name, 'w') as file:
            file.write(content)

        return selected_genes

    else:
        gene_sig_exp.error(":blue[Please ] :red[upload ] :blue[a gene list]")


def clinical_gene_finder(data_set_option, df, gene_sig_name, user_id, clinic_st, clinic_df, subtype_exp,
                         gene_finder_cols):
    p_correction_kruskal = None
    kruskal_alpha = None
    sig_gene = None
    significant_results = None
    selected_genes = None
    random_state = None
    n_neighbors = None
    st_annotation = None
    column_name = None
    n_genes = None
    df_kw = None

    gene_sig_exp = st.expander(":blue[Advanced Settings]")
    selectbox_options = ["George-SCLC", "Anish-SCLC", 'Jiang-SCLC', "George-LCNEC", "Fernandez-Carcinoid",
                         "Alcala-Carcinoid", "Liu-SCLC", 'Rousseaux-Mixed(non-NE)']

    if data_set_option == "TCGA":
        dataset_index = 8
        selectbox_options = ["George-SCLC", "Anish-SCLC", 'Jiang-SCLC', "George-LCNEC", "Fernandez-Carcinoid",
                             "Alcala-Carcinoid", "Liu-SCLC", 'Rousseaux-Mixed(non-NE)', "TCGA"]

    elif data_set_option == "George-SCLC":
        dataset_index = 0
    elif data_set_option == "Anish-SCLC":
        dataset_index = 1
    elif data_set_option == "Jiang-SCLC":
        dataset_index = 2
    elif data_set_option == "George-LCNEC":
        dataset_index = 3
    elif data_set_option == "Fernandez-Carcinoid":
        dataset_index = 4
    elif data_set_option == "Alcala-Carcinoid":
        dataset_index = 5
    elif data_set_option == "Liu-SCLC":
        dataset_index = 6
    else:
        dataset_index = 7

    dataset_selection = gene_sig_exp.selectbox(
        ":blue[Choose a ] :red[Database ] :blue[to Serve as the Basis for the Analysis]",
        options=selectbox_options, index=dataset_index, key='ds_select2', help=dsc.describe_ds_select2())

    if dataset_selection == "George-SCLC":
        data, _, clinic_st = dl.george_sclc_data_load()
        data.index = data.index.str.upper()

    elif dataset_selection == "Anish-SCLC":
        data, _, clinic_st = dl.anish_sclc_data_load()
        data.index = data.index.str.upper()

    elif dataset_selection == "Jiang-SCLC":
        data, _, clinic_st = dl.jiang_sclc_data_load()
        data.index = data.index.str.upper()

    elif dataset_selection == "George-LCNEC":
        data, _, clinic_st = dl.george_lcnec_data_load()
        data.index = data.index.str.upper()

    elif dataset_selection == "Fernandez-Carcinoid":
        data, _, clinic_st = dl.carcinoid_data_load()
        data.index = data.index.str.upper()

    elif dataset_selection == "Liu-SCLC":
        data, _, clinic_st = dl.liu_data_load()
        data.index = data.index.str.upper()

    elif dataset_selection == "Alcala-Carcinoid":
        data, _, clinic_st = dl.alcala_data_load(False)
        data.index.rename("GeneSymbol", inplace=True)
        data.index = data.index.str.upper()

    elif dataset_selection == "Rousseaux-Mixed(non-NE)":
        data, _, clinic_st = dl.rousseaux_data_load(False, None)
        data.index = data.index.str.upper()

    elif dataset_selection == "TCGA":
        data = df

    else:
        st.info(":blue[Choose ] :red[ a Dataset]")
        data = None
        return None

    if dataset_selection == data_set_option:
        data = df.copy()

    selected_genes_anno = data.copy()
    selected_genes_anno = selected_genes_anno.transpose()
    if dataset_selection == 'Anish-SCLC':
        st_annotation = dl.data_load_anno_anish(selected_genes_anno)
        column_name = "NE"
    elif dataset_selection == 'George-SCLC':
        st_annotation = dl.data_load_anno_george(selected_genes_anno)
        column_name = "NE"
    elif dataset_selection == 'George-LCNEC':
        st_annotation = dl.data_load_anno_george_lcnec(selected_genes_anno)
        column_name = "Molecular_Subtypes"
    elif dataset_selection == 'Fernandez-Carcinoid':
        st_annotation = dl.data_load_anno_carcinoid(selected_genes_anno)
        column_name = "cluster_LNET"
    elif dataset_selection == 'Liu-SCLC':
        st_annotation = dl.data_load_anno_liu(selected_genes_anno)
        column_name = "NE"
    elif dataset_selection == 'Alcala-Carcinoid':
        st_annotation = dl.data_load_anno_alcala(selected_genes_anno)
        column_name = "Histopathology"
    elif dataset_selection == 'Rousseaux-Mixed(non-NE)':
        st_annotation = dl.data_load_anno_rousseaux(selected_genes_anno)
        column_name = "NE"
    elif dataset_selection == 'Jiang-SCLC':
        st_annotation = dl.data_load_anno_jiang(selected_genes_anno)
        column_name = "NE"


    if dataset_selection != "TCGA":

        data, merged_df = subtype_select_gf(data.copy().transpose(), st_annotation.copy(), data.copy().transpose(),
                                            dataset_selection, clinic_st.copy(), subtype_exp, key='subtype_ne')
    else:
        data, merged_df = subtypes_selector_tcga_gf(data.copy(), clinic_df.copy(), clinic_st.copy(), subtype_exp,
                                                    key='subtype_tcga')



    # Dropping duplicate columns from st_annotation
    merged_df = merged_df.loc[:, ~merged_df.columns.str.endswith('_drop')]

    # Toggle for gene list using
    gene_list_opt = gene_finder_cols[2].toggle(":blue[Use ] :red[Gene List]", key='gl_toggle2',
                                        help=dsc.describe_gl_toggle2())

    if gene_list_opt:
        gene_finder_cols[2].info(":blue[You can upload your file on the ] :red[sidebar!]")

    if gene_list_opt:
        gene_sig_exp.info(":blue[You can check your gene list in ] :red[Gene Handling Tab]")
        data = gl.app_maker(user_session_id=user_id, option_dataset=dataset_selection, df=data)

    if data is not None:
        # Separating out the features (assuming your data is just the expression values)
        data.dropna(axis=1, how='all', inplace=True)
        data = data.transpose()
        # Store the column names before converting data to numpy array
        column_names = data.columns
        data.dropna(axis=0, how='any', inplace=True)

    else:
        st.stop()

    # Dropping columns with only one unique value (including NaNs)
    merged_df = merged_df.loc[:, merged_df.nunique(dropna=False) > 1]

    # Dropping columns with all NaN values
    merged_df = merged_df.dropna(axis=1, how='all')

    if dataset_selection == "TCGA":
        column_name = merged_df.columns[-1]

    merged_df_index = merged_df.columns.get_loc(column_name)

    clinic_detail = gene_sig_exp.selectbox(":blue[Choose a ] :red[Clinical Detail ] :blue[for the Analysing]",
                                           options=merged_df.columns, index=merged_df_index, key="clin_detail",
                                           help=dsc. describe_clin_detail())

    merged_df = merged_df[[clinic_detail]]
    merged_df.dropna(how='any', axis=0, inplace=True)

    merged_df = merged_df[merged_df.index.isin(data.index)]
    data = data[data.index.isin(merged_df.index)]

    methods_selection = gene_sig_exp.selectbox(":blue[Select a ] :red[Method for Analysing]", key="meth_sele",
                                               options=["Kruskal Wallis", "Artificial Neural Network", "Nearest Centroid",
                                                        "Naive Bayes"],
                                               index=1, help=dsc.describe_method_sb())

    if methods_selection != "Kruskal Wallis":
        if data.shape[1] < 2500:
            max_genes = int(data.shape[1])
            def_comp_num = int(data.shape[1] / 2)
        else:
            max_genes = 2500
            def_comp_num = 1000
        n_genes = gene_sig_exp.number_input(":blue[Select the ] :red[Number of Genes]", min_value=5, max_value=max_genes,
                                            value=def_comp_num, key="n_genes_nir", help=dsc.describe_n_genes_nir())

    if methods_selection == "Kruskal Wallis":
        top_var = gene_sig_exp.number_input(
            ":blue[Select the ] :red[Top Variable Genes]",
            min_value=100,
            max_value=10000,
            value=5000,
            key="top_var",
            help=dsc.describe_num_com_ni_std()
        )
        gene_sig_exp.warning(
            ":blue[Kruskal-Wallis Test ] :red[Run Time Increases with the Number of Genes to be Tested.]"
        )

        df_kw = data.copy().transpose()
        # Calculate the standard deviation for each gene
        df_kw["STD"] = data.std(axis=0)

        # Sort the genes by standard deviation in descending order and assign back to df
        df_kw = df_kw.sort_values(ascending=False, by=["STD"])

        # Select the top 'top_var' genes
        df_kw = df_kw.head(top_var)

        # Drop the 'STD' column as it's no longer needed
        df_kw.drop(columns=['STD'], inplace=True)
        df_kw = df_kw.transpose()

    group_select = gene_sig_exp.multiselect(":blue[Select ] :red[Clinical Features]",
                                            options=merged_df[clinic_detail].unique(),
                                            key='group_ms', default=merged_df[clinic_detail].unique()[:2],
                                            max_selections=8, help=dsc.describe_group_ms())

    if methods_selection == "Kruskal Wallis":
        p_correction_kruskal = gene_sig_exp.selectbox(":blue[Select the ] :red[P Value Correction Method]",
                                                      options=["fdr_bh", "fdr_by", "bonferroni", "holm"], index=0,
                                                      key="p_correction_kruskal",
                                                      help=dsc.describe_p_correction_kruskal())
        kruskal_alpha = gene_sig_exp.select_slider(":blue[Select the ] :red[Alpha Value]", options=[0.001, 0.01, 0.05],
                                                   value=0.05, key="kruskal_alpha", help=dsc.describe_kruskal_alpha())

    random_state = gene_sig_exp.number_input(":blue[Set the ] :red[Random State ] :blue[of Your Work]",
                                             min_value=1, max_value=2 ** 32 - 1, value=42,
                                             key="random_state_gene_sign_clinic",
                                             help=dsc. describe_random_state_gene_sign2())

    # Check if each selected group has at least 3 occurrences
    valid_selection = True
    for feature in group_select:
        if merged_df[clinic_detail].value_counts().get(feature, 0) < 3:
            st.warning(f":red[The Feature ] '{feature}' :blue[Has Less Than 3 Occurrences.]")
            valid_selection = False
            st.stop()

    # Filter the dataframe if all selected features are valid
    if valid_selection and group_select:
        merged_df = merged_df[merged_df[clinic_detail].isin(group_select)]
    else:
        st.warning("Please select clinical features with at least 3 occurrences.")

    if merged_df[clinic_detail].nunique() < 2 or len(group_select) < 2:
        st.warning(":red[Not Enough Data ] :blue[for Analysis! Please Change the ] :red[Settings! ] :blue[You Can ] "
                   ":blue[Adjust the ] :red[Clinical Data, ] :blue[or Select Different ] :red[Groups.]")
        st.stop()
    else:
        balanced_toggle = st.toggle(
            ":red[Balance ] :blue[the Dataset]",
            help=":red[Enable this toggle to downsample unbalanced groups, ] "
                 ":blue[ensuring all groups have an equal number of samples.]"
        )

        if balanced_toggle:
            # Create an expander for advanced balanced settings
            balanced_expander = st.expander(
                ":blue[Advanced Setting of ] :red[Balanced Features]",
                expanded=False
            )
            # 1) Choose the number of top variable genes (slider or select_slider)
            top_genes_balanced = balanced_expander.select_slider(
                ":blue[Number of ] :red[Top Variable Genes]",
                options=[1000, 2000, 3000, 4000, 5000, 10000],
                value=5000,
                help="Select how many top variable genes to keep before downsampling."
            )

            # 2) Toggle whether to use PCA for dimension reduction
            use_pca_balanced = balanced_expander.toggle(
                ":red[Use PCA ] :blue[for Dimension Reduction]",
                value=True,
                help="Check to enable PCA on the selected genes."
            )

            # If PCA is enabled, let the user select the number of components
            if use_pca_balanced:
                n_components_balanced = balanced_expander.slider(
                    ":blue[Number of ] :red[PCA Components]",
                    min_value=1,
                    max_value=10,
                    value=2,
                    help="Select how many principal components to keep."
                )
            else:
                n_components_balanced = None  # No PCA if not checked

            # 3) Select the representative method (Centroid or Medoid)
            representative_method_balanced = balanced_expander.selectbox(
                ":red[Representative ] :blue[Method]",
                options=["Centroid", "Medoid"],
                index=0,
                help="Choose how to represent each group's center when downsampling."
            )

            # 4) Select the distance metric
            dist_methods_balanced = balanced_expander.selectbox(
                ":red[Distance ] :blue[Metric]",
                options=["Euclidean", "Pearson"],
                index=0,
                help="How to measure distance from each sample to the representative center."
            )

            var_reduction_toggle = balanced_expander.toggle(":blue[Reduce Number of ] :red[Genes by Variance]",
                                                          key="var_reduction_toggle", value=False)

            if var_reduction_toggle:
                gene_number_reduction = balanced_expander.number_input(":blue[Reeducated Number of ] :red[Genes]",
                                                                       key="gene_number_reduction", max_value=20000,
                                                                       min_value=2500, value=5000,
                                                                       help=":blue[Select how many genes to ] :red[reduce by variation.]")
            else:
                gene_number_reduction = None

            if methods_selection == "Kruskal Wallis":
                df_kw, merged_df = balanc_dataset(expression=df_kw,
                                                 merged_df=merged_df,
                                                 top_genes=top_genes_balanced,
                                                 use_pca=use_pca_balanced,
                                                 n_components=n_components_balanced,
                                                 representative_method=representative_method_balanced,  # or "Medoid"
                                                 distance_metric=dist_methods_balanced,       # or "Pearson"
                                                 random_state=random_state,
                                                 varience_reduction=var_reduction_toggle,
                                                 gene_reduction = gene_number_reduction
                 )
            else:
                data, merged_df = balanc_dataset(expression=data,
                                                merged_df=merged_df,
                                                top_genes=top_genes_balanced,
                                                use_pca=True,
                                                n_components=2,
                                                representative_method=representative_method_balanced,  # or "Medoid" "Centroid"
                                                distance_metric=dist_methods_balanced,  # or "Pearson" "Euclidean"
                                                random_state=random_state,
                                                varience_reduction=var_reduction_toggle,
                                                gene_reduction=gene_number_reduction
                 )

        if methods_selection == "Kruskal Wallis":
            significant_results, sig_gene = kruskal_wallis_test(df_kw, merged_df, p_correction_kruskal, kruskal_alpha)

            if len(sig_gene) == 0:
                return None
                #st.stop()
            else:
                if len(sig_gene) > 2500:
                    max_gene = 2500
                else:
                    max_gene = significant_results.shape[0] - 1
                n_genes = gene_sig_exp.number_input(
                    ":blue[Select the ] :red[Number of Genes]",
                    min_value=0,
                    max_value=max_gene,
                    value=int((significant_results.shape[0] - 1)/2),
                    key="n_genes_nir",
                    help=dsc.describe_n_genes_nir()
                )
                significant_results = significant_results.head(n_genes)
                sig_gene = sig_gene[:n_genes]
                sig_cols = gene_sig_exp.columns(3)
                sig_cols[1].write(significant_results)

        elif methods_selection == "Nearest Centroid":
            sig_gene = clanc(n_genes, merged_df, data, random_state)

        elif methods_selection == "Artificial Neural Network":
            sig_gene = mlp_scipy(n_genes, merged_df, data, random_state)

        elif methods_selection == "Naive Bayes":
            sig_gene = catBN(n_genes, merged_df, data, random_state)

        else:
            st.warning(":blue[Select a ] :red[Method for Analysing ] :blue[the Dataset!]")

        if len(sig_gene) == 0:
            st.warning(":red[Not Found Significant Genes! ] :blue[Change the Advanced Settings.]")
            return None
            #st.stop()
        gene_names = sig_gene.values
        selected_genes = df[df.index.isin(gene_names)]
        selected_genes.dropna(how='all', axis=1, inplace=True)
        gene_sig_exp.info(f":blue[You selected: ] :red[{str(len(selected_genes))} ] :blue[genes]")
        gene_sig_exp.info(
            ':blue[The expected number may vary as different databases do ] :red[NOT ] '
            ':blue[contain the same genes!]')
        gl.col_maker(selected_genes)

        gene_list_file = pd.DataFrame(selected_genes.index)
        gene_list_file_dl = gene_list_file.to_csv(index=False, header=None)

        down_load_add_pdf = st.download_button(
            label="Download Chosen Gene List as CSV",
            key='download_add_btn_clin',
            data=gene_list_file_dl,
            file_name='gene_list.csv',
            mime='text/csv',
        )

        return selected_genes
